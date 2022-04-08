#include "EMReader/emdata_.h"
#include "EMReader/util_func.h"
#include "EMReader/DataReader.h"
#include <cuda_runtime.h>
#include <cufft.h>
#include "GPU_func.cuh"
#include <malloc.h>
#include <cmath>
#include <iostream>

using namespace std;

void readEMData(Parameters *para, EulerData *euler)
{
    readEulerData(para->eulerf,euler);
	readSNRWeight(para->snr,&para->a,&para->b,&para->b2,&para->bfactor,&para->bfactor2,&para->bfactor3);
}

//add prefix of inlst to t. 
//.exp  inlst:../test.lsh  t:a.mrc => t:../a.mrc
void addPrefix(char *inlst, char *t)
{
    int l = strlen(inlst);
    char prefix[l+50];
    memset(prefix, 0, (l+50)*sizeof(char));
    int ed = -1;
    for(int i = l-1;i>=0;i--) if(inlst[i] == '/') {ed = i; break;}
    strncpy(prefix,inlst,ed+1);
    strcat(prefix,t);
    strcpy(t,prefix);
#ifdef DEBUG
    printf("--raw_image        %s\n",t);
#endif
}

void readRawImage(emdata *i2d_emdata,Parameters *para,int n, int *nn, char* t)
{
    vector<string> pairs;
    int return_code = readInLst_and_consturctPairs(para->inlst,t,&pairs,nn,n);
    if( return_code < 0) printf("Error %d occured in readRawImage!\n",return_code);

    //add prefix of inlst to t. 
    //.exp  inlst:../test.lsh  t:a.mrc => t:../a.mrc
    addPrefix(para->inlst,t);
    i2d_emdata->readImage(t,*nn);
    parsePairs(pairs, &para->defocus, &para->dfdiff, &para->dfang);
    
    para->dfu=para->defocus+para->dfdiff; //defocus is minus, so abs(dfu) < abs(dfv)
    para->dfv=para->defocus-para->dfdiff;
    para->lambda=12.2639/sqrt(para->energy*1000.0+0.97845*para->energy*para->energy);
    para->ds=1/(para->apix * para->padding_size);
}


//read template(float *)
//covert template from float* to cufftcomplex *
void readAndPaddingTemplate(Parameters *para,cufftComplex *h_templates,int N, double *sigmas)
{
    int padded_template_size = para->padding_size * para->padding_size;
    emdata *tp = new emdata();
    for(int J=0;J<N;J++)
    {
        tp->readImage(para->temp2d,J);
        float *data = tp->getData();
        if(para->padding_size < tp->header.nx || para->padding_size < tp->header.ny)
        {
            printf("Padded size is smaller than template.nx /ny\n");
            exit(-1) ;
        }
        int sx = (para->padding_size - tp->header.nx)/2;
        int sy = (para->padding_size - tp->header.ny)/2;
        for(int j=0;j<tp->header.ny;j++)
            for(int i=0;i<tp->header.nx;i++)
            {
                long long index = padded_template_size*J + (sy+j)*para->padding_size + (sx+i);
                h_templates[index].x = data[i+j*tp->header.nx];
            }
    }
	para->template_x = tp->header.nx;
    para->template_y = tp->header.ny;
    para->template_z = tp->header.nz;
    //para->overlap = para->template_x*0.13+1;
    //free heap memory
    if(tp!=NULL) delete tp;
    free(tp);
}

void cudaAllocTemplateMem(int N, int nx, int ny, Parameters *para,float **h_reduction_buf,float **d_reduction_buf,float **d_means,double **sigmas,double **d_sigmas,cufftComplex **h_templates,
    cufftComplex **d_templates,cufftComplex **CCG,cudaStream_t *stream,float **ra, float **rb, cufftHandle *plan_for_temp)
{
    int tmp = (para->padding_size - para->overlap);
    //num of blocks in x,y axis
    int block_x = (nx-para->overlap) / tmp;
    if((nx-para->overlap) % tmp > 0 ) block_x ++; 
    int block_y = (ny-para->overlap) / tmp;
    if((nx-para->overlap) % tmp > 0 ) block_y ++; 
    para->block_x = block_x;
    para->block_y = block_y;
    
    //N = max{num of tmplates ,  num of subimgs }
    N = max(N,block_x*block_y);

    //Number of Pixels for each padded template
    long long padded_template_size = para->padding_size*para->padding_size;

    //All padded templates (complex) At CPU,GPU
    *h_templates = (cufftComplex *)malloc(sizeof(cufftComplex)*padded_template_size*N);
    memset(*h_templates,0,sizeof(cufftComplex)*padded_template_size*N);
    CUDA_CALL(  cudaMalloc(d_templates,sizeof(cufftComplex)*padded_template_size*N)  );
    
    //Memory alloc for CCG
    CUDA_CALL(  cudaMalloc(CCG,sizeof(cufftComplex)*padded_template_size*N)  );

    //Store sigma value for all templates
    *sigmas = (double *)malloc(sizeof(double)*N);
    CUDA_CALL(  cudaMalloc(d_sigmas,sizeof(double)*N)  );
    
    //Cuda Stream
	cudaStreamCreate(stream);

    //Temp buffer for whiten
    CUDA_CALL(  cudaMalloc(ra,N*(RA_SIZE)*sizeof(float))  );
    CUDA_CALL(  cudaMalloc(rb,N*(RA_SIZE)*sizeof(float))  );

    //Buffer for reduction
    long long buf_size = max((long long)RA_SIZE*RA_SIZE/BLOCK_SIZE, 4*padded_template_size*N/BLOCK_SIZE);
    *h_reduction_buf = (float *)malloc(buf_size*sizeof(float));
    CUDA_CALL(  cudaMalloc(d_reduction_buf,buf_size*sizeof(float))  );

    //Store mean for every template
    CUDA_CALL(  cudaMalloc(d_means,sizeof(float)*N)  );

    /*
	cufftMakePlanMany(cufftHandle plan, int rank, int *n, int *inembed,
	int istride, int idist, int *onembed, int ostride,
	int odist, cufftType type, int batch, size_t *workSize);
	 */
	const int rank = 2;//维数
    int n[rank] = { para->padding_size, para->padding_size };//n*m
    int *inembed = n;//输入的数组size
    int istride = 1;//数组内数据连续，为1
    int idist = n[0] * n[1];//1个数组的内存大小
    int *onembed = n;//输出是一个数组的size
    int ostride = 1;//每点DFT后数据连续则为1
    int odist = n[0] * n[1];//输出第一个数组与第二个数组的距离，即两个数组的首元素的距离
    int batch = N;//批量处理的批数
    //采用cufftPlanMany方法
    
    //FFT handler for all templates
    CUFFT_CALL(cufftPlanMany(plan_for_temp, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, batch));//针对多信号同时进行FFT
    //Binding stream to plan
    CUFFT_CALL(cufftSetStream(*plan_for_temp, *stream));

}

void edgeNormalize(int N,double *sigmas,double *d_sigmas,cufftComplex *mask,cufftComplex *h_templates,cufftComplex *d_templates,float *h_buf,float *d_buf,
    float *N_buffer,Parameters para,cudaStream_t *stream)
{
    int l = para.padding_size;
    long long padded_template_size = l*l;
    int block_num = padded_template_size*N/BLOCK_SIZE;

    float r = l/(float)2-2;
    float up_bound = 0, low_bound = 0;
    if(r > 1)
    {
        up_bound = (r+1)*(r+1);
        low_bound = (r-1)*(r-1);
    }

    //Clear to {0}
    clear_image<<<block_num,BLOCK_SIZE,0,*stream>>>(mask);

    //em = cirMean()
    //generate Mask, count numbor of no-zero digits
    generate_mask<<<block_num,BLOCK_SIZE,BLOCK_SIZE*sizeof(float),*stream>>>(l,mask,r,d_buf,up_bound,low_bound);
    CUDA_CHECK();

    float N_buffer_host[2*N];
    CUDA_CALL(cudaMemcpyAsync(h_buf, d_buf,sizeof(float)*padded_template_size*N/BLOCK_SIZE, cudaMemcpyDeviceToHost, *stream));
    memset(N_buffer_host,0,N*sizeof(float));
    for(int k=0;k<padded_template_size*N/BLOCK_SIZE;k++)
    {
        int id = k/(padded_template_size/BLOCK_SIZE);
        N_buffer_host[id] += h_buf[k];
        //Number for every no-zero digits
    }

    //Calculate dot of mask and all templates
    CUDA_CALL(  cudaMemcpyAsync(N_buffer, N_buffer_host, sizeof(float)*N, cudaMemcpyHostToDevice, *stream)  );
    multiCount_dot<<<block_num,BLOCK_SIZE,BLOCK_SIZE*sizeof(float),*stream>>>(l,mask,d_templates,N_buffer,d_buf);
    CUDA_CHECK();
    CUDA_CALL(cudaMemcpyAsync(h_buf, d_buf,sizeof(float)*padded_template_size*N/BLOCK_SIZE, cudaMemcpyDeviceToHost, *stream));
    memset(N_buffer_host,0,N*sizeof(float));
    for(int k=0;k<padded_template_size*N/BLOCK_SIZE;k++)
    {
        int id = k/(padded_template_size/BLOCK_SIZE);
        N_buffer_host[id] += h_buf[k];
        //Dot result:em
    }
    
    CUDA_CALL(  cudaMemcpyAsync(d_templates, h_templates, sizeof(cufftComplex)*padded_template_size*N, cudaMemcpyHostToDevice, *stream)  );
    UpdateSigma<<<block_num,BLOCK_SIZE,BLOCK_SIZE*sizeof(float)*2,*stream>>>(d_templates,d_buf);
    CUDA_CHECK();
    CUDA_CALL(  cudaMemcpyAsync(h_buf, d_buf,2*sizeof(float)*padded_template_size*N/BLOCK_SIZE, cudaMemcpyDeviceToHost, *stream));
    
    //put em on GPU
    CUDA_CALL(  cudaMemcpyAsync(N_buffer, N_buffer_host, sizeof(float)*N, cudaMemcpyHostToDevice, *stream)  );

    memset(N_buffer_host,0,2*N*sizeof(float));
    for(int k=0;k<padded_template_size*N/BLOCK_SIZE;k++)
    {
        int id = k/(padded_template_size/BLOCK_SIZE);
        N_buffer_host[2*id] += h_buf[2*k];
        //sum of value
        N_buffer_host[2*id+1] += h_buf[2*k+1];
        //sum of value^2
    }

    for(int i=0;i<N;i++) {
        double mean = N_buffer_host[2*i] / (double)(l*l);
        sigmas[i] = sqrt(N_buffer_host[2*i+1] / (double)(l*l) - mean*mean);
        if(sigmas[i]<0 || !finite(sigmas[i])) sigmas[i] = 0 ;
    }

    //data[i]=(data[i]-em)/s;
    CUDA_CALL(  cudaMemcpyAsync(d_sigmas, sigmas, sizeof(double)*N, cudaMemcpyHostToDevice, *stream)  );
    scale_each<<<block_num,BLOCK_SIZE,0,*stream>>>(l,d_templates,N_buffer,d_sigmas);
    CUDA_CHECK();
}


void handleTemplate(int N, float *ra, float *rb,float *h_buf,float *d_buf,float *d_means,
    cufftComplex *h_templates,cufftComplex *d_templates,Parameters *para,cudaStream_t *stream, cufftHandle *plan_for_temp)
{

//***************************************************************
// 1. apply whitening filter
// 2. apply mask
// 3. low pass filter and apply weighting function
// 4. normalize
//***************************************************************

    long long padded_template_size = para->padding_size*para->padding_size;
    int block_num = padded_template_size*N/BLOCK_SIZE;

// **************************************************************
// apply whitening filter and do ift
// input: Padded IMAGE (Real SPACE)
// output: IMAGE_whiten (Fourier SPACE in RI)
// **************************************************************
//    CUDA_CALL(    ); => CUDA Exeception handler

    // Inplace FFT
    CUFFT_CALL(  cufftExecC2C(*plan_for_temp, d_templates, d_templates, CUFFT_FORWARD)  );
    // CUFFT will enlarge VALUE to N times. Restore it
    scale<<<block_num,BLOCK_SIZE,0,*stream>>>(d_templates,padded_template_size);
    CUDA_CHECK();

    //Whiten at fourier space
    //contain ri2ap
    SQRSum_by_circle<<<block_num,BLOCK_SIZE,0,*stream>>>(d_templates,ra,rb,para->padding_size,para->padding_size);
    CUDA_CHECK();

    //contain ap2ri
    whiten_Tmp<<<block_num,BLOCK_SIZE,0,*stream>>>(d_templates,ra,rb,para->padding_size);
    CUDA_CHECK();

// **************************************************************
// apply mask 
// input: whiten_IMAGE (Fourier SPACE in RI)
// output: masked_whiten_IMAGE (Fourier SPACE in RI)
// **************************************************************
    CUFFT_CALL(  cufftExecC2C(*plan_for_temp, d_templates, d_templates, CUFFT_INVERSE)  );
    apply_mask<<<block_num,BLOCK_SIZE,0,*stream>>>(d_templates,para->d_m,para->edge_half_width,para->padding_size);
    CUDA_CHECK();
    CUFFT_CALL(  cufftExecC2C(*plan_for_temp, d_templates, d_templates, CUFFT_FORWARD)  );
    // CUFFT will enlarge VALUE to N times. Restore it
    scale<<<block_num,BLOCK_SIZE,0,*stream>>>(d_templates,padded_template_size);
    CUDA_CHECK();

// **************************************************************
// 1. lowpass
// 2. apply weighting function
// 3. normlize
// input: masked_whiten_IMAGE (Fourier SPACE in RI) 
// output: PROCESSED_IMAGE (Fourier SPACE in AP)
// **************************************************************
    //contain ri2ap
    apply_weighting_function<<<block_num,BLOCK_SIZE,0,*stream>>>(d_templates,*para);
    CUDA_CHECK();

    compute_area_sum_ofSQR<<<block_num,BLOCK_SIZE,2*BLOCK_SIZE*sizeof(float),*stream>>>(d_templates,d_buf,para->padding_size,para->padding_size);
    CUDA_CHECK();
    CUDA_CALL(cudaMemcpyAsync(h_buf, d_buf,2*sizeof(float)*padded_template_size*N/BLOCK_SIZE, cudaMemcpyDeviceToHost, *stream));
    cudaStreamSynchronize(*stream);
    //After Reduction -> compute mean for each image
    float infile_mean[N],counts[N];
    memset(infile_mean,0,N*sizeof(float));
    memset(counts,0,N*sizeof(float));
    for(int k=0;k<padded_template_size*N/BLOCK_SIZE;k++)
    {
        int id = k/(padded_template_size/BLOCK_SIZE);
        infile_mean[id] += h_buf[2*k];
        counts[id] += h_buf[2*k+1];
    }
    for(int k=0;k<N;k++) infile_mean[k] = sqrtf(infile_mean[k]/(counts[k]*counts[k]));
    //Do Normalization with computed infile_mean[]
    CUDA_CALL(  cudaMemcpyAsync(d_means, infile_mean, sizeof(float)*N, cudaMemcpyHostToDevice, *stream)  );
    //Contain ap2ri
    normalize<<<block_num,BLOCK_SIZE,0,*stream>>>(d_templates,para->padding_size,para->padding_size,d_means);
    CUDA_CHECK();

}

void cudaAllocImageMem(float **d_image,cufftComplex **d_rotated_image,cufftComplex **rotated_splitted_image,cudaStream_t *stream,
    cufftHandle *plan_for_image,cufftHandle *plan_for_whole_IMG,int nx,int ny,int N,Parameters *para)
{
    //X Y Size for padded IMG
    int ix = para->block_x*para->padding_size;
    int iy = para->block_y*para->padding_size;

    CUDA_CALL(  cudaMalloc(d_image,nx*ny*sizeof(float))  );
    CUDA_CALL(  cudaMalloc(d_rotated_image,ix*iy*sizeof(cufftComplex))  );
    CUDA_CALL(  cudaMalloc(rotated_splitted_image,ix*iy*sizeof(cufftComplex))  );

    /*
	cufftMakePlanMany(cufftHandle plan, int rank, int *n, int *inembed,
	int istride, int idist, int *onembed, int ostride,
	int odist, cufftType type, int batch, size_t *workSize);
	 */
	const int rank = 2;//维数
    int n[rank] = { para->padding_size, para->padding_size };//n*m
    int *inembed = n;//输入的数组size
    int istride = 1;//数组内数据连续，为1
    int idist = n[0] * n[1];//1个数组的内存大小
    int *onembed = n;//输出是一个数组的size
    int ostride = 1;//每点DFT后数据连续则为1
    int odist = n[0] * n[1];//输出第一个数组与第二个数组的距离，即两个数组的首元素的距离
    int batch = para->block_x*para->block_y;//批量处理的批数
    //采用cufftPlanMany方法
    
    //FFT handler for all sub images
    CUFFT_CALL(cufftPlanMany(plan_for_image, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, batch));//针对多信号同时进行FFT
    //Binding stream to plan
    CUFFT_CALL(cufftSetStream(*plan_for_image, *stream));

    int n2[rank] = { nx, ny };//n*m
    inembed = n2;//输入的数组size
    istride = 1;//数组内数据连续，为1
    idist = n2[0] * n2[1];//1个数组的内存大小
    onembed = n2;//输出是一个数组的size
    ostride = 1;//每点DFT后数据连续则为1
    odist = n2[0] * n2[1];//输出第一个数组与第二个数组的距离，即两个数组的首元素的距离
    batch = 1;//批量处理的批数
    //FFT handler for single whole images
    CUFFT_CALL(cufftPlanMany(plan_for_whole_IMG, rank, n2, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, batch));//针对多信号同时进行FFT
    //Binding stream to plan
    CUFFT_CALL(cufftSetStream(*plan_for_whole_IMG, *stream));
}  

void init_d_image(Parameters para,cufftComplex *filter,float *d_image, float*ra, float *rb, emdata *image, int nx, int ny, cudaStream_t *stream,cufftHandle *plan_for_whole_IMG)
{
    //Translate origin of Image to (0,0)
    image->rotate(0);
    //Put Image on GPU
    CUDA_CALL(  cudaMemcpyAsync(d_image, image->getData(), sizeof(float)*nx*ny, cudaMemcpyHostToDevice, *stream)  );
    if(para.phase_flip == 1)
    {
        int block_num = nx*ny/BLOCK_SIZE+1;
        float2Complex<<<block_num,BLOCK_SIZE,0,*stream>>>(filter,d_image,nx,ny);
        //fft inplace
        CUFFT_CALL(cufftExecC2C(*plan_for_whole_IMG, filter, filter, CUFFT_FORWARD));
        scale<<<block_num,BLOCK_SIZE,0,*stream>>>(filter,nx*ny);
        CUDA_CHECK();

        //phase flip
        do_phase_flip<<<block_num,BLOCK_SIZE,0,*stream>>>(filter,para,nx,ny);
        CUDA_CHECK();

        //Whiten at fourier space
        clear_float<<<RA_SIZE/BLOCK_SIZE+1,BLOCK_SIZE,0,*stream>>>(ra);
        clear_float<<<RA_SIZE/BLOCK_SIZE+1,BLOCK_SIZE,0,*stream>>>(rb);
        //contain ri2ap
        SQRSum_by_circle<<<block_num,BLOCK_SIZE,0,*stream>>>(filter,ra,rb,nx,ny,1);
        CUDA_CHECK();

        // 1. whiten
        // 2. low pass
        // 3. weight
        // 4. ap2ri
        whiten_filetr_weight_Img<<<block_num,BLOCK_SIZE,0,*stream>>>(filter,ra,rb,nx,ny,para);
        CUDA_CHECK();

        //ifft inplace
        CUFFT_CALL(cufftExecC2C(*plan_for_whole_IMG, filter, filter, CUFFT_INVERSE));
        Complex2float<<<block_num,BLOCK_SIZE,0,*stream>>>(d_image,filter,nx,ny);

    }

}

void split_normalize_image(float *d_image,cufftComplex *d_rotated_image,float *h_buf,float *d_buf,float *d_means, Parameters para, cudaStream_t *stream, int nx, int ny, cufftHandle *image_plan)
{
    int l = para.padding_size;
    // Init d_rotated_imge to all {0}
    int ix = para.block_x;
    int iy = para.block_y;
    int blockIMG_num = ix*iy*l*l / BLOCK_SIZE;
    clear_image<<<blockIMG_num,BLOCK_SIZE,0,*stream>>>(d_rotated_image);
    CUDA_CHECK();
    // split Image into blocks with overlap  
    split_IMG<<<blockIMG_num,BLOCK_SIZE,0,*stream>>>(d_image,d_rotated_image,nx,ny,para.padding_size,para.block_x,para.overlap);	
    CUDA_CHECK();

    // do normalize to all subIMGs
    if(para.phase_flip == 1)
    {
        //Inplace FFT
        CUFFT_CALL(  cufftExecC2C(*image_plan, d_rotated_image, d_rotated_image, CUFFT_FORWARD)  );
        //Scale IMG to normel size
        scale<<<blockIMG_num,BLOCK_SIZE,0,*stream>>>(d_rotated_image,l*l);
        CUDA_CHECK();

        ri2ap<<<blockIMG_num,BLOCK_SIZE,0,*stream>>>(d_rotated_image);

        compute_area_sum_ofSQR<<<blockIMG_num,BLOCK_SIZE,2*BLOCK_SIZE*sizeof(float),*stream>>>(d_rotated_image,d_buf,l,l);
        CUDA_CHECK();
        CUDA_CALL(cudaMemcpyAsync(h_buf, d_buf,2*sizeof(float)*blockIMG_num, cudaMemcpyDeviceToHost, *stream));
        cudaStreamSynchronize(*stream);

        int N_IMG = para.block_x*para.block_y;
        //After Reduction -> compute mean for each image
        float infile_mean[N_IMG],counts[N_IMG];
        memset(infile_mean,0,N_IMG*sizeof(float));
        memset(counts,0,N_IMG*sizeof(float));
        for(int k=0;k<blockIMG_num;k++)
        {
            int id = k/( (l*l)/BLOCK_SIZE );
            infile_mean[id] += h_buf[2*k];
            counts[id] += h_buf[2*k+1];
        }
        for(int k=0;k<N_IMG;k++) infile_mean[k] = sqrtf(infile_mean[k]/(counts[k]*counts[k]));
        
        //Do Normalization with computed infile_mean[]
        CUDA_CALL(  cudaMemcpyAsync(d_means, infile_mean, sizeof(float)*N_IMG, cudaMemcpyHostToDevice, *stream)  );
        //Contain ap2ri
        normalize<<<blockIMG_num,BLOCK_SIZE,0,*stream>>>(d_rotated_image,l,l,d_means);
        CUDA_CHECK();

        //Inplace IFT
        CUFFT_CALL(  cufftExecC2C(*image_plan, d_rotated_image, d_rotated_image, CUFFT_INVERSE)  );
    }
    
}

void rotateImage(cufftComplex *splitted_image,cufftComplex *rotated_splitted_image,Parameters para, float e, cudaStream_t *stream)
{
    // Init d_rotated_imge to all {0}
    int ix = para.block_x;
    int iy = para.block_y;
    int blockIMG_num = ix*iy*para.padding_size*para.padding_size / BLOCK_SIZE;
    // rotate subIMG with angle "e"
    rotate_subIMG<<<blockIMG_num,BLOCK_SIZE,0,*stream>>>(splitted_image,rotated_splitted_image,e,para.padding_size);
    CUDA_CHECK();

}

void pickPartcles(cufftComplex *CCG,cufftComplex *d_templates,cufftComplex *rotated_splitted_image,float *h_buf,float *d_buf, Parameters para, 
    cufftHandle *template_plan, cufftHandle *image_plan,cudaStream_t *stream,int N,float *scores, int nx, int ny, float euler3)
{           
    int l = para.padding_size;
    long long padded_template_size = l*l;
    int blockGPU_num = padded_template_size*N/BLOCK_SIZE;
    int blockIMG_num = padded_template_size*para.block_x*para.block_y/BLOCK_SIZE;

    //peak,sum of data[i],sum of data[i]^2
    float peaks[N],pos[N],sums[N],sum2s[N];

    //find MAX score need initialize
    memset(scores,0,sizeof(float)*3*N);

    //Inplace FFT
    CUFFT_CALL(  cufftExecC2C(*image_plan, rotated_splitted_image, rotated_splitted_image, CUFFT_FORWARD)  );

    //Scale IMG to normel size
    scale<<<blockIMG_num,BLOCK_SIZE,0,*stream>>>(rotated_splitted_image,padded_template_size);
    CUDA_CHECK();

    //compute score for each block
    for(int j=0;j<para.block_y;j++)
    {
        for(int i=0;i<para.block_x;i++)
        {
            //find peak need initialize
            memset(peaks,0,sizeof(float)*N);
            memset(pos,0,sizeof(float)*N);
            memset(sums,0,sizeof(float)*N);
            memset(sum2s,0,sizeof(float)*N);
            //compute CCG
            compute_corner_CCG<<<blockGPU_num,BLOCK_SIZE,0,*stream>>>(CCG,d_templates,rotated_splitted_image,l,i+j*para.block_x);
            CUDA_CHECK();
            //Inplace IFT
            CUFFT_CALL(  cufftExecC2C(*template_plan, CCG, CCG, CUFFT_INVERSE)  );
            //Ingore padded 0 at raw IMG
            int x_bound = nx - i*(l-para.overlap);
            int y_bound = ny - j*(l-para.overlap);
            //find peak(position) and get sum of data,data^2
            get_peak_and_SUM<<<blockGPU_num,BLOCK_SIZE,4*BLOCK_SIZE*sizeof(float),*stream>>>(CCG,d_buf,l,para.d_m,x_bound,y_bound);
            CUDA_CHECK();
            CUDA_CALL(cudaMemcpyAsync(h_buf, d_buf, 4*sizeof(float)*padded_template_size*N/BLOCK_SIZE, cudaMemcpyDeviceToHost, *stream));
            cudaStreamSynchronize(*stream);

            //After Reduction -> compute mean for each image
            for(int k=0;k<(padded_template_size*N)/BLOCK_SIZE;k++)
            {
                int id = k/(padded_template_size/BLOCK_SIZE);
                if(peaks[id] < h_buf[4*k])
                {
                    peaks[id] = h_buf[4*k];
                    pos[id] = h_buf[4*k+1];
                }
                sums[id] += h_buf[4*k+2];
                sum2s[id] += h_buf[4*k+3];
            }

            //Update global score with local-block score for each template
            for(int J=0;J<N;J++) 
            {
                float ra = sums[J]-peaks[J];
                float rb = sum2s[J]-peaks[J]*peaks[J];
                float rc = padded_template_size - 1;
                float sd = sqrt(rb/rc - (ra/rc)*(ra/rc));
                float score;
                if(sd == 0) score = 0;
                else score = peaks[J]/sqrt(rb/rc - (ra/rc)*(ra/rc));
                
                int cx = (int)pos[J]%l;
                int cy = (int)pos[J]/l;

                //Rotate (cx,cy) to its soriginal angle
                float centerx =  i*(l-para.overlap) + (cx-l/2)*cos(euler3*PI/180)+(cy-l/2)*sin(euler3*PI/180)+l/2; // centerx
                float centery =  j*(l-para.overlap) + (cy-l/2)*cos(euler3*PI/180)-(cx-l/2)*sin(euler3*PI/180)+l/2; // centery
                //if(J==0) printf("%d %d %f %f %f %d %d =>(%d %d)\n",i,j,centerx,centery,score,cx,cy,i*(l-para.overlap),j*(l-para.overlap));

                if(scores[3*J] < score)
                {                    
                    //float Ny = para.d_m;
                    //if(cy-Ny/3>=0 && cx-Ny/3>=0 && cy+Ny/3<=ny && cx+Ny/3<=nx)
                    if(centerx>=0 && centerx<nx &&centery>=0 &&centery<ny)
                    {
                        scores[3*J] = score;
                        scores[3*J+1] = centerx;
                        scores[3*J+2] = centery;
                        //if(centerx <Ny/3 || centery<Ny/3 || centerx>(ny-Ny/3) || centery>(ny-Ny/3)) scores[3*k]=0;
                    }
                }

            }
        }
    }
}

void writeScoreToDisk(int N_tmp,float *scores,Parameters para,EulerData euler,FILE *fp, int nn, char *t, int nx, int ny, float euler3)
{
    for(int J=0;J<N_tmp;J++)
    {
        float score = scores[3*J];
        float centerx = scores[3*J+1];
        float centery = scores[3*J+2];

        if(score > para.thres)
        {
            fprintf(fp, "%d\t%s\tdefocus=%f\tdfdiff=%f\tdfang=%f\teuler=%f,%f,%f\tcenter=%f,%f\tscore=%f\n",
                    nn,t,(-1)*para.defocus,para.dfdiff,para.dfang,euler.euler1[J],euler.euler2[J],euler3,centerx,centery,score
            );
        }
    }
}

int main(int argc, char *argv[])
{
    //Timer
    time_t first, second;  
    first=time(NULL);

    //euler.length, number of template
    int N_tmp;

    //Print Help Message if no para input
    if(argc==1){
		printHelpMsg(); 
		return 0;
	}

    //Store all parameters(some from input, some from computation)
    Parameters para;
    //Stroe Euler data
    EulerData euler;
    
    readParameters(argc,argv,&para);
#ifdef DEBUG
    para.printAllPara();
#endif
    readEMData(&para,&euler);
    N_tmp = euler.length; 

    //IMG id
    int nn=0;
    //IMG filename
    char t[MAXPATHLEN+1];
    //Used to write res
    FILE *fp=fopen(para.outlst,"wb");
    // Set GPU device ID
    CUDA_CALL(  cudaSetDevice(para.device_id)  );

//*************************************************
// Process All Templates
// 1. Alloc memory
// 2. Read Template
// 3. Preprocess Template
//*************************************************
    //Store templates at CPU/GPU
    cufftComplex *h_templates,*d_templates;
    //Store res = F(template)F'(Sliptted_padded_IMG)
    cufftComplex *CCG;
    //Temp buffer for whiten
    float *ra,*rb;
    //Buffer for reduction
    float *h_reduction_buf,*d_reduction_buf;
    //Store mean for every template at GPU
    float *d_means;
    //Cuda Stream
    cudaStream_t stream;
    //cufft handler for template FFT
    cufftHandle plan_for_temp;
    //Store sigma value for all templates
    double *sigmas;
    double *d_sigmas;

    // READ FIRST IMG TO GET nx & ny
    //Image to be searched 
    emdata *image = new emdata();
    //Read Image
    readRawImage(image,&para,para.first,&nn,t);
    //size of first IMG
    int nx = image->header.nx;
    int ny = image->header.ny;
#ifdef DEBUG
        printf("IMG size:           %d %d\n",nx,ny);
        printf("Number of template: %d\n",N_tmp);
#endif
    cudaAllocTemplateMem(N_tmp,nx,ny,&para,&h_reduction_buf,&d_reduction_buf,&d_means,&sigmas,&d_sigmas,
        &h_templates,&d_templates,&CCG,&stream,&ra,&rb,&plan_for_temp);
    //free heap memory
    if(image!=NULL) delete image;

    //Loop for all Images. (Last - First) normally is 1;
    for(int n=para.first;n<para.last;n++)
	{
        if(n != para.first)
        {
            //Image to be searched 
            image = new emdata();
            //Read Image
            readRawImage(image,&para,n,&nn,t);
        }
        //Reading Template
        readAndPaddingTemplate(&para,h_templates,N_tmp,sigmas);

        //Process edge normlization
        //To avoid alloc memory,use CCG as tempary buf
        edgeNormalize(N_tmp,sigmas,d_sigmas,CCG,h_templates,d_templates,h_reduction_buf,d_reduction_buf,d_means,para,&stream);

        //whiten, apply mask, appy weighting ..
        handleTemplate(N_tmp,ra,rb,h_reduction_buf,d_reduction_buf,d_means,
            h_templates,d_templates,&para,&stream,&plan_for_temp);

//*************************************************
// Process Image to be searched 
// 0. Read
// 1. Rotate
// 2. Split
// 3. doFFT
//*************************************************   
        float *d_image;
        cufftComplex *d_rotated_image;
        cufftComplex *rotated_splitted_image;
        //cufft handler for sub-IMGs and whole-IMG FFT
        cufftHandle plan_for_image,plan_for_whole_IMG;
        
        cudaAllocImageMem(&d_image,&d_rotated_image,&rotated_splitted_image,&stream,&plan_for_image,&plan_for_whole_IMG,nx,ny,N_tmp,&para);
        // 1.Put Image on GPU 2.phaseflip
        init_d_image(para,rotated_splitted_image,d_image,ra,rb,image,nx,ny,&stream,&plan_for_whole_IMG);
        // split Image into blocks with overlap
        split_normalize_image(d_image,d_rotated_image,h_reduction_buf,d_reduction_buf,d_means,para,&stream,nx,ny,&plan_for_image);

        //Scores : [1-N]->socre  [N+1-2N]->cx+cy*padding_size
        float *scores = new float[N_tmp*3];    
        for(float euler3=0.0;euler3<360.0;euler3+=para.phi_step)
        {	 
#ifdef DEBUG
            printf("Now euler3 => %f / 360.0\n",euler3);
#endif 
            // rotate Image  
            rotateImage(d_rotated_image,rotated_splitted_image,para,euler3,&stream);

            //***************************************************
            //Pick particles from IMG with template
            //1. Calculate ccg
            //2. Find peak 
            //3. Calculate variance
            //4. Output score
            //***************************************************

            pickPartcles(CCG,d_templates,rotated_splitted_image,h_reduction_buf,d_reduction_buf,para,&plan_for_temp,&plan_for_image,&stream,N_tmp,scores,nx,ny,euler3);
            cudaStreamSynchronize(stream);
            writeScoreToDisk(N_tmp,scores,para,euler,fp,nn,t,nx,ny,euler3);
        }
        cufftDestroy(plan_for_whole_IMG);
        cufftDestroy(plan_for_image);
        cudaFree(rotated_splitted_image);
        cudaFree(d_rotated_image);
        cudaFree(d_image);
        delete []scores;
        
        //free heap memory
        if(image!=NULL) delete image;
    }

    cufftDestroy(plan_for_temp);
    cudaStreamDestroy(stream);
    cudaFree(d_reduction_buf);
    cudaFree(d_templates);
    cudaFree(d_sigmas);
    cudaFree(d_means);
    cudaFree(CCG);
    cudaFree(ra);
    cudaFree(rb);
    free(h_reduction_buf);
    free(h_templates);
    free(sigmas);
    fclose(fp);

    //Timer
    second=time(NULL);  
    printf("Total consumed time is: %f seconds\n",difftime(second,first)); 

    return 0;
}
