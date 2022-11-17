#include "EMReader/emdata_.h"
#include "EMReader/util_func.h"
#include "EMReader/DataReader.h"
#include <cuda_runtime.h>
#include <cufft.h>
#include "GPU_func.cuh"
#include <malloc.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sys/stat.h>

using namespace std;

//for debug
void out_mean_var_matrix(Parameters para,cufftComplex* ccg_sum, cudaStream_t *stream)
{
    int l = para.padding_size;   
    int data_size = l*l*para.block_x*para.block_y;
    cufftComplex *h_ccg_sum = (cufftComplex *)malloc(sizeof(cufftComplex)*data_size);
    CUDA_CALL( cudaMemcpyAsync(h_ccg_sum, ccg_sum, sizeof(cufftComplex)*data_size, cudaMemcpyDeviceToHost, *stream) );
    CUDA_CALL( cudaStreamSynchronize(*stream) );

    for(int i=0;i<para.block_x;i++)
    {
        for(int j=0;j<para.block_y;j++)
        {
            int k = i*para.block_y + j;
            if(k > 1) break; 
            string path = "debug/debug_mean_var/";
            string mean_path = path + "mean" + to_string(k) + ".txt";
            string var_path = path + "var" + to_string(k) + ".txt";

            ofstream meanFile;	
            meanFile.open(mean_path);
            ofstream varFile;	
            varFile.open(var_path);		
            for(int x = 0;x<l;x++)
            {
                for(int y=0;y<l;y++) 
                {
                    meanFile <<setprecision( 10 )<< h_ccg_sum[l*l*k + x*l+y].x/l/l <<" ";
                    varFile <<setprecision( 10 )<< h_ccg_sum[l*l*k + x*l+y].y/l/l <<" ";
                }
                meanFile<<endl;
                varFile<<endl;
            }            
            meanFile.close();	          
            varFile.close();	
        }
    }
    free(h_ccg_sum);
}

//for debug
void write_ccg(Parameters para,cufftComplex* d_ccg, int N, float e, cudaStream_t *stream)
{
    int l = para.padding_size;   
    int data_size = N*l*l;
    cufftComplex *h_ccg;
    h_ccg = (cufftComplex *)malloc(sizeof(cufftComplex)*data_size);

    CUDA_CALL( cudaMemcpyAsync(h_ccg, d_ccg, sizeof(cufftComplex)*data_size, cudaMemcpyDeviceToHost, *stream) );
    CUDA_CALL( cudaStreamSynchronize(*stream) );
    for(int i=0;i<N;i++)
    {
        string path = "debug/ccg_all_angle/ccg_"+ to_string(e) + "_" + to_string(i) + ".txt";
        ofstream outFile;	
        outFile.open(path);	
        for(int x = 0;x<l;x++)
        {
            for(int y=0;y<l;y++) 
            {
                outFile <<setprecision( 10 )<< h_ccg[l*l*i + x*l + y].x/(l*l) <<" ";
            }
            outFile<<endl;
        }                      
        outFile.close();	
    }
    free(h_ccg);
}

//for debug
void write_temp(Parameters para,cufftComplex* d_ccg, int N, float e, cudaStream_t *stream)
{
    int l = para.padding_size;   
    int data_size = N*l*l;
    cufftComplex *h_buf;
    h_buf = (cufftComplex *)malloc(sizeof(cufftComplex)*data_size);
    CUDA_CALL( cudaMemcpyAsync(h_buf, d_ccg, sizeof(cufftComplex)*data_size, cudaMemcpyDeviceToHost, *stream) );
    CUDA_CALL( cudaStreamSynchronize(*stream) );
    for(int i=0;i<1;i++)
    {
        string path = "debug/template_rotate/temp_"+ to_string(e) + "_" + to_string(i) + ".txt";
        ofstream outFile;	
        outFile.open(path);	
        for(int x = 0;x<l;x++)
        {
            for(int y=0;y<l;y++) 
            {
                outFile <<setprecision( 10 )<< h_buf[l*l*i + x*l + y].x <<" ";
            }
            outFile<<endl;
        }                      
        outFile.close();	
    }
    free(h_buf);
}


//for debug
//check GPU mem leaking
void checkGPUMem()
{
    size_t avail(0);
    size_t total(0);
    cudaMemGetInfo(&avail,&total);
    printf("Avail:%zu / Total:%zu \n",avail,total);    
}

//add prefix of inlst to t. 
//.exp  inlst:../test.lsh  t:a.mrc => t:../a.mrc
void addPrefix(char *inlst, char *t, int header)
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
    if(header == 0) printf("--raw_image        %s\n",t);
#endif
}

void readRawImage(emdata *i2d_emdata,Parameters *para,int n, int *nn, char* t, int header = 0)
{
    vector<string> pairs;
    int return_code = readInLst_and_consturctPairs(para->inlst,t,&pairs,nn,n);
    if( return_code < 0) printf("Error %d occured in readRawImage!\n",return_code);

    //add prefix of inlst to t. 
    //.exp  inlst:../test.lsh  t:a.mrc => t:../a.mrc
    addPrefix(para->inlst, t, header);
    i2d_emdata->readImage(t,*nn);
    parsePairs(pairs, &para->defocus, &para->dfdiff, &para->dfang);
    
    para->dfu=para->defocus+para->dfdiff; //defocus is minus, so abs(dfu) < abs(dfv)
    para->dfv=para->defocus-para->dfdiff;
    para->lambda=12.2639/sqrt(para->energy*1000.0+0.97845*para->energy*para->energy);
    para->ds=1/(para->apix * para->padding_size);
}

void adjustPaddingSize(Parameters *para)
{
    emdata *tp = new emdata();
    //read first temp
    tp->readImage(para->temp2d,0);
    
    #ifdef DEBUG
    printf("Temp size:         %d %d\n",tp->header.nx,tp->header.ny);
    #endif
    
    if(para->padding_size < tp->header.nx || para->padding_size < tp->header.ny)
    {
        printf("Padding size (%d) is smaller than template.nx/ny=(%d,%d)\n",para->padding_size, tp->header.nx ,tp->header.ny);
        para->padding_size = 0;
        while(para->padding_size < tp->header.nx || para->padding_size < tp->header.ny) para->padding_size += 32;
        printf("Aujust padding size to %d\n",para->padding_size);
    }

    //free heap memory
    if(tp!=NULL) delete tp;
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
       // double mean = 0, sigma = 0;
       /* if(para->norm_type == 1)
        {
            for(int j=0;j<tp->header.ny;j++)
                for(int i=0;i<tp->header.nx;i++)
                {
                    double cur = data[i+j*tp->header.nx];
                    mean = mean + cur;
                    sigma = sigma + (cur*cur);
                }
            mean = mean/(tp->header.ny*tp->header.nx);
            sigma = sqrtf(sigma/(tp->header.ny*tp->header.nx) - mean*mean);
        }*/
        


        int sx = (para->padding_size - tp->header.nx)/2;
        int sy = (para->padding_size - tp->header.ny)/2;
        for(int j=0;j<tp->header.ny;j++)
            for(int i=0;i<tp->header.nx;i++)
            {
                long long index = padded_template_size*J + (sy+j)*para->padding_size + (sx+i);
                float cur = data[i+j*tp->header.nx];
              //  if(para->norm_type == 1) cur = sigma>0? (cur-mean)/sigma:cur;
                h_templates[index].x = cur;
            }
    }
	para->template_x = tp->header.nx;
    para->template_y = tp->header.ny;
    para->template_z = tp->header.nz;
    //free heap memory
    if(tp!=NULL) delete tp;
}

void cudaAllocTemplateMem(int N, int nx, int ny,vector<int> &block_off_x, vector<int> & block_off_y, Parameters *para,float **h_reduction_buf,float **d_reduction_buf,float **d_means,double **sigmas,double **d_sigmas,cufftComplex **h_templates,
    cufftComplex **d_templates,cufftComplex **CCG,cufftComplex **CCG_sum,cufftComplex **CCG_buf,cudaStream_t *stream,float **ra, float **rb, cufftHandle *plan_for_temp)
{
    int interval = (para->padding_size - para->overlap);
    //num of blocks in x,y axis
    int off_x = 0, off_y = 0;
    int block_x = 1, block_y = 1;
    block_off_x.push_back(0);
    block_off_y.push_back(0);
    while(off_x + para->padding_size < nx) 
    {
        off_x += interval;
        if(off_x + para->padding_size >= nx)  off_x = nx - para->padding_size;
        block_off_x.push_back(off_x);
        block_x ++;
    }

    while(off_y + para->padding_size < ny) 
    {
        off_y += interval;
        if(off_y + para->padding_size >= ny)  off_y = ny - para->padding_size;
        block_off_y.push_back(off_y);
        block_y ++;
    }	
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
    //Memory of n&n2
    CUDA_CALL(  cudaMalloc(CCG_sum,sizeof(cufftComplex)*padded_template_size*block_x*block_y)  );
    CUDA_CALL(  cudaMalloc(CCG_buf,sizeof(cufftComplex)*padded_template_size*N)  );

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
    cufftComplex *d_templates,Parameters *para,cudaStream_t *stream, cufftHandle *plan_for_temp)
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
// 1. lowpass
// 2. apply weighting function
// 3. normlize
// input: masked_whiten_IMAGE (Fourier SPACE in RI) 
// output: PROCESSED_IMAGE (Fourier SPACE in AP)
// **************************************************************
    //contain ri2ap
    apply_weighting_function<<<block_num,BLOCK_SIZE,0,*stream>>>(d_templates,*para);
    CUDA_CHECK();
    //ap2ri<<<block_num,BLOCK_SIZE,0,*stream>>>(d_templates);
    CUDA_CHECK();
   // if(para->norm_type == 0)
   // {
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

        // **************************************************************
        // apply mask 
        // input: normed_IMAGE (Fourier SPACE in RI)
        // **************************************************************
        CUFFT_CALL(  cufftExecC2C(*plan_for_temp, d_templates, d_templates, CUFFT_INVERSE)  );
    //    apply_mask<<<block_num,BLOCK_SIZE,0,*stream>>>(d_templates,para->d_m,para->edge_half_width,para->padding_size);
    //    CUDA_CHECK();
   // }
   /*else if(para->norm_type == 1)
    {
        ap2ri<<<block_num,BLOCK_SIZE,0,*stream>>>(d_templates);
        CUFFT_CALL(  cufftExecC2C(*plan_for_temp, d_templates, d_templates, CUFFT_INVERSE)  );
        apply_mask<<<block_num,BLOCK_SIZE,0,*stream>>>(d_templates,para->d_m,para->edge_half_width,para->padding_size);
        CUDA_CHECK();
    }*/
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

    int n2[rank] = { ny, nx };//n*m
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

    float *h_img = image->getData();
    float sum = 0, sum_s2 = 0;
    for(int i=0;i<nx * ny;i++)
    {
        float cur = h_img[i];
        sum += cur /nx / ny;
        sum_s2 += (cur*cur /nx/ny);
    }
    float avg = sum ;
    float var = sqrt(sum_s2 - avg*avg);

    int up_bound = avg + 6*var;
    int low_bound = avg - 6*var;
    if(sum_s2 > avg*avg)
        for(int i=0;i<nx*ny;i++)
            if(h_img[i] > up_bound || h_img[i] < low_bound) h_img[i] = avg;
    if(para.invert == 1)
        for(int i=0;i<nx*ny;i++) h_img[i] = (-1) * h_img[i];
    //Put Image on GPU
    CUDA_CALL(  cudaMemcpyAsync(d_image, image->getData(), sizeof(float)*nx*ny, cudaMemcpyHostToDevice, *stream)  );

    
    int block_num = nx*ny/BLOCK_SIZE+1;
    float2Complex<<<block_num,BLOCK_SIZE,0,*stream>>>(filter,d_image,nx,ny);
    //fft inplace
    CUFFT_CALL(cufftExecC2C(*plan_for_whole_IMG, filter, filter, CUFFT_FORWARD));
    scale<<<block_num,BLOCK_SIZE,0,*stream>>>(filter,nx*ny);
    CUDA_CHECK();

    if(para.phase_flip == 1)
    {
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
    }
    
    // 0Hz -> 0
    set_0Hz_to_0_at_RI<<<block_num,BLOCK_SIZE,0,*stream>>>(filter,nx,ny);
    CUDA_CHECK();

    //ifft inplace
    CUFFT_CALL(cufftExecC2C(*plan_for_whole_IMG, filter, filter, CUFFT_INVERSE));
    Complex2float<<<block_num,BLOCK_SIZE,0,*stream>>>(d_image,filter,nx,ny);

}

void split_normalize_image(float *d_image,cufftComplex *d_rotated_image,float *h_buf,float *d_buf,float *d_means, Parameters para, cudaStream_t *stream, vector<int> &block_off_x, vector<int> &block_off_y, int nx, int ny, cufftHandle *image_plan)
{
    int l = para.padding_size;
    // Init d_rotated_imge to all {0}
    int ix = para.block_x;
    int iy = para.block_y;
    int blockIMG_num = ix*iy*l*l / BLOCK_SIZE;
    clear_image<<<blockIMG_num,BLOCK_SIZE,0,*stream>>>(d_rotated_image);
    CUDA_CHECK();

    // cpy data to GPU
    int h_off_x[block_off_x.size()], h_off_y[block_off_y.size()];
    for(int i = 0; i < block_off_x.size(); i ++) h_off_x[i] = block_off_x[i];
    for(int i = 0; i < block_off_y.size(); i ++) h_off_y[i] = block_off_y[i];
    int *d_off_x, *d_off_y;
    CUDA_CALL(  cudaMalloc(&d_off_x,sizeof(int)*block_off_x.size())  );
    CUDA_CALL(  cudaMalloc(&d_off_y,sizeof(int)*block_off_y.size())  );
    CUDA_CALL(  cudaMemcpyAsync(d_off_x, h_off_x,sizeof(int)*block_off_x.size(), cudaMemcpyHostToDevice, *stream)  );
    CUDA_CALL(  cudaMemcpyAsync(d_off_y, h_off_y,sizeof(int)*block_off_y.size(), cudaMemcpyHostToDevice, *stream)  );

    // split Image into blocks with overlap  
    split_IMG<<<blockIMG_num,BLOCK_SIZE,0,*stream>>>(d_image,d_rotated_image,d_off_x,d_off_y,nx,ny,para.padding_size,para.block_x,para.overlap);	
    CUDA_CHECK();

    // free buf
    cudaFree(d_off_x);
    cudaFree(d_off_y);

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

void rotateTemplate(cufftHandle *plan_for_temp, cufftComplex *d_templates,cufftComplex *rotated_templates,float *h_buf,float *d_buf,float *d_means,Parameters para, int N, float e, cudaStream_t *stream)
{
    // Init d_rotated_imge to all {0}
    int l = para.padding_size;
    long long padded_template_size = l*l;
    int blockGPU_num = padded_template_size*N/BLOCK_SIZE;
    // rotate subIMG with angle "e"
    rotate_subIMG<<<blockGPU_num,BLOCK_SIZE,0,*stream>>>(d_templates,rotated_templates,e,l);
    CUDA_CHECK();
    apply_mask<<<blockGPU_num,BLOCK_SIZE,0,*stream>>>(rotated_templates,para.d_m,para.edge_half_width,para.padding_size);
    CUDA_CHECK();
    //debug
   // write_temp(para,rotated_templates,N,e,stream);
    if(para.norm_type == 1)
	{
        compute_sum_sqr<<<blockGPU_num,BLOCK_SIZE,2*BLOCK_SIZE*sizeof(float),*stream>>>(rotated_templates,d_buf,para.padding_size,para.padding_size);
        CUDA_CHECK();
        CUDA_CALL(cudaMemcpyAsync(h_buf, d_buf,2*sizeof(float)*padded_template_size*N/BLOCK_SIZE, cudaMemcpyDeviceToHost, *stream));
        cudaStreamSynchronize(*stream);
        
        float infile_mean[N],infile_sqr[N];
        memset(infile_mean,0,N*sizeof(float));
        memset(infile_sqr,0,N*sizeof(float));
        for(int k=0;k<padded_template_size*N/BLOCK_SIZE;k++)
        {
            int id = k/(padded_template_size/BLOCK_SIZE);
            infile_mean[id] += h_buf[2*k];
            infile_sqr[id] += h_buf[2*k+1];
        }  
        for(int k=0;k<N;k++) 
        {
            infile_mean[k] = (infile_mean[k]/(padded_template_size));
            infile_sqr[k] = infile_sqr[k]/padded_template_size-infile_mean[k]*infile_mean[k];
        }
        //Do Normalization with computed infile_mean[]
        CUDA_CALL(  cudaMemcpyAsync(d_means, infile_mean, sizeof(float)*N, cudaMemcpyHostToDevice, *stream)  );
        substract_by_mean<<<blockGPU_num,BLOCK_SIZE,0,*stream>>>(rotated_templates,l,l,d_means);
        CUDA_CHECK(); 
        CUDA_CALL(  cudaMemcpyAsync(d_means, infile_sqr, sizeof(float)*N, cudaMemcpyHostToDevice, *stream)  );
        divided_by_var<<<blockGPU_num,BLOCK_SIZE,0,*stream>>>(rotated_templates,para.padding_size,para.padding_size,d_means);
        //normalize<<<block_num,BLOCK_SIZE,0,*stream>>>(rotated_templates,para.padding_size,para.padding_size,d_means);
        CUDA_CHECK();  
    }		
        //Inplace FFT
    CUFFT_CALL(  cufftExecC2C(*plan_for_temp, rotated_templates, rotated_templates, CUFFT_FORWARD)  );
    //Scale IMG to normal size
    scale<<<blockGPU_num,BLOCK_SIZE,0,*stream>>>(rotated_templates,padded_template_size);
    CUDA_CHECK();
}

void clear_ccg_sum(cufftComplex *CCG_sum, Parameters para, cudaStream_t *stream)
{
    int block_num = para.padding_size*para.padding_size*para.block_x*para.block_y / BLOCK_SIZE;
    clear_image<<<block_num,BLOCK_SIZE,0,*stream>>>(CCG_sum);
}

void compute_CCG_sum(cufftComplex *CCG,cufftComplex *CCG_sum,cufftComplex *d_templates,cufftComplex *rotated_splitted_image,
    Parameters para, cufftHandle *template_plan,cudaStream_t *stream,int N, float euler3, int N_euler)
{           
    int l = para.padding_size;
    long long padded_template_size = l*l;
    int blockGPU_num = padded_template_size*N/BLOCK_SIZE;

    //compute score for each block
    for(int j=0;j<para.block_y;j++)
    {
        for(int i=0;i<para.block_x;i++)
        {
            //compute CCG
            compute_corner_CCG<<<blockGPU_num,BLOCK_SIZE,0,*stream>>>(CCG,d_templates,rotated_splitted_image,l,i+j*para.block_x);
            CUDA_CHECK();
            //Inplace IFT
            CUFFT_CALL(  cufftExecC2C(*template_plan, CCG, CCG, CUFFT_INVERSE)  );
            //debug
            //if(i ==0 && j == 0) write_ccg(para,CCG,N,euler3,stream);
            //compute avg/vairance
            add_CCG_to_sum<<<l*l/BLOCK_SIZE,BLOCK_SIZE,0,*stream >>>(CCG_sum,CCG,l,N,i+j*para.block_x);
            CUDA_CHECK();
        }
    }
}

void compute_CCG_mean(cufftComplex *CCG_sum, Parameters para,int N_tmp, int N_euler, cudaStream_t *stream)
{
    int l = para.padding_size;
    long long padded_template_size = l*l;
    int blockIMG_num = padded_template_size*para.block_x*para.block_y/BLOCK_SIZE;
    set_CCG_mean<<<blockIMG_num,BLOCK_SIZE>>>(CCG_sum,para.padding_size,N_tmp,N_euler);
}

void fft_on_img(cufftComplex *img, cufftHandle *image_plan, Parameters para,cudaStream_t *stream)
{
    int l = para.padding_size;
    long long padded_template_size = l*l;
    int blockIMG_num = padded_template_size*para.block_x*para.block_y/BLOCK_SIZE;

    //Inplace FFT
    CUFFT_CALL(  cufftExecC2C(*image_plan, img, img, CUFFT_FORWARD)  );
    //Scale IMG to normel size
    scale<<<blockIMG_num,BLOCK_SIZE,0,*stream>>>(img,padded_template_size);
    CUDA_CHECK();
}

void pickPartcles_CCG_NORM(cufftComplex *CCG,cufftComplex *CCG_sum,cufftComplex *d_templates,cufftComplex *rotated_splitted_image,float *h_buf,float *d_buf, Parameters para, 
    cufftHandle *template_plan, cudaStream_t *stream,int N,vector<vector<float> > &score_info,vector<int> &block_off_x,vector<int> &block_off_y,
     int nx, int ny, float euler3)
{           
    int l = para.padding_size;
    long long padded_template_size = l*l;
    int blockGPU_num = padded_template_size*N/BLOCK_SIZE;

    //compute score for each block
    for(int j=0;j<para.block_y;j++)
    {
        for(int i=0;i<para.block_x;i++)
        {
            //compute CCG
            compute_corner_CCG<<<blockGPU_num,BLOCK_SIZE,0,*stream>>>(CCG,d_templates,rotated_splitted_image,l,i+j*para.block_x);
            CUDA_CHECK();
            //Inplace IFT
            CUFFT_CALL(  cufftExecC2C(*template_plan, CCG, CCG, CUFFT_INVERSE)  );
            //update CCG with avg/vairance
            update_CCG<<<blockGPU_num,BLOCK_SIZE,0,*stream >>>(CCG_sum,CCG,l,i+j*para.block_x);
            CUDA_CHECK();

            // find peak in each block
            get_peak_pos<<<blockGPU_num,BLOCK_SIZE,2*BLOCK_SIZE*sizeof(float),*stream>>>(CCG,d_buf,l,para.d_m);
            CUDA_CHECK();
            CUDA_CALL(cudaMemcpyAsync(h_buf, d_buf, 2*sizeof(float)*padded_template_size*N/BLOCK_SIZE, cudaMemcpyDeviceToHost, *stream));
            cudaStreamSynchronize(*stream);

            //After Reduction -> compute mean for each image
            for(int k=0;k<(padded_template_size*N)/BLOCK_SIZE;k++)
            {
                int J = k/(padded_template_size/BLOCK_SIZE);
                if(h_buf[2*k] >= para.thres )
                {
                    float score = h_buf[2*k];
                    int centerx = block_off_x[i] + (int)h_buf[2*k+1] % l;
                    int centery = block_off_y[j] + (int)h_buf[2*k+1] / l;
                    //printf("J:%d  k:%d local_k:%lld i:%d j:%d score:%f  pos:%f\n",J,k,k%(padded_template_size/BLOCK_SIZE),i,j,h_buf[2*k],h_buf[2*k+1]);
                    if(centerx>=para.d_m && centerx<nx-para.d_m &&centery>=para.d_m &&centery<ny-para.d_m)
                    {
                        score_info[0].push_back(score);
                        score_info[1].push_back(centerx);
                        score_info[2].push_back(centery);
                        score_info[3].push_back(J);
                    }
                }
            }
           
        }
    }
}


void pickPartcles_IMG_NORM(cufftComplex *CCG,cufftComplex *d_templates,cufftComplex *rotated_splitted_image,float *h_buf,float *d_buf, Parameters para, 
    cufftHandle *template_plan, cufftHandle *image_plan,cudaStream_t *stream,int N,vector<vector<float> > &score_info, vector<int> &block_off_x,vector<int> &block_off_y,int nx, int ny, float euler3)
{           
    int l = para.padding_size;
    long long padded_template_size = l*l;
    int blockGPU_num = padded_template_size*N/BLOCK_SIZE;

    //peak,sum of data[i],sum of data[i]^2
    float peaks[N],pos[N],sums[N],sum2s[N];

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
            //find peak(position) and get sum of data,data^2
            get_peak_and_SUM<<<blockGPU_num,BLOCK_SIZE,4*BLOCK_SIZE*sizeof(float),*stream>>>(CCG,d_buf,l,para.d_m);
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
                if(sd > 0) score = peaks[J]/sqrt(rb/rc - (ra/rc)*(ra/rc));
                
                float centerx =  block_off_x[i] + (int)pos[J]%l; 
                float centery =  block_off_y[j] + (int)pos[J]/l; 

                if(score >= para.thres)
                {                    
                    if(centerx>=para.d_m && centerx<nx-para.d_m &&centery>=para.d_m &&centery<ny-para.d_m)
                    {
                        score_info[0].push_back(score);
                        score_info[1].push_back(centerx);
                        score_info[2].push_back(centery);
                        score_info[3].push_back(J);
                    }
                }
            }
        }
    }
}

void pickPartcles_CCG_IMG_NORM(cufftComplex *CCG,cufftComplex *CCG_sum,cufftComplex *d_templates,cufftComplex *rotated_splitted_image,float *h_buf,float *d_buf, Parameters para, 
    cufftHandle *template_plan, cudaStream_t *stream,int N,vector<vector<float> > &score_info,vector<int> &block_off_x,vector<int> &block_off_y,
     int nx, int ny, float euler3)
{           
    int l = para.padding_size;
    long long padded_template_size = l*l;
    int blockGPU_num = padded_template_size*N/BLOCK_SIZE;

    //peak,sum of data[i],sum of data[i]^2
    float peaks[N],pos[N],sums[N],sum2s[N];

    //compute score for each block
    for(int j=0;j<para.block_y;j++)
    {
        for(int i=0;i<para.block_x;i++)
        {
            //compute CCG
            compute_corner_CCG<<<blockGPU_num,BLOCK_SIZE,0,*stream>>>(CCG,d_templates,rotated_splitted_image,l,i+j*para.block_x);
            CUDA_CHECK();
            //Inplace IFT
            CUFFT_CALL(  cufftExecC2C(*template_plan, CCG, CCG, CUFFT_INVERSE)  );
            //update CCG with avg/vairance
            update_CCG<<<blockGPU_num,BLOCK_SIZE,0,*stream >>>(CCG_sum,CCG,l,i+j*para.block_x);
            CUDA_CHECK();

            //find peak(position) and get sum of data,data^2
            get_peak_and_SUM<<<blockGPU_num,BLOCK_SIZE,4*BLOCK_SIZE*sizeof(float),*stream>>>(CCG,d_buf,l,para.d_m);
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
                float score = 0;
                if(sd > 0) score = peaks[J]/sqrt(rb/rc - (ra/rc)*(ra/rc));
                
                int cx = (int)pos[J]%l;
                int cy = (int)pos[J]/l;

                //Rotate (cx,cy) to its soriginal angle
                float centerx =  block_off_x[i] + cx; // centerx
                float centery =  block_off_y[j] + cy; // centery

                if(score >= para.thres)
                {                    
                    if(centerx>=para.d_m && centerx<nx-para.d_m &&centery>=para.d_m &&centery<ny-para.d_m)
                    {
                        score_info[0].push_back(score);
                        score_info[1].push_back(centerx);
                        score_info[2].push_back(centery);
                        score_info[3].push_back(J);
                    }
                }
            }
           
        }
    }
}


void writeScoreToDisk(vector<vector<float> > &score_info,Parameters para,EulerData euler,FILE *fp, int nn, char *t, float euler3)
{
    for(int i=0;i<score_info[0].size();i++)
    {
        float score = score_info[0][i];
        float centerx = score_info[1][i];
        float centery = score_info[2][i];
        int J = (int)score_info[3][i];
  //      printf( "%d\t%s\tdefocus=%f\tdfdiff=%f\tdfang=%f\teuler=%f,%f,%f\tcenter=%f,%f\tscore=%f\n",
  //      nn,t,(-1)*para.defocus,para.dfdiff,para.dfang,euler.euler1[J],euler.euler2[J],euler3,centerx,centery,score
//);
        fprintf(fp, "%d\t%s\tdefocus=%f\tdfdiff=%f\tdfang=%f\teuler=%f,%f,%f\tcenter=%f,%f\tscore=%f\n",
                nn,t,(-1)*para.defocus,para.dfdiff,para.dfang,euler.euler1[J],euler.euler2[J],euler3,centerx,centery,score
        );
    }
}



/**
 * @brief pathConvert_Double2Single  "\\" -> "/"
 * @param s
 */
 inline void pathConvert_Double2Single ( std::string& s )
 {
     std::string::size_type pos = 0;
     while ( ( pos = s.find('\\', pos) ) !=  std::string::npos )
     {
         s.replace( pos, 1,"/" );
         pos = pos + 2;
     }
 }
 
inline bool create_outpath(char * path)
{
    int ret=0;
    std::string tempStr(path);
    std::string forler;
 
    pathConvert_Double2Single ( tempStr );
    for(int pos=0;pos<tempStr.size();pos++)
    {
        if(tempStr[pos] != '/') continue;
        forler =  tempStr.substr( 0, pos);
        ret = mkdir( forler.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if(!ret)
            return false;
    }
    return true;
}

int main(int argc, char *argv[])
{
    //Timer
    time_t first, second;  

    //euler.length, number of template
    int N_tmp;

    //Print Help Message if no para input
    if(argc==1 || (argv[1][0] == '-' && argv[1][1] == 'h')){
		printHelpMsg(); 
		return 0;
    }

    //Store all parameters(some from input, some from computation)
    Parameters para;
    //Stroe Euler data
    EulerData euler;
    
    //readParameters(argc,argv,&para);
    readConfig(argv[1],&para);
    para.printAllPara();

    readEulerData(para.eulerf,&euler);
    N_tmp = euler.length; 

    //IMG id
    int nn=0;
    //IMG filename
    char t[MAXPATHLEN+1];
    //create write path
    int rtn = create_outpath(para.outlst);
    if(rtn < 0) 
    {
        printf("Out file create failed! %s\n",para.outlst);
        return 1;
    }
    //Used to write res
    FILE *fp=fopen(para.outlst,"wb");
    if(fp==NULL)
    {
        printf("Out file create failed! %s\n",para.outlst);
        printf("Please confirm the path exists!\n");
        return 1;
    }
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
    //Store sum of CCG
    cufftComplex *CCG_sum, *CCG_buf;
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
    //Store the offset of each block
    vector<int > block_off_x;
    vector<int > block_off_y;

    // READ FIRST IMG TO GET nx & ny
    //Image to be searched 
    emdata *first_image = new emdata();
    //Read first Image
    readRawImage(first_image,&para,para.first,&nn,t,1);
    //size of first IMG
    int nx = first_image->header.nx;
    int ny = first_image->header.ny;
#ifdef DEBUG
        printf("IMG size:          %d %d\n",nx,ny);
        printf("Number of template:%d\n",N_tmp);
#endif
    // check size of template, if nx/ny < padding_size , adjust padding size
    // granatee 1.padding_size > tmp.nx  2.padding_size > tmp.ny  3.padding_size is 32*K 
    adjustPaddingSize(&para);

    cudaAllocTemplateMem(N_tmp,nx,ny,block_off_x,block_off_y,&para,&h_reduction_buf,&d_reduction_buf,&d_means,&sigmas,&d_sigmas,
        &h_templates,&d_templates,&CCG,&CCG_sum,&CCG_buf,&stream,&ra,&rb,&plan_for_temp);
    
    //Reading Template
    readAndPaddingTemplate(&para,h_templates,N_tmp,sigmas);

    //Process edge normlization
    //To avoid alloc memory,use CCG as tempary buf
    edgeNormalize(N_tmp,sigmas,d_sigmas,CCG,h_templates,d_templates,h_reduction_buf,d_reduction_buf,d_means,para,&stream);

    //whiten, apply mask, appy weighting ..
    handleTemplate(N_tmp,ra,rb,h_reduction_buf,d_reduction_buf,d_means,
        d_templates,&para,&stream,&plan_for_temp);

    //free heap memory
    if(first_image!=NULL) delete first_image;
    
    //pointer to store image data
    float *d_image;
    cufftComplex *d_rotated_image;
    cufftComplex *rotated_splitted_image;
    //cufft handler for sub-IMGs and whole-IMG FFT
    cufftHandle plan_for_image,plan_for_whole_IMG;
    //pointer of image to be searched 
    emdata *image = new emdata();;
    
    //timer
    first=time(NULL);
    //flag for the first loop
    bool firstOrNot = true;
    //Loop for all Images. (Last - First) normally is 1;
    for(int n=para.first;n<para.last;n++)
    {
//*************************************************
// Process Image to be searched 
// 0. Read
// 1. Rotate
// 2. Split
// 3. doFFT
//*************************************************   
        //Read Image
        readRawImage(image,&para,n,&nn,t);
        if(firstOrNot) cudaAllocImageMem(&d_image,&d_rotated_image,&rotated_splitted_image,&stream,&plan_for_image,&plan_for_whole_IMG,nx,ny,N_tmp,&para);
        firstOrNot = false;
        // 1.Put Image on GPU 2.phaseflip 3.Replace Hot pixels With mean
        init_d_image(para,rotated_splitted_image,d_image,ra,rb,image,nx,ny,&stream,&plan_for_whole_IMG); 
        // split Image into blocks with overlap
        split_normalize_image(d_image,d_rotated_image,h_reduction_buf,d_reduction_buf,d_means,para,&stream,block_off_x,block_off_y,nx,ny,&plan_for_image);

        // fft on image
        fft_on_img(d_rotated_image,&plan_for_image,para,&stream);
        // normalize at all angles
        if( para.norm_type == 1 )
        {
            // h_templates used as a buffer of CCG result
            clear_ccg_sum(CCG_sum,para,&stream);
            
            int n_euler = 0;  
            for(float euler3=0.0;euler3<360.0;euler3+=para.phi_step) n_euler ++ ;

            printf("\nCompute the average and variance of all CCGs\n");
            // compute CCG sum
            for(float euler3=0.0;euler3<360.0;euler3+=para.phi_step)
            {   
                printf("CCG_sum of CC&C^2 in euler3:  %f / 360.0\n",euler3);
                rotateTemplate(&plan_for_temp,d_templates,CCG_buf,h_reduction_buf,d_reduction_buf,d_means,para,N_tmp,euler3,&stream);
               // rotateImage(d_rotated_image,rotated_splitted_image,para,euler3,&stream);
               // fft_on_img(rotated_splitted_image,&plan_for_image,para,&stream);
                //compute CCG sum from IMG with template
                //1. Calculate ccg
                //2. Calculate avg and variance
                //***************************************************

               // compute_CCG_sum(CCG,CCG_sum,CCG_buf,d_rotated_image,para,&plan_for_temp,&stream,N_tmp,euler3,n_euler);
               compute_CCG_sum(CCG,CCG_sum,CCG_buf,d_rotated_image,para,&plan_for_temp,&stream,N_tmp,euler3,n_euler);
            }
            
            compute_CCG_mean(CCG_sum,para, N_tmp,n_euler, &stream);
            // debug
            // out_mean_var_matrix(para,CCG_sum,&stream);break;

            printf("\nUpdate CCGs and compute scores\n");
            // compute CCG & (CCG-avg)/sqrt(s2)
            for(float euler3=0.0;euler3<360.0;euler3+=para.phi_step)
            {
                //Score_info: [0]->score  [1]->cx  [2]->cy [3]->id of tmp 
                vector<vector<float> > score_info(4,vector<float>());	

                printf("Now euler3:  %f / 360.0\n",euler3);
                // rotate tmp at CCG buf
                rotateTemplate(&plan_for_temp,d_templates,CCG_buf,h_reduction_buf,d_reduction_buf,d_means,para,N_tmp,euler3,&stream);
                //rotateImage(d_rotated_image,rotated_splitted_image,para,euler3,&stream);
                //fft_on_img(rotated_splitted_image,&plan_for_image,para,&stream);
                //***************************************************
                //Pick particles from IMG with template
                //1. Calculate ccg
                //2. (CCG-avg)/sqrt(s2)
                //3. Output score
                //***************************************************

                //h_templates used as a buffer of CCG result
               // pickPartcles_CCG_NORM(CCG,CCG_sum,CCG_buf,d_rotated_image,h_reduction_buf,d_reduction_buf,
               //     para,&plan_for_temp,&stream,N_tmp,score_info,block_off_x,block_off_y,nx,ny,euler3);
                pickPartcles_CCG_NORM(CCG,CCG_sum,CCG_buf,d_rotated_image,h_reduction_buf,d_reduction_buf,
                    para,&plan_for_temp,&stream,N_tmp,score_info,block_off_x,block_off_y,nx,ny,euler3);
                cudaStreamSynchronize(stream);

                float cur_e =  360.0f-euler3;
                if(cur_e >= 360) cur_e -= 360;
                writeScoreToDisk(score_info,para,euler,fp,nn,t,cur_e);
            }

        }else if(para.norm_type == 0)
        {
            for(float euler3=0.0;euler3<360.0;euler3+=para.phi_step)
            {	 
                //Score_info: [0]->score  [1]->cx  [2]->cy [3]->id of tmp 
                vector<vector<float> > score_info(4,vector<float>());	

                printf("Now euler3:  %f / 360.0\n",euler3);
                // rotate Image  
                //rotateTemplate(&plan_for_temp,d_templates,CCG_buf,para,N_tmp,euler3,&stream);
                rotateTemplate(&plan_for_temp,d_templates,CCG_buf,h_reduction_buf,d_reduction_buf,d_means,para,N_tmp,euler3,&stream);
    
                //***************************************************
                //Pick particles from IMG with template
                //1. Calculate ccg
                //2. Find peak 
                //3. Calculate variance
                //4. Output score
                //***************************************************
    
                pickPartcles_IMG_NORM(CCG,CCG_buf,d_rotated_image,h_reduction_buf,d_reduction_buf,para,&plan_for_temp,&plan_for_image,
                    &stream,N_tmp,score_info,block_off_x,block_off_y,nx,ny,euler3);
                cudaStreamSynchronize(stream);
                
                float cur_e =  360.0f-euler3;
                if(cur_e >= 360) cur_e -= 360;
                writeScoreToDisk(score_info,para,euler,fp,nn,t,cur_e);
            }

        }
        else{
            printf("Unknown norm type!\n");
            exit(-1);
        }
        
  
        //check avail mem, for debugging
        //checkGPUMem();
        printf("------------------------Image %d finished--------------------------------\n",n);
    }
    //Timer
    second=time(NULL);  
    printf("Total consumed time is: %f seconds\n",difftime(second,first)); 
    
    cufftDestroy(plan_for_whole_IMG);
    cufftDestroy(plan_for_image);
    cudaFree(rotated_splitted_image);
    cudaFree(d_rotated_image);
    cudaFree(d_image);
    //free heap memory
    if(image!=NULL) delete image;

    cufftDestroy(plan_for_temp);
    cudaStreamDestroy(stream);
    cudaFree(d_reduction_buf);
    cudaFree(d_templates);
    cudaFree(d_sigmas);
    cudaFree(d_means);
    cudaFree(CCG);
    cudaFree(CCG_sum);
    cudaFree(ra);
    cudaFree(rb);
    free(h_reduction_buf);
    free(h_templates);
    free(sigmas);
    fclose(fp);


    return 0;
}
