#include "emdata_.h"
#include "util_func.h"
#include "DataReader.h"
#include <cuda_runtime.h>
#include <cufft.h>
#include "GPU_func.cuh"
#include <malloc.h>

void readEMData(Parameters *para, EulerData *euler)
{
    readEulerData(para->eulerf,euler);
	readSNRWeight(para->snr,&para->a,&para->b,&para->b2,&para->bfactor,&para->bfactor2,&para->bfactor3);
}

void readRawImage(emdata *i2d_emdata,Parameters *para,int n, int *nn, char* t)
{
    vector<string> pairs;
    int return_code = readInLst_and_consturctPairs(para->inlst,t,&pairs,nn,n);
    if( return_code < 0) printf("Error %d occured in readRawImage!\n",return_code);
    
    i2d_emdata->readImage(t,*nn);
    parsePairs(t,*nn,pairs, &para->defocus, &para->dfdiff, &para->dfang);
    
    para->dfu=para->defocus+para->dfdiff; //defocus is minus, so abs(dfu) < abs(dfv)
    para->dfv=para->defocus-para->dfdiff;
    para->lambda=12.2639/sqrt(para->energy*1000.0+0.97845*para->energy*para->energy);
    para->ds=1/(para->apix * i2d_emdata->header.ny);
}

//read template(float *)
//covert template from float* to cufftcomplex *
void readAndPaddingTemplate(Parameters *para,cufftComplex *h_templates,int N)
{
    int padded_template_size = para->padding_size * para->padding_size;
    emdata *tp = new emdata();
    for(int J=0;J<N;J++)
    {
        tp->readImage(para->temp2d,J);
        float *data = tp->getData();
        if(para->padding_size < tp->header.nx)
        {
            printf("Padded size is smaller than template\n");
            return ;
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
    para->overlap = para->template_x*0.13+1;
    free(tp);
}

void cudaAllocTemplateMem(int N, Parameters para,float **h_reduction_buf,float **d_reduction_buf,float **d_means,cufftComplex **h_templates,
    cufftComplex **d_templates,cufftComplex **CCG,float **ra, float **rb, cufftHandle *plan_for_temp)
{
    //Number of Pixels for each padded template
    long long padded_template_size = para.padding_size*para.padding_size;

    //All padded templates (complex) At CPU,GPU
    *h_templates = (cufftComplex *)malloc(sizeof(cufftComplex)*padded_template_size*N);
    memset(*h_templates,0,sizeof(cufftComplex)*padded_template_size*N);
    CUDA_CALL(  cudaMalloc(d_templates,sizeof(cufftComplex)*padded_template_size*N)  );
    
    //Memory alloc for CCG
    CUDA_CALL(  cudaMalloc(CCG,sizeof(cufftComplex)*padded_template_size*N)  );

    //Temp buffer for whiten
    CUDA_CALL(  cudaMalloc(ra,N*(para.padding_size/2)*sizeof(float))  );
    CUDA_CALL(  cudaMalloc(rb,N*(para.padding_size/2)*sizeof(float))  );

    //Buffer for reduction
    *h_reduction_buf = (float *)malloc(4*sizeof(float)*padded_template_size*N/1024);
    CUDA_CALL(  cudaMalloc(d_reduction_buf,4*sizeof(float)*padded_template_size*N/1024)  );

    //Store mean for every template
    CUDA_CALL(  cudaMalloc(d_means,sizeof(float)*N)  );

    /*
	cufftMakePlanMany(cufftHandle plan, int rank, int *n, int *inembed,
	int istride, int idist, int *onembed, int ostride,
	int odist, cufftType type, int batch, size_t *workSize);
	 */
	const int rank = 2;//维数
    int n[rank] = { para.padding_size, para.padding_size };//n*m
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

}

void processTemplate(int N, float *ra, float *rb,float *h_inner_sum,float *d_inner_sum,float *d_means,
    cufftComplex *h_templates,cufftComplex *d_templates,Parameters *para,cufftHandle *plan_for_temp)
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

    CUDA_CALL(  cudaMemcpy(d_templates, h_templates, sizeof(cufftComplex)*padded_template_size*N, cudaMemcpyHostToDevice)  );
    // Inplace FFT
    CUFFT_CALL(  cufftExecC2C(*plan_for_temp, d_templates, d_templates, CUFFT_FORWARD)  );
    // CUFFT will enlarge VALUE to N times. Restore it
    scale<<<block_num,BLOCK_SIZE>>>(d_templates,padded_template_size);
    CUDA_CHECK();
    CUDA_CALL(  cudaMemcpy(h_templates, d_templates, sizeof(cufftComplex)*padded_template_size*N, cudaMemcpyDeviceToHost)  );

    //Whiten at fourier space
    SQRSum_by_circle<<<block_num,BLOCK_SIZE>>>(d_templates,ra,rb,para->padding_size);
    CUDA_CHECK();
    whiten_Img<<<block_num,BLOCK_SIZE>>>(d_templates,ra,rb,para->padding_size);
    CUDA_CHECK();

// **************************************************************
// apply mask 
// input: whiten_IMAGE (Fourier SPACE in RI)
// output: masked_whiten_IMAGE (Fourier SPACE in RI)
// **************************************************************
    CUFFT_CALL(  cufftExecC2C(*plan_for_temp, d_templates, d_templates, CUFFT_INVERSE)  );
    apply_mask<<<block_num,BLOCK_SIZE>>>(d_templates,para->d_m,para->edge_half_width,para->padding_size);
    CUDA_CHECK();
    CUFFT_CALL(  cufftExecC2C(*plan_for_temp, d_templates, d_templates, CUFFT_FORWARD)  );
    // CUFFT will enlarge VALUE to N times. Restore it
    scale<<<block_num,BLOCK_SIZE>>>(d_templates,padded_template_size);
    CUDA_CHECK();

// **************************************************************
// 1. lowpass
// 2. apply weighting function
// 3. normlize
// input: masked_whiten_IMAGE (Fourier SPACE in RI) 
// output: PROCESSED_IMAGE (Fourier SPACE in AP)
// **************************************************************

    apply_weighting_function<<<block_num,BLOCK_SIZE>>>(d_templates,*para);
    CUDA_CHECK();
    compute_area_sum_ofSQR<<<block_num,BLOCK_SIZE,2*BLOCK_SIZE*sizeof(float)>>>(d_templates,d_inner_sum,para->padding_size);
    CUDA_CHECK();
    CUDA_CALL(cudaMemcpy(h_inner_sum, d_inner_sum,2*sizeof(float)*padded_template_size*N/1024, cudaMemcpyDeviceToHost));
    //After Reduction -> compute mean for each image
    float infile_mean[N],counts[N];
    memset(infile_mean,0,N*sizeof(float));
    memset(counts,0,N*sizeof(float));
    for(int k=0;k<padded_template_size*N/1024;k++)
    {
        int id = k/(padded_template_size/1024);
        infile_mean[id] += h_inner_sum[2*k];
        counts[id] += h_inner_sum[2*k+1];
    }
    for(int k=0;k<N;k++) infile_mean[k] /= counts[k];
    //Do Normalization with computed infile_mean[]
    CUDA_CALL(  cudaMemcpy(d_means, infile_mean, sizeof(float)*N, cudaMemcpyHostToDevice)  );
    normalize<<<block_num,BLOCK_SIZE>>>(d_templates,N,d_means);
    CUDA_CHECK();

}

void cudaAllocImageMem(float **d_image,cufftComplex **d_rotated_image,
    cufftHandle *plan_for_image,int nx,int ny,Parameters *para)
{
    int tmp = (para->padding_size - para->overlap);
    //num of blocks in x,y axis
    int block_x = (nx-para->overlap) / tmp;
    if((nx-para->overlap) % tmp > 0 ) block_x ++; 
    int block_y = (ny-para->overlap) / tmp;
    if((nx-para->overlap) % tmp > 0 ) block_y ++; 
    //X Y Size for padded IMG
    int ix = block_x*para->padding_size;
    int iy = block_y*para->padding_size;
    para->block_x = block_x;
    para->block_y = block_y;

    CUDA_CALL(  cudaMalloc(d_image,nx*ny*sizeof(float))  );
    CUDA_CALL(  cudaMalloc(d_rotated_image,ix*iy*sizeof(cufftComplex))  );

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
    int batch = block_x*block_y;//批量处理的批数
    //采用cufftPlanMany方法
    
    //FFT handler for all templates
    CUFFT_CALL(cufftPlanMany(plan_for_image, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, batch));//针对多信号同时进行FFT

}  

void init_d_image(float *d_image, emdata *image, int nx, int ny)
{
    //Translate origin of Image to (0,0)
    image->rotate(0);
    //Put Image on GPU
    CUDA_CALL(  cudaMemcpy(d_image, image->getData(), sizeof(float)*nx*ny, cudaMemcpyHostToDevice)  );
}

void processImage(float *d_image,cufftComplex *d_rotated_image,Parameters para, float e, int nx, int ny)
{
    int block_num = nx*ny/BLOCK_SIZE;
    //Init d_rotated_imge to all {0}
    int ix = para.block_x;
    int iy = para.block_y;
    int blockIMG_num = ix*iy*para.padding_size*para.padding_size / BLOCK_SIZE;
    //printf("%d %d %lld %lld\n",ix,iy,block_num,blockIMG_num);
    clear_image<<<blockIMG_num,BLOCK_SIZE>>>(d_rotated_image);
    CUDA_CHECK();
    // 1.rotate Image  2.split Image into blocks with overlap
    rotate_and_split<<<block_num,BLOCK_SIZE>>>(d_image,d_rotated_image,e,nx,ny,para.padding_size,para.block_x,para.overlap);
    CUDA_CHECK();	
}

void pickPartcles(cufftComplex *CCG,cufftComplex *d_templates,cufftComplex *d_rotated_image,Parameters para, 
    cufftHandle *template_plan, cufftHandle *image_plan,int N,float *h_buf,float *d_buf,float *scores)
{           

    long long padded_template_size = para.padding_size*para.padding_size;
    int blockGPU_num = padded_template_size*N/BLOCK_SIZE;
    int blockIMG_num = padded_template_size*para.block_x*para.block_y/BLOCK_SIZE;

    //peak,sum of data[i],sum of data[i]^2
    float peaks[N],pos[N],sums[N],sum2s[N];

    //find MAX score need initialize
    memset(scores,0,sizeof(float)*2*N);

    //Inplace FFT
    CUFFT_CALL(  cufftExecC2C(*image_plan, d_rotated_image, d_rotated_image, CUFFT_FORWARD)  );

    //Scale IMG to normel size
    scale<<<blockIMG_num,BLOCK_SIZE>>>(d_rotated_image,para.padding_size*para.padding_size);
    CUDA_CHECK();

    //compute score for each block
    for(int j=0;j<para.block_y;j++)
        for(int i=0;i<para.block_x;i++)
        {
            //find peak need initialize
            memset(peaks,0,sizeof(float)*N);
            memset(pos,0,sizeof(float)*N);
            memset(sums,0,sizeof(float)*N);
            memset(sum2s,0,sizeof(float)*N);
            //printf("%d %d \n",i,j);
            //compute CCG
            //printf("%d %d %d %d\n",blockGPU_num,BLOCK_SIZE,i+j*para.block_x,para.block_x);
	        compute_CCG<<<blockGPU_num,BLOCK_SIZE>>>(CCG,d_templates,d_rotated_image,para.padding_size,i+j*para.block_x);
            CUDA_CHECK();
            cudaDeviceSynchronize();
            
            //Inplace IFT
            CUFFT_CALL(  cufftExecC2C(*template_plan, CCG, CCG, CUFFT_INVERSE)  );
            //find peak(position) and get sum of data,data^2
            get_peak_and_SUM<<<blockGPU_num,BLOCK_SIZE,4*BLOCK_SIZE*sizeof(float)>>>(CCG,d_buf,para.padding_size,para.d_m);
            CUDA_CHECK();
            CUDA_CALL(cudaMemcpy(h_buf, d_buf, 4*sizeof(float)*padded_template_size*N/1024, cudaMemcpyDeviceToHost));

            //After Reduction -> compute mean for each image
            for(int k=0;k<padded_template_size*N/1024;k++)
            {
                int id = k/(padded_template_size/1024);
                if(peaks[id] < h_buf[4*k])
                {
                    peaks[id] = h_buf[4*k];
                    pos[id] = h_buf[4*k+1];
                }
                sums[id] += h_buf[4*k+2];
                sum2s[id] += h_buf[4*k+3];
            }
            
            //Update global score with local-block score for each template
            for(int k=0;k<N;k++) 
            {
                float ra = sums[k];
                float rb = sum2s[k];
                float rc = padded_template_size - 1;
                float score = peaks[k]/sqrt(rb/rc - (ra/rc)*(ra/rc));
                if(scores[2*k] < score)
                {
                    scores[2*k] = score;
                    scores[2*k+1] = pos[k];
                }
            }
        }

}

void writeScoreToDisk(float *scores,Parameters para,EulerData euler,FILE *fp, int nn, char *t, int nx, int ny, float euler3)
{
    for(int J=0;J<euler.length;J++)
    {
        int cx = (int)scores[2*J+1] % para.padding_size;
        int cy = (int)scores[2*J+1] / para.padding_size;
        //Rotate (cx,cy) to its soriginal angle
        float centerx = (cx-nx/2)*cos(euler3*PI/180)+(cy-ny/2)*sin(euler3*PI/180)+nx/2; // centerx
		float centery = (cy-ny/2)*cos(euler3*PI/180)-(cx-nx/2)*sin(euler3*PI/180)+ny/2; // centery
        float score = scores[2*J];

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
    readEMData(&para,&euler);

#ifdef DEBUG
    para.printAllPara();
#endif

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
    //cufft handler for template FFT
    cufftHandle plan_for_temp;
    cudaAllocTemplateMem(euler.length,para,&h_reduction_buf,&d_reduction_buf,&d_means,
        &h_templates,&d_templates,&CCG,&ra,&rb,&plan_for_temp);

    //Loop for all Images. (Last - First) normally is 1;
    for(int n=para.first;n<para.last;n++)
	{
        //Image to be searched 
        emdata *image = new emdata();
        //Read Image
        readRawImage(image,&para,n,&nn,t);
        //Reading Template
        readAndPaddingTemplate(&para,h_templates,euler.length);

        processTemplate(euler.length,ra,rb,h_reduction_buf,d_reduction_buf,d_means,
            h_templates,d_templates,&para,&plan_for_temp);

//*************************************************
// Process Image to be searched 
// 0. Read
// 1. Rotate
// 2. Split
// 3. doFFT
//*************************************************   
        float *d_image;
        cufftComplex *d_rotated_image;
        //cufft handler for IMG FFT
        cufftHandle plan_for_image;
        //size of IMG
        int nx = image->header.nx;
        int ny = image->header.ny;

#ifdef DEBUG
        printf("IMG size: %d %d\n",nx,ny);
#endif
        
        cudaAllocImageMem(&d_image,&d_rotated_image,&plan_for_image,nx,ny,&para);
        //Put Image on GPU
        init_d_image(d_image,image,nx,ny);

        //Scores : [1-N]->socre  [N+1-2N]->cx+cy*padding_size
        float *scores = new float[euler.length*2];    
        for(float euler3=0.0;euler3<360.0;euler3+=para.phi_step)
        {	 
#ifdef DEBUG
            printf("Now euler3 => %f / 360.0\n",euler3);
#endif
            processImage(d_image,d_rotated_image,para,euler3,nx,ny);
            //***************************************************
            //Pick particles from IMG with template
            //1. Calculate ccg
            //2. Find peak 
            //3. Calculate variance
            //4. Output score
            //***************************************************

            pickPartcles(CCG,d_templates,d_rotated_image,para,&plan_for_temp,&plan_for_image,euler.length,h_reduction_buf,d_reduction_buf,scores);
            //writeScoreToDisk(scores,para,euler,fp,nn,t,nx,ny,euler3);
        }
        cufftDestroy(plan_for_image);
        cudaFree(d_rotated_image);
        cudaFree(d_image);
        delete scores;
        
    }

    cufftDestroy(plan_for_temp);
    cudaFree(d_reduction_buf);
    cudaFree(d_templates);
    cudaFree(d_means);
    cudaFree(CCG);
    cudaFree(ra);
    cudaFree(rb);
    free(h_reduction_buf);
    free(h_templates);
    fclose(fp);

    //Timer
    second=time(NULL);  
    printf("Consumed time is: %f seconds\n",difftime(second,first)); 

    return 0;
}
