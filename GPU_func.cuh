#ifndef GPU_FUNC_H
#define GPU_FUNC_H
#include "cuda_runtime.h"
#include "cufft.h"
#include <device_launch_parameters.h>
#include "EMReader/DataReader.h"

/*
cufftResult_t :
    CUFFT_SUCCESS        = 0,  //  The cuFFT operation was successful
    CUFFT_INVALID_PLAN   = 1,  //  cuFFT was passed an invalid plan handle
    CUFFT_ALLOC_FAILED   = 2,  //  cuFFT failed to allocate GPU or CPU memory
    CUFFT_INVALID_TYPE   = 3,  //  No longer used
    CUFFT_INVALID_VALUE  = 4,  //  User specified an invalid pointer or parameter
    CUFFT_INTERNAL_ERROR = 5,  //  Driver or internal cuFFT library error
    CUFFT_EXEC_FAILED    = 6,  //  Failed to execute an FFT on the GPU
    CUFFT_SETUP_FAILED   = 7,  //  The cuFFT library failed to initialize
    CUFFT_INVALID_SIZE   = 8,  //  User specified an invalid transform size
    CUFFT_UNALIGNED_DATA = 9,  //  No longer used
    CUFFT_INCOMPLETE_PARAMETER_LIST = 10, //  Missing parameters in call
    CUFFT_INVALID_DEVICE = 11, //  Execution of a plan was on different GPU than plan creation
    CUFFT_PARSE_ERROR    = 12, //  Internal plan database error 
    CUFFT_NO_WORKSPACE   = 13  //  No workspace has been provided prior to plan execution
    CUFFT_NOT_IMPLEMENTED = 14, // Function does not implement functionality for parameters given.
    CUFFT_LICENSE_ERROR  = 15, // Used in previous versions.
    CUFFT_NOT_SUPPORTED  = 16  // Operation is not supported for parameters given.
*/
#define EPS 1e-8
#define PI (3.141592654)
#define BLOCK_SIZE 1024
#define RA_SIZE 5000
#ifdef DEBUG
#define CUDA_CALL(F)  if( (F) != cudaSuccess ) \
  {printf("Error %s at %s:%d\n", cudaGetErrorString(cudaGetLastError()), \
   __FILE__,__LINE__); exit(-1);}
#define CUFFT_CALL(F)  if( (F) != CUFFT_SUCCESS ) \
  {printf("Error id:%d at %s:%d\n", F, \
   __FILE__,__LINE__); exit(-1);}
#define CUDA_CHECK()  if( (cudaPeekAtLastError()) != cudaSuccess ) \
  {printf("Error %s at %s:%d\n", cudaGetErrorString(cudaGetLastError()), \
   __FILE__,__LINE__-1); exit(-1);}
#else
#define CUDA_CALL(F) (F)
#define CUDA_CHECK()
#endif
__global__ void UpdateSigma(cufftComplex *d_templates,float *d_buf);
__global__ void generate_mask(int l,cufftComplex *mask,float r,float *d_buf,float up,float low);
__global__ void multiCount_dot(int l,cufftComplex *mask,cufftComplex *d_templates,float *constants,float *res);
__global__ void scale_each(int l,cufftComplex *d_templates,float *ems,double *d_sigmas);
__global__ void SQRSum_by_circle(cufftComplex *data, float *ra, float *rb, int nx, int ny, int mode=0);
__global__ void whiten_Tmp(cufftComplex *data, float *ra, float *rb, int l);
__global__ void whiten_filetr_weight_Img(cufftComplex *data, float *ra, float *rb, int nx, int ny, Parameters para);
__global__ void apply_mask(cufftComplex *data,float d_m,float edge_half_width,int l);
__global__ void apply_weighting_function(cufftComplex *data,Parameters para);
__global__ void compute_area_sum_ofSQR(cufftComplex *data,float *res,int nx, int ny,int mode = 0);
__global__ void normalize(cufftComplex *d_templates,int nx, int ny,float *means);
__global__ void normalize_Img(cufftComplex *data,int nx, int ny,float mean);
__global__ void rotate_IMG(float *d_image,float *d_rotated_image,float e,int nx,int ny);
__global__ void split_IMG(float *Ori,cufftComplex *IMG, int nx,int ny,int l,int bx,int overlap);
__global__ void compute_corner_CCG(cufftComplex *CCG, cufftComplex *Tl, cufftComplex *IMG, int l, int block_id);
__global__ void get_peak_and_SUM(cufftComplex *odata,float *res,int l,float d_m,int x_bound, int y_bound);
__global__ void scale(cufftComplex *data,int l2);
__global__ void clear_image(cufftComplex *data);
__global__ void Complex2float(float *f, cufftComplex *c, int nx, int ny);
__global__ void float2Complex(cufftComplex *c, float *f, int nx, int ny);
__global__ void do_phase_flip(cufftComplex *filter, Parameters para, int nx, int ny);

__device__ float CTF_AST(int x1, int y1, int nx, int ny, float apix, float ds, float dfu, float dfv, float dfdiff, float dfang ,float lambda, float cs, float ampconst, int mode);

void cudaMemoryTest();
#endif
