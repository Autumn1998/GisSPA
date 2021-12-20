#include "GPU_func.cuh"

__global__ void UpdateSigma(cufftComplex *d_templates,float *d_buf)
{
	extern __shared__ float sdata[];
    // each thread loads one element from global to shared mem
    int tid = threadIdx.x;
    long long  i = blockIdx.x*blockDim.x + threadIdx.x;

	sdata[tid] = d_templates[i].x;
	sdata[tid+blockDim.x] = d_templates[i].x*d_templates[i].x;
	__syncthreads();

	for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
		if(tid<s)
		{
			//sum of data[i] & data[i]^2
			sdata[tid] += sdata[tid + s]; 
			sdata[tid+blockDim.x] += sdata[tid + blockDim.x + s];
		}
		__syncthreads();
	}
	if(tid==0){
		d_buf[blockIdx.x*2] = sdata[0];
		d_buf[blockIdx.x*2+1] = sdata[blockDim.x];
	}
}

__global__  void  generate_mask(int l,cufftComplex *mask,float r,float *res,float up,float low)
{
	extern __shared__ float sdata[];
    // each thread loads one element from global to shared mem
    int tid = threadIdx.x;
    long long  i = blockIdx.x*blockDim.x + threadIdx.x;
	int image_size = l*l;
	int local_id = i % image_size;
    int x = local_id % l;
    int y = local_id / l;

	//Dis^2 between (x,y) and center (l/2,l/2)
	float rr = (x-l/2)*(x-l/2) + (y-l/2)*(y-l/2);
	if(rr>=low && rr<=up)
	{
		mask[i].x = 1;
		mask[i].y = 0;
	}

	//reduction for number of non-zero digits
	sdata[tid] = mask[i].x;
	__syncthreads();

	for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
		if(tid<s)	sdata[tid] += sdata[tid+s];
		__syncthreads();
	}
	if(tid==0) res[blockIdx.x] = sdata[0];

}

__global__  void  multiCount_dot(int l,cufftComplex *mask,cufftComplex *d_templates,float *constants,float *res)
{
	extern __shared__ float sdata[];
    // each thread loads one element from global to shared mem
    int tid = threadIdx.x;
    long long  i = blockIdx.x*blockDim.x + threadIdx.x;
	int image_size = l*l;
	int image_id = i/image_size;

	//Multi constant 1/non-zeros
	if(constants[image_id] != 0) mask[i].x *= 1.0/constants[image_id];

	//reduction for dot
	sdata[tid] = mask[i].x * d_templates[i].x;
	__syncthreads();

	for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
		if(tid<s)	sdata[tid] += sdata[tid+s];
		__syncthreads();
	}
	//use res to store dot
	if(tid==0) res[blockIdx.x] = sdata[0];

}

__global__  void  scale_each(int l,cufftComplex *d_templates,float *ems,double *d_sigmas)
{
    long long  i = blockIdx.x*blockDim.x + threadIdx.x;
	int image_size = l*l;
	int image_id = i/image_size;

	if(d_sigmas[image_id]-0 < EPS && d_sigmas[image_id]-0 >-EPS ) return;
	d_templates[i].x = (d_templates[i].x - ems[image_id])/d_sigmas[image_id];
}

//mode = 0 (default)  for template
//mode = 1 for raw image
__global__  void  SQRSum_by_circle(cufftComplex *data, float *ra, float *rb, int nx, int ny, int mode)
{
    // i <==> global ID
    long long  i = blockIdx.x*blockDim.x + threadIdx.x;
    int image_size = nx*ny;
	if(mode == 1 && i>=image_size) return;
    int local_id = i % image_size;
    int x = local_id % nx;
    int y = local_id / nx;

    float tmp;
	// ri2ap
	tmp=hypotf(data[i].x, data[i].y);
	if (data[i].x==0 && data[i].y==0) 
        data[i].y=0;
	else data[i].y=atan2(data[i].y,data[i].x);
	data[i].x=tmp;

	if(x>nx/2) return;

    //calculate the number of point with fixed distance ('r') from center 
	int r = floor( hypotf(min(y,ny-y) ,min(x,nx-x)) + 0.5) - 1;

	if (r < max(nx,ny)/2 && r >= 0) {
		//Add offset
		r+= RA_SIZE*(i/image_size);
		atomicAdd(&ra[r],data[i].x*data[i].x);
		atomicAdd(&rb[r],1.0);
	}

}

__global__  void  whiten_Tmp(cufftComplex *data, float *ra, float *rb, int l)
{
    // i <==> global ID
    long long  i = blockIdx.x*blockDim.x + threadIdx.x;
    int image_size = l*l;  
    int local_id = i % image_size;
    int x = local_id % l;
    int y = local_id / l;
	int r = floor( hypotf(min(y,l-y) ,min(x,l-x)) + 0.5) - 1;

	if (r < l/2 && r >= 0) {
		//Add offset
		r+= RA_SIZE*(i/image_size);
		float fb_infile=ra[r]/rb[r];
		data[i].x=data[i].x/(float)sqrt(fb_infile);
	}

	//ap2ri
	float tmp=data[i].x*sinf(data[i].y);
	data[i].x=data[i].x*cosf(data[i].y);
	data[i].y=tmp;

}


__global__  void  whiten_filetr_weight_Img(cufftComplex *data, float *ra, float *rb, int nx, int ny, Parameters para)
{
    // i <==> global ID
    long long  i = blockIdx.x*blockDim.x + threadIdx.x;
    int image_size = nx*ny;  
	if(i>=image_size) return;
    int local_id = i % image_size;
    int x = local_id % nx;
    int y = local_id / nx;

	float rf = hypotf(min(y,ny-y) ,min(x,nx-x));
	int r = floor( rf + 0.5) - 1;
	float ss=rf*para.ds;
	int l = max(nx,ny);

	float v,signal,Ncurve;
    //apply weighting function
	if( r < l/2 && r >= 0){
		v=CTF_AST(x,(y+ny/2)%ny,nx,ny,para.apix,para.ds,para.dfu,para.dfv,para.dfdiff,para.dfang,para.lambda,para.cs,para.ampconst,2);
		signal=(exp(para.bfactor*ss*ss+para.bfactor2*ss+para.bfactor3));
		Ncurve=exp(para.a*ss*ss+para.b*ss+para.b2);

		data[i].x=data[i].x*sqrt((signal*v*v+Ncurve)/signal)/sqrt(ra[r]/rb[r]);
		if(r>(l*para.apix/6)) data[i].x=data[i].x*exp(-100*ss*ss);
	}

	// low pass
	if (r<l*para.apix/para.highres && r >= l*para.apix/para.lowres) {}
	else if(r>=l*para.apix/para.highres && r<l*para.apix/para.highres+8){
		data[i].x=data[i].x*(0.5*cosf(PI*(r-l*para.apix/para.highres)/(2*8))+0.5);
	}
	else if(r>=(l*para.apix/para.lowres-8) && r<l*para.apix/para.lowres && r>=0){
		data[i].x=data[i].x*(0.5*cosf(PI*(l*para.apix/para.lowres-r)/(2*8))+0.5);
	}
	else
		data[i].x=0;

    //apply weighting function
	if( r < l/2 && r >= 0){
		signal/=(para.kk+1);
		Ncurve/=signal;
		data[i].x=data[i].x*sqrt(1/(Ncurve+para.kk*v*v ));
	}
}

__global__ void normalize_Img(cufftComplex *data,int nx, int ny,float mean)
{
    // i <==> global ID
    long long  i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= nx*ny) return;

	if(mean!=0)	data[i].x=data[i].x/mean;

	//ap2ri
	float tmp=data[i].x*sinf(data[i].y);
	data[i].x=data[i].x*cosf(data[i].y);
	data[i].y=tmp;
}

__global__ void apply_mask(cufftComplex *data,float d_m,float edge_half_width,int l)
{	
    // i <==> global ID
    long long  i = blockIdx.x*blockDim.x + threadIdx.x;
    int image_size = l*l;
    int local_id = i % image_size;
    int x = local_id % l;
    int y = local_id / l;

    float r=hypotf(x-l/2,y-l/2);
	if( r > (d_m/2+2*edge_half_width)){
			data[i].x=0;
	}else if (r >= d_m/2){
			float d=0.5*cosf(PI*(r-d_m/2)/(2*edge_half_width))+0.5;
			data[i].x *=d;
	}

}

__global__ void apply_weighting_function(cufftComplex *data,Parameters para)
{
    int l = para.padding_size;
    // i <==> global ID
    long long  i = blockIdx.x*blockDim.x + threadIdx.x;
    int image_size = l*l;
    int local_id = i % image_size;
    int x = local_id % l;
    int y = local_id / l;

    float tmp;
	//ri2ap
	tmp=hypotf(data[i].x, data[i].y);
	if (data[i].x==0 && data[i].y==0) 
        data[i].y=0;
	else data[i].y=atan2(data[i].y,data[i].x);
	data[i].x=tmp;

    // low pass
	float r = hypotf(min(y,l-y) ,min(x,l-x) );
	int  r_round = floor(r + 0.5) - 1;
	if (r_round<l*para.apix/para.highres && r_round >= 0) {}
	else if(r_round>=l*para.apix/para.highres && r_round<l*para.apix/para.highres+8){
		data[i].x=data[i].x*(0.5*cosf(PI*(r_round-l*para.apix/para.highres)/(2*8))+0.5);
	}
	else if(r_round>=(l*para.apix/para.lowres-8) && r_round<l*para.apix/para.lowres && r_round>=0){
		data[i].x=data[i].x*(0.5*cosf(PI*(l*para.apix/para.lowres-r_round)/(2*8))+0.5);
	}
	else
		data[i].x=0;
	float ss=r*para.ds;

	float v,signal,Ncurve;
    //apply weighting function
	if( r_round < l/2 && r_round >= 0){
		v=CTF_AST(x,(y+l/2)%l,l,l,para.apix,para.ds,para.dfu,para.dfv,para.dfdiff,para.dfang,para.lambda,para.cs,para.ampconst,2);
		signal=(exp(para.bfactor*ss*ss+para.bfactor2*ss+para.bfactor3))/(para.kk+1);
		Ncurve=exp(para.a*ss*ss+para.b*ss+para.b2)/signal;
		//euler_w[x]=1.68118*ss;
		data[i].x=data[i].x*v*sqrt(1/(Ncurve+para.kk*v*v ));
	}
	
}

__device__ float CTF_AST (int x1, int y1,int nx, int ny, float apix, float ds, float dfu, float dfv, float dfdiff, float dfang ,float lambda, float cs, float ampconst, int mode){
	float v,ss,ag,gamma,df_ast;
	if(x1 > nx/2) 
	{
		x1 = nx-x1;
		y1 = ny-y1;
	}
	ss = ( x1*x1/(float)(nx*nx) + y1*y1/(float)(ny*ny) )/ (apix*apix);
	//ss = hypotf((float)x1,(float)y1-ny/2)*ds*hypotf((float)x1,(float)y1-ny/2)*ds;
	ag=atan2(float(y1-ny/2),float(x1));
	df_ast=0.5*(dfu+dfv+2*dfdiff*cosf(2*(dfang*PI/180-ag)));
	gamma=-2*PI*(cs*2.5e6*lambda*lambda*lambda*ss*ss+df_ast*5000.0*lambda*ss);
	if (mode==0){
		v=(sqrtf(1.0-ampconst*ampconst)*sinf(gamma)+ampconst*cosf(gamma))>0?1.0:-1.0;		//do flipphase
	}
	else if (mode==2){
		v=fabs(sqrtf(1.0-ampconst*ampconst)*sinf(gamma)+ampconst*cosf(gamma));		//	return abs ctf value
	}
	else{
		v=(sqrtf(1.0-ampconst*ampconst)*sinf(gamma)+ampconst*cosf(gamma));		//	return ctf value
	}
	
	return v;
}

//mode = 0 (default)  for template
//mode = 1 for raw image
__global__ void compute_area_sum_ofSQR(cufftComplex *data,float *res,int nx, int ny, int mode)
{
	extern __shared__ float sdata[];
    // each thread loads one element from global to shared mem
    int tid = threadIdx.x;
    long long  i = blockIdx.x*blockDim.x + threadIdx.x;
    int image_size = nx*ny;
    int local_id = i % image_size;
    int x = local_id % nx;
    int y = local_id / nx;
	int r = floor( hypotf(min(y,ny-y) ,min(x,nx-x)) + 0.5) - 1;
	int l = max(nx,ny);

	if (r < l/2 && r >= 0 && x<=l/2 && mode ==0) {
		sdata[tid] = data[i].x*data[i].x;
		sdata[tid+blockDim.x] = 1;
	}
	else if(mode == 1 && i<nx*ny)
	{
		sdata[tid] = data[i].x*data[i].x;
		sdata[tid+blockDim.x] = 1;
	}else
	{
		sdata[tid]=0;
		sdata[tid+blockDim.x] = 0;
	}
	__syncthreads();
	
	if (tid < 512) { sdata[tid] += sdata[tid + 512]; sdata[tid+blockDim.x] += sdata[tid +blockDim.x+ 512];} __syncthreads();
	if (tid < 256) { sdata[tid] += sdata[tid + 256]; sdata[tid+blockDim.x] += sdata[tid +blockDim.x+ 256];} __syncthreads();
	if (tid < 128) { sdata[tid] += sdata[tid + 128]; sdata[tid+blockDim.x] += sdata[tid +blockDim.x+ 128];} __syncthreads();
	if (tid < 64) { sdata[tid] += sdata[tid + 64]; sdata[tid+blockDim.x] += sdata[tid +blockDim.x+ 64];} __syncthreads();

	if(tid < 32)
	{
		sdata[tid] += sdata[tid + 32]; sdata[tid+blockDim.x] += sdata[tid +blockDim.x+ 32];
		sdata[tid] += sdata[tid + 16]; sdata[tid+blockDim.x] += sdata[tid +blockDim.x+ 16];
		sdata[tid] += sdata[tid + 8]; sdata[tid+blockDim.x] += sdata[tid +blockDim.x+ 8];
		sdata[tid] += sdata[tid + 4]; sdata[tid+blockDim.x] += sdata[tid +blockDim.x+ 4];
		sdata[tid] += sdata[tid + 2]; sdata[tid+blockDim.x] += sdata[tid +blockDim.x+ 2];
		sdata[tid] += sdata[tid + 1]; sdata[tid+blockDim.x] += sdata[tid +blockDim.x+ 1];
	}

	// write result for this block 
	if (tid == 0) {
		res[2*blockIdx.x] = sdata[0];
		res[2*blockIdx.x+1] = sdata[blockDim.x];
	}
}

__global__ void normalize(cufftComplex *data,int nx, int ny,float *means)
{
    // i <==> global ID
    long long  i = blockIdx.x*blockDim.x + threadIdx.x;
    int image_size = nx*ny;
	int template_id = i / image_size;

	if(means[template_id]!=0)	data[i].x=data[i].x/means[template_id];

	//ap2ri
	float tmp=data[i].x*sinf(data[i].y);
	data[i].x=data[i].x*cosf(data[i].y);
	data[i].y=tmp;
}

__global__ void rotate_IMG(float *d_image,float *d_rotated_image,float e,int nx,int ny)
{
	float cose=cos(e*PI/180);
	float sine=sin(e*PI/180);
	long long  id = blockIdx.x*blockDim.x + threadIdx.x;
	int j = id/nx, i = id%nx;
	float y = j-ny/2, x = i-nx/2;
	if(i>=nx || j>=ny) return;
	
	//Res of rotation from (x,y) 
	float res = 0;

	//(x,y) rotate e with (nx/2,ny/2) (clockwise) 
	float x2 = (cose*x+sine*y)+nx/2;
	float y2 = (-sine*x+cose*y)+ny/2;

	//Ouf of boundary after rotation
	if (x2<0||x2>nx-1.0||y2<0||y2>ny-1.0) res=0;
	else
	{
		int ii,jj;
		int k0,k1,k2,k3;
		float t,u,p0,p1,p2,p3;
		ii=floor(x2);
		jj=floor(y2);
		k0=ii+jj*nx;
		k1=k0+1;
		k2=k0+nx+1;
		k3=k0+nx;

		//handle situation when ii,jj are out of boundary
		if (ii==nx-1) { k1--; k2--; }
		if (jj==ny-1) { k2-=nx; k3-=nx; }
		t=(x2-(float)ii);
		u=(y2-(float)jj);
		float tt=1.0-t;
		float uu=1.0-u;

		//bilinear interpolation of raw data (i,j)(i+1,j)(i,j+1)(i+1,j+1)
		p0=d_image[k0]*tt*uu;
		p1=d_image[k1]*t*uu;
		p3=d_image[k3]*tt*u;
		p2=d_image[k2]*t*u;
		res=p0+p1+p2+p3;

	}

	// res <=> data[i+j*nx] after rotation
	d_rotated_image[id] = res;

}

__global__ void split_IMG(float *Ori,cufftComplex *IMG, int nx,int ny,int l,int bx,int overlap)
{
    long long  i = blockIdx.x*blockDim.x + threadIdx.x;
	int image_size = l*l;
	int image_id = i/image_size;
    int local_id = i % image_size;
    int x = local_id % l;
    int y = local_id / l;

	int tmp = l - overlap;

	int area_x_id = image_id%bx;
	int area_y_id = image_id/bx;
	int ori_x = area_x_id*tmp + x;
	int ori_y = area_y_id*tmp + y;

	if(ori_x>=nx || ori_y>=ny) return;
	IMG[i].x = Ori[ori_x + ori_y*nx];
}

//Tl = template(template has been predefined by C++)
__global__ void compute_corner_CCG(cufftComplex *CCG, cufftComplex *Tl, cufftComplex *IMG, int l, int block_id)
{
	//On this function,block means subimage splitted from IMG, not block ON GPU
	long long  i = blockIdx.x*blockDim.x + threadIdx.x;
	//Area of rectangle, l^2
	int l2 = l*l;

	//Local id corresponding to splitted IMG 
	int local_id = i%l2;
	int local_x = local_id%l;
	int local_y = local_id/l;

	int off = block_id * l2;

	//Global ID in IMG
	int j = local_x + local_y*l + off;

	//CCG[i] = IMG'[i]*template[i]
	// ' means conjugate
	CCG[i].x = (IMG[j].x*Tl[i].x+IMG[j].y*Tl[i].y);
	CCG[i].y = (IMG[j].y*Tl[i].x-IMG[j].x*Tl[i].y);

	//Move center to around
	int of = (l/2)%2,st;
	if( of == local_y%2) st = 1; else st = 0;
	if( (local_x-st)%2 == 0 ) 
	{
		CCG[i].x *= -1;
		CCG[i].y *= -1;
	}


}

//"MAX" reduction for *odata : return max{odata[i]},i
//"SUM" reduction for *odata : return sum{odata[i]},sum{odata[i]^2}
__global__ void get_peak_and_SUM(cufftComplex *odata,float *res,int l,float d_m,int x_bound,int y_bound)
{
	extern __shared__ float sdata[];
    // each thread loads one element from global to shared mem
    int tid = threadIdx.x;
    long long  i = blockIdx.x*blockDim.x + threadIdx.x;
	int image_size = l*l;
    int local_id = i % image_size;
    int x = local_id % l;
    int y = local_id / l;

	sdata[tid] = odata[i].x;
	if(x>=x_bound || y>=y_bound ||x<d_m/4 || x>l-d_m/4 || y<d_m/4 || y>l-d_m/4 ) sdata[tid] = 0;
	sdata[tid+blockDim.x] = local_id;
	sdata[tid+2*blockDim.x] = odata[i].x;
	sdata[tid+3*blockDim.x] = odata[i].x*odata[i].x;
	__syncthreads();

	for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
		if(tid<s)
		{
			//find max
			if(sdata[tid+s]>sdata[tid]){
				sdata[tid] = sdata[tid+s];
				sdata[tid+blockDim.x] = sdata[tid+blockDim.x+s];
			}
			//sum of data[i] & data[i]^2
			sdata[tid+2*blockDim.x] += sdata[tid+2*blockDim.x + s]; 
			sdata[tid+3*blockDim.x] += sdata[tid+3*blockDim.x + s];
		}
		__syncthreads();
	}
	if(tid==0){
		res[blockIdx.x*4] = sdata[0];
		res[blockIdx.x*4+1] = sdata[blockDim.x];
		res[blockIdx.x*4+2] = sdata[2*blockDim.x];
		res[blockIdx.x*4+3] = sdata[3*blockDim.x];
	}

}

// CUFFT will enlarge VALUE to N times. Restore it
__global__ void scale(cufftComplex *data,int l2)
{
	long long  i = blockIdx.x*blockDim.x + threadIdx.x;
	data[i].x /= l2;
	data[i].y /= l2;
}

__global__ void clear_image(cufftComplex *data)
{
	long long  i = blockIdx.x*blockDim.x + threadIdx.x;
	data[i].x = 0;
	data[i].y = 0;
}

__global__ void Complex2float(float *f, cufftComplex *c, int nx, int ny)
{
	long long  i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= nx*ny) return;
	f[i] = c[i].x;
}

__global__ void float2Complex(cufftComplex *c, float *f, int nx, int ny)
{
	long long  i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= nx*ny) return;
	c[i].x = f[i] ;
	c[i].y = 0 ;
}

__global__ void do_phase_flip(cufftComplex *filter, Parameters para, int nx, int ny)
{
    long long  i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= nx*ny) return;
    int x = i % nx;
    int y = i / nx;
	float v=CTF_AST(x,(y+ny/2)%ny,nx,ny,para.apix,para.ds,para.dfu,para.dfv,para.dfdiff,para.dfang,para.lambda,para.cs,para.ampconst,0);

	filter[i].x *= v;
	filter[i].y *= v;

}

void cudaMemoryTest()
{
    const unsigned int N = 1048576;
    const unsigned int bytes = N * sizeof(int);
    int *h_a = (int*)malloc(bytes);
    int *d_a;
    CUDA_CALL(cudaMalloc((int**)&d_a, bytes));

    memset(h_a, 0, bytes);
    CUDA_CALL(cudaMemcpy(d_a, h_a, bytes, cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(h_a, d_a, bytes, cudaMemcpyDeviceToHost));
    printf("Test finished.\n");
}