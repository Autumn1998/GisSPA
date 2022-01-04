#ifndef DATAREADER_H
#define DATAREADER_H
#include <string>
#include <vector>
#include <cstdio>
#include "util_func.h"
#include "emdata_.h"
#include <sstream>   // string streams
#include <fstream>
#define MAXPATHLEN 1024

//**** Struct contains all possible parameters in argv[]
struct Parameters
{	
	char inlst[512],outlst[512],temp2d[512],eulerf[512],snr[512],sigma_path[512];
	
	float apix,phi_step,first,last;
	float kk=3,energy=300,cs=2.7;

	float highres=8,lowres=30,thres;
	float d_m;
	float edge_half_width=4;
	float lambda,dfu,dfv,ds;
	float ampconst=0.07;
	float defocus, dfdiff, dfang;
	float a,b,b2,bfactor,bfactor2,bfactor3;
	int device_id = 6;
	int padding_size=256;
	int overlap = 0;
	int template_x = 0;
	int template_y = 0;
	int template_z = 0;
	int phase_flip = 0;
	float ds_Img = 0;

	int block_x = 0;
	int block_y = 0;

#ifdef DEBUG
	void printAllPara()
	{
		printf("--input            %s\n",inlst);
		printf("--output           %s\n",outlst);
		printf("--template         %s\n",temp2d);
		printf("--eulerfile        %s\n",eulerf);
		printf("--weight           %s\n",snr);
		printf("--angpix           %f\n",apix);
		printf("--phistep          %f\n",phi_step);
		printf("--first            %f\n",first);
		printf("--last             %f\n",last);
		printf("--kk               %f\n",kk);
		printf("--energy           %f\n",energy);
		printf("--cs               %f\n",cs);
		printf("--Highres          %f\n",highres);
		printf("--Lowres           %f\n",lowres);
		printf("--threshold        %f\n",thres);
		printf("--diameter         %f\n",d_m);
		printf("--padding_size     %d\n",padding_size);
		printf("--device_id        %d\n",device_id);
		printf("--phase_flip       %d\n",phase_flip);
		printf("--overlap          %d\n",overlap);
		//printf("--sigma file       %s\n",sigma_path);
	}
#endif
};

struct EulerData
{
	float euler1[20000];
	float euler2[20000];
	int length = 0;
};
//num of blocks in x,y axis
struct IMG
{
	int blockx,blocky;
};

//**** Parse info from argv[]
void readParameters(int argc, char *argv[], Parameters *p);

void readEulerData(char *eulerf, EulerData *euler);

void readSNRWeight(char *snr,float* a,float* b,float* b2,float* bfactor,float* bfactor2,float* bfactor3);

int readInLst_and_consturctPairs(char * inlst, char *t, vector<string> *pairs, int *nn, int n);

void parsePairs( vector<string> pairs, float *defocus, float *dfdiff, float *dfang);

#endif