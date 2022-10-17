#ifndef DATAREADER_H
#define DATAREADER_H
#include <string>
#include <vector>
#include <cstdio>
#include "util_func.h"
#include "emdata_.h"
#include <sstream>   // string streams
#include <fstream>
#include <vector>
#define MAXPATHLEN 1024

//**** Struct contains all possible parameters in argv[]
struct Parameters
{	
	char inlst[512]="  ",temp2d[512]="  ",eulerf[512]="  ";
	float apix=-1,phi_step=-1;
	float kk=-1,energy=-1,cs=-1;
	float highres=-1,lowres=-1;
	float d_m=-1;
	float edge_half_width=4;
	float lambda,dfu,dfv,ds;
	float ampconst=0.07;
	float defocus, dfdiff, dfang;
	float a=-9.32, b=2.65, b2=0.01908, bfactor=-78.7757, bfactor2=-12.9121, bfactor3=1.28732;
	int template_x = 0;
	int template_y = 0;
	int template_z = 0;
	int block_x = 0;
	int block_y = 0;
	int norm_type = 0;

	//optional para
	int first=-1,last=-1;
	char outlst[512] = "./outlst";
	float thres = 7;
	int padding_size=320;
	int phase_flip = 0;
	int device_id = 0;
	//will be d_m after inputting of d_M if user do not make it fixed
	int overlap = -1;
	int replace_hot = 0; 
	int maximum_n_sigma = 2;

	void printAllPara()
	{
		printf("Parameter requested:\n");
		printf("--input            %s\n",inlst);
		printf("--template         %s\n",temp2d);
		printf("--eulerfile        %s\n",eulerf);
		printf("--angpix           %f\n",apix);
		printf("--phistep          %f\n",phi_step);
		printf("--kk               %f\n",kk);
		printf("--energy           %f\n",energy);
		printf("--cs               %f\n",cs);
		printf("--Highres          %f\n",highres);
		printf("--Lowres           %f\n",lowres);
		printf("--diameter         %f\n",d_m);

		printf("Parameter optional:\n");
		printf("--threshold        %f\n",thres);
		printf("--output           %s\n",outlst);
		printf("--first            %d\n",first);
		printf("--last             %d\n",last);
		printf("--window_size      %d\n",padding_size);
		printf("--GPU_ID           %d\n",device_id);
		printf("--phase_flip       %d\n",phase_flip);
		printf("--overlap          %d\n",overlap);
		printf("--norm_type        %d\n",norm_type);
		printf("--replace_hot      %d\n",replace_hot);
		printf("--max_n_sigma      %d\n",maximum_n_sigma);

		printf("Parameter fixed:\n");
		printf("--a                %f\n",a);
		printf("--b                %f\n",b);
		printf("--b2               %f\n",b2);
		printf("--bfactor          %f\n",bfactor);
		printf("--bfactor2         %f\n",bfactor2);
		printf("--bfactor3         %f\n\n",bfactor3);
	}
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

//**** Parse from argv
void readParameters(int argc, char *argv[], Parameters *p);

//**** Parse from config file
void readConfig(char *confil_path, Parameters *p);

// Check all requested para
void checkRequestPara(Parameters *P);

void readEulerData(char *eulerf, EulerData *euler);

void readSNRWeight(char *snr,float* a,float* b,float* b2,float* bfactor,float* bfactor2,float* bfactor3);

int readInLst_and_consturctPairs(char * inlst, char *t, vector<string> *pairs, int *nn, int n);

void parsePairs( vector<string> pairs, float *defocus, float *dfdiff, float *dfang);

#endif
