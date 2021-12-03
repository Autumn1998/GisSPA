#include "DataReader.h"

//**** Parse info from argv[]
void readParameters(int argc, char *argv[], Parameters *p)
{
	string option;
	for(int i=1;i<argc;i++){
		option=argv[i];
		if(option.compare("--input")==0){
			i++;
			strcpy(p->inlst, argv[i]);
			continue;
		}
		option=argv[i];
		if(option.compare("--output")==0){
			i++;
			strcpy(p->outlst, argv[i]);
			continue;
		}
                if(option.compare("--template")==0){
			i++;
			strcpy(p->temp2d, argv[i]);
			continue;
		}		
		if(option.compare("--eulerfile")==0){
			i++;
			strcpy(p->eulerf, argv[i]);
			continue;
		}
		if(option.compare("--weight")==0){
			i++;
			strcpy(p->snr, argv[i]);
			continue;
		}
		if(option.compare("--kk")==0){
			i++;
			p->kk=atof(argv[i]);
			continue;
		}
		if(option.compare("--angpix")==0){
			i++;
			p->apix=atof(argv[i]);
			continue;
		}
		if(option.compare("--phistep")==0){
			i++;
			p->phi_step=atof(argv[i]);
			continue;
		}
		if(option.compare("--first")==0){
			i++;
			p->first=atoi(argv[i]);
			continue;
		}
		if(option.compare("--last")==0){
			i++;
			p->last=atoi(argv[i]);
			continue;
		}
		if(option.compare("--energy")==0){
			i++;
			p->energy=atof(argv[i]);
			continue;
		}
		if(option.compare("--cs")==0){
			i++;
			p->cs=atof(argv[i]);
			continue;
		}
		if(option.compare("--Highres")==0){
			i++;
			p->highres=atof(argv[i]);
			continue;
		}
		if(option.compare("--Lowres")==0){
			i++;
			p->lowres=atof(argv[i]);
			continue;
		}
		if(option.compare("--diameter")==0){
			i++;
			p->d_m=atof(argv[i]);
			continue;
		}
		if(option.compare("--threshold")==0){
			i++;
			p->thres=atof(argv[i]);
			continue;
		}
		if(option.compare("--device")==0){
			i++;
			p->device_id=atoi(argv[i]);
			continue;
		}
		if(option.compare("--padding_size")==0){
			i++;
			p->padding_size=atoi(argv[i]);
			continue;
		}
		printf("Undefined option: %s .Abort.\n",argv[i]);
	}
	if(p->padding_size*p->padding_size % 1024 != 0)
	{
		printf("Padded IMG size must be N*1024\n");
		int ori = p->padding_size;
		for(int s = 128;s<2048;s+=32)
			if(s > ori) 
			{
				p->padding_size = s;
				break;
			}
		printf("Adjust Padded IMG size from %d to %d\n",ori,p->padding_size);
	}
}


void readEulerData(char *eulerf, EulerData *euler)
{
	std::ifstream eulerfile(eulerf);
    std::string buf;
    int cnt=0;
	while(!eulerfile.eof()){	
	    getline(eulerfile,buf,'\n');
		if(buf[0] == '#'|| buf == "") continue;
		sscanf(buf.c_str(),"%f %f",&euler->euler1[cnt],&euler->euler2[cnt]);
		cnt++;
	}
	euler->length = cnt;
}

void readSNRWeight(char *snr,float* a,float* b,float* b2,float* bfactor,float* bfactor2,float* bfactor3)
{
	std::ifstream snr_parm(snr);
	string buf;
	while(!snr_parm.eof()){
	    getline(snr_parm, buf,'\n');
		if(buf[0] == '#'|| buf == "") continue;
		sscanf(buf.c_str(),"%f\t%f\t%f\t%f\t%f\t%f",a,b,b2,bfactor,bfactor2,bfactor3);
	}
}

int readInLst_and_consturctPairs(char * inlst, char *t, vector<string> *pairs, int *nn, int n)
{
//***************************************************************
// Read from inlst. consturct pairs
//***************************************************************
	char buf[MAXPATHLEN+10],u[256];
	char *tmp; // To prevent the warning
	FILE *in = fopen(inlst,"rb");
	if (!in) { printf("'%s' not found\n",inlst); return (-3); }
	tmp = fgets(buf,MAXPATHLEN+10,in);
    if(tmp == NULL) return 0;
	if (strncmp(buf,"#LSX",4)==0) {
		int step=0;
		tmp = fgets(buf,MAXPATHLEN+10,in);
		if(tmp == NULL) return 0;
		tmp = fgets(buf,MAXPATHLEN+10,in);
		if(tmp == NULL) return 0;
		sscanf(buf+1," %d",&step);
		if (step==0) {fclose(in); return -1; }
		int tl=ftell(in);
		if (n>0) fseek(in,step*n,SEEK_CUR);
		if (!fgets(buf,MAXPATHLEN+10,in)) { fclose(in); return -2; }
		if ((ftell(in)-tl)%step !=0) {
			printf("Error in #LSX file detected. Please rerun lstfast.py\n");
			return(-3);
		}
		u[0]='\0';
		sscanf(buf," %d %s %[^\n]",nn,t,u);
		fclose(in);
	}
	else 
	{
		if (n<0) n=0;
		t[0]=0;
		int i;
		for (i=0; i<=n; i++) {
			u[0]='\0';
			if (!fgets(buf,MAXPATHLEN+10,in)) break;
			if (buf[0]=='#') i--;
			else sscanf(buf," %d %s %[^\n]",nn,t,u);
		}
		fclose(in);
		if (i<=n) return (-1);
	}
	Split(u,*pairs,0,0);
	return 1;
}

void parsePairs(char * t, int nn, vector<string> pairs, float *defocus, float *dfdiff, float *dfang)
{

//***************************************************************
// Get value from consturcted pairs
//***************************************************************

	for(int i=0; i<pairs.size(); i++){
		//此处直接使用sscanf写入%f,可能会损失精度， 发现91.5999998写入为91.6
		if(strncmp(pairs[i].c_str(),"defocus=",8)==0){
			if(sscanf(pairs[i].substr(8).c_str(),"%f",defocus)==1){
				*defocus = -*defocus;
			}
		}
		else if(strncmp(pairs[i].c_str(),"dfdiff=",7)==0){
			sscanf(pairs[i].substr(7).c_str(),"%f",dfdiff);
		}
		else if(strncmp(pairs[i].c_str(),"dfang=",6)==0){
			sscanf(pairs[i].substr(6).c_str(),"%f",dfang);
		}
	}
    
}
