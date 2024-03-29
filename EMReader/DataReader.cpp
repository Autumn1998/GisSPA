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
		if(option.compare("--norm_type")==0){
		        i++;
		        p->norm_type=atoi(argv[i]);
		        return;
	       }
	        if(option.compare("--invert")==0){
	                i++;
	                p->invert=atoi(argv[i]);
	                return; 
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
		if(option.compare("--phase_flip")==0){
			i++;
			p->phase_flip=atoi(argv[i]);
			continue;
		}
		printf("Undefined option: %s .Abort.\n",argv[i]);
	}
	p->overlap =  p->d_m ;
	if(p->padding_size*p->padding_size % 1024 != 0 || (p->padding_size <= p->overlap))
	{
		printf("Padded IMG size must be N*1024 and larger than overlap{ = 1.5*d_m}\n");
		int tmp = p->padding_size;
		int ori = max(p->padding_size,p->overlap);
		for(int s = 128;s<2048*4;s+=32)
			if(s > ori) 
			{
				p->padding_size = s;
				break;
			}
		printf("Adjust Padded IMG size from %d to %d\n",tmp,p->padding_size);
	}
}

// Parse line to para
void parse(string line, Parameters *p)
{
	string option;
	char value[200];
	int i;
	// get the type of line
	for(i=0;i<line.size();i++) if(line[i] == ' ' || line[i] == '=') break;
	option = line.substr(0,i);
	for(;i<line.size();i++) if(line[i]!=' ' && line[i]!='=') break;
	strcpy(value, line.substr(i,line.size()-i).c_str());

	// parse value and option
	if(option.compare("input")==0){
		i++;
		strcpy(p->inlst, value);
		return;
	}
	if(option.compare("template")==0){
		i++;
		strcpy(p->temp2d, value);
		return;
	}		
	if(option.compare("eulerfile")==0){
		i++;
		strcpy(p->eulerf, value);
		return;
	}
	if(option.compare("kk")==0){
		i++;
		p->kk=atof(value);
		return;
	}
	if(option.compare("angpix")==0){
		i++;
		p->apix=atof(value);
		return;
	}
	if(option.compare("phistep")==0){
		i++;
		p->phi_step=atof(value);
		return;
	}
	if(option.compare("energy")==0){
		i++;
		p->energy=atof(value);
		return;
	}
	if(option.compare("cs")==0){
		i++;
		p->cs=atof(value);
		return;
	}
	if(option.compare("Highres")==0){
		i++;
		p->highres=atof(value);
		return;
	}
	if(option.compare("Lowres")==0){
		i++;
		p->lowres=atof(value);
		return;
	}
	if(option.compare("diameter")==0){
		i++;
		p->d_m=atof(value);
		return;
	}
	if(option.compare("output")==0){
		i++;
		strcpy(p->outlst, value);
		return;
	}
	if(option.compare("first")==0){
		i++;
		p->first=atoi(value);
		return;
	}
	if(option.compare("last")==0){
		i++;
		p->last=atoi(value);
		return;
	}
	if(option.compare("threshold")==0){
		i++;
		p->thres=atof(value);
		return;
	}
	if(option.compare("GPU_ID")==0){
		i++;
		p->device_id=atoi(value);
		return;
	}
	if(option.compare("window_size")==0){
		i++;
		p->padding_size=atoi(value);
		return;
	}
	if(option.compare("phase_flip")==0){
		i++;
		p->phase_flip=atoi(value);
		return;
	}
	if(option.compare("overlap")==0){
		i++;
		p->overlap=atoi(value);
		return;
	}
	if(option.compare("norm_type")==0){
		i++;
		p->norm_type=atoi(value);
		return;
	}
	if(option.compare("invert")==0){
	       i++;
	       p->invert=atoi(value);
	       return; 
        }
	printf("Undefined option: %s .Abort.\n",value);
}

// Parse config file
void readConfig(char *path, Parameters *p)
{
	std::ifstream conf(path);
	if (!conf)
	{
		printf("Error => Open config file failed : %s\n\n",path);
		return;
	}
    std::string buf;
	while(!conf.eof()){	
	    getline(conf,buf,'\n');
		if(buf[0] == '#'|| buf == "") continue;
		parse(buf,p);
	}
	
	checkRequestPara(p);
} 

void readEulerData(char *eulerf, EulerData *euler)
{
	std::ifstream eulerfile(eulerf);
	if (!eulerfile)
	{
		printf("Error => Open euler file failed : %s\n\n",eulerf);
		return;
	}
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
	if (!snr_parm)
	{
		printf("Error => Open snr file failed : %s\n\n",snr);
		return;
	}
	string buf;
	while(!snr_parm.eof()){
	    getline(snr_parm, buf,'\n');
		if(buf[0] == '#'|| buf == "") continue;
		sscanf(buf.c_str(),"%f\t%f\t%f\t%f\t%f\t%f",a,b,b2,bfactor,bfactor2,bfactor3);
	}
}

void checkRequestPara(Parameters *p)
{
	if(p->inlst[0] == ' ') printf("Error : lst/star file is requested.\n"); 
	if(p->temp2d[0] == ' ') printf("Error : template file is requested.\n");
	if(p->eulerf[0] == ' ') printf("Error : euler file is requested.\n");
	if(p->apix < 0 ) printf("Error : angpix is requested.\n");
	if(p->phi_step < 0 ) printf("Error : phistep is requested.\n");
	if(p->kk < 0 ) printf("Error : kk is requested.\n");
	if(p->energy < 0 ) printf("Error : energy is requested.\n");
	if(p->cs < 0 ) printf("Error : cs is requested.\n");
	if(p->highres < 0 ) printf("Error : highres is requested.\n");
	if(p->lowres < 0 ) printf("Error : lowres is requested.\n");
	if(p->d_m < 0 ) printf("Error : diameter is requested.\n");
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

void parsePairs( vector<string> pairs, float *defocus, float *dfdiff, float *dfang)
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
