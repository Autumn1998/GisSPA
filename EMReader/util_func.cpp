#include "util_func.h"

// return a random number between lo and hi
float frand(float lo, float hi)
{
    static int f=1;

    if (f) { srandom(time(0)); f=0; }

    return ((float)random()/2147483647.0*(hi-lo)+lo);
}

float grand(float mean,float sigma)
{
    float x,y,r,f;

    do {
        x=frand(-1.0,1.0);
        y=frand(-1.0,1.0);
        r=x*x+y*y;
    } while (r>1.0||r==0);
    f=sqrt(-2.0*log(r)/r);

    return x*f*sigma+mean;
}

int Split(const char* str, vector<string>& results, const char* delim, bool empties)
{
	char* pstr = const_cast<char*>(str);
	char* r = NULL;
	if (delim) r = strstr(pstr, delim);
	else {
		r = NULL;
		for(int i=0; i<strlen(pstr); i++) {
			if(isspace(pstr[i])) {
				r = &pstr[i];
				break;
			}
		}
	}
	int dlen = 1;
	if (delim) dlen = int(strlen(delim));	// Cast used to avoid compiler warning

	while( r != NULL )
	{
		char* cp = new char[(r-pstr)+1];
		memcpy(cp, pstr, (r-pstr));
		cp[(r-pstr)] = '\0';

		if( strlen(cp) > 0 || empties )
		{
			string s(cp);
			results.push_back(s);
		}

		delete[] cp;
		pstr = r + dlen;
		if (delim) r = strstr(pstr, delim);
		else {
			r = NULL;
			for(int i=0; i<strlen(pstr); i++) {
				if(isspace(pstr[i])) {
					r = &pstr[i];
					break;
				}
			}
		}
	}

	if( strlen(pstr) > 0 || empties )
	{
		results.push_back(string(pstr));
	}

	return int(results.size()); // Cast used to avoid compiler warning
}

void printHelpMsg()
{
	printf("This program detect targets with orientations and tanslations.\n");
	printf("How to use this problem:  \"./main config_file\"  \n");
	printf("HDF5 and CUDA lib are needed to run\n");
	printf("All parameters should be set on the config file.\n");
	printf("All parameters are listed as below which can be set at config file.\n");
	printf("(for .eg,  \"input  =  /Data/inputfile\" )");
	printf("Requested Parameters:\n");
	printf("input            = input filtered images lstfile with ctf information\n");
	printf("template         = input 2D projections templates in .hdf format\n");
	printf("eulerfile        = euler file with euler values\n");
	printf("angpix           = input pixel size in angstroms\n");
	printf("phistep          = inplane rotation sampling\n");
	printf("kk               = overlapping density parameter, default is 3.\n");
	printf("energy           = accerlerating voltage in kV.\n");
	printf("cs               = spherical aberration in um.\n");
	printf("Highres          = high resolution cut \n");
	printf("Lowres           = low resolution cut\n");
	printf("diameter         = target diameter in pixels\n");
	printf("Optional Parameters:\n");
        printf("norm_type        = CCG norm type, 0/ 1\n");
        printf("invert           = 0 / 1, choose whether to invert contrast, 0 for skipping and 1 for invert\n");
	printf("threshold        = cc threshold value, only output score beyond this value\n");
	printf("output           = output lstfile filaname\n");
	printf("first            = the first image id to process.\n");
	printf("last             = the last image id to process.\n");
	printf("window_size      = the window size which is splitted from raw IMG.\n");
	printf("GPU_ID           = ID of GPU device.\n");
	printf("phase_flip       = Whether use phase filp operation(1) or not(0).\n");
	printf("overlap          = size of overlap between diff window.\n");
	
}

bool is_little_endian() {	// check if this machine is little endian
	int one=1;
	char* bytep=(char*)(&one);
	if(bytep[0]==1 && bytep[1]==0 && bytep[2]==0 && bytep[3]==0) return true;
	else	return false; 
}

bool is_big_endian() {		// check if this machine is big endian
	return ! is_little_endian();
}
