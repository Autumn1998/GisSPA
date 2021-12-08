#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <string>
#include <string.h>
using namespace std;

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sys/types.h>

#if defined(__sparc)
#include <ieeefp.h>
#endif

#define PI (3.141592654)
#define RTOD (57.29577951)


#define PEAKMAX 2

//#define MSB_FIRST

#define ERR_EXIT	0
#define ERR_CRIT	1
#define ERR_ERROR	2
#define ERR_WARN	3
#define ERR_MSG		4
#define ERR_INFO	5
#define ERR_DEBUG	6

#define DBUG printf

#ifndef MAXPATHLEN
#define MAXPATHLEN	1024
#endif


#ifdef AIX
// This is basically finite(), but shouldn't crash on an alpha
// if exponent is 255, number is infinte or nan
inline int goodf(float *f) { 
    // the first is abnormal zero the second is +-inf or NaN
    if (finite(*f)) return 1;
    return 0;
    }
#else
// This is basically finite(), but shouldn't crash on an alpha
// if exponent is 255, number is infinte or nan
inline int goodf(float *f) { 

	//This is the old way to judge a good float, which cause problem on 
	//Fedora Core 64 bit system. Now use isinff() and isnanf() on Linux.
#if defined(_WIN32) || defined(__APPLE__) || defined(__sgi) || defined(aix) || defined(_AIX)
	// the first is abnormal zero the second is +-inf or NaN
	if ((((int *)f)[0] & 0x7f800000)==0 || (((int *)f)[0] & 0x7f800000)==255) return 0;
#else
	if(isinff(*f) || isnanf(*f)) return 0;
#endif	//_WIN32 || __APPLE__ || IRIX
	return 1;
}
#endif	//AIX

// return a random number between lo and hi
float frand(float lo, float hi);
float grand(float mean,float sigma);

int Split(const char* str, vector<string>& results, const char* delim, bool empties);
void printHelpMsg();

bool is_big_endian();
bool is_little_endian();
#endif