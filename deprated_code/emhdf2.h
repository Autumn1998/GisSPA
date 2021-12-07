#ifndef _emhdf2_h_
#define _emhdf2_h_
#include "emdata_.h"
#include <hdf5.h>
#include <string>
#include <map>

using namespace std;

class emdataHeader;
class Euler;
class EMObject;

class EMhdf2 {
    public:
        EMhdf2(const char* hdf_filename, const char* rw_mode="rw");
        ~EMhdf2();
        emdataHeader * read_header(int image_index, Euler& euler);
        int read_data(float* data, int image_index);
                
    private:
        hid_t file;
        hid_t group;
        hid_t accprop;
        hid_t simple_space;
        string filename;   
        bool initialized;
        const char * rwmode;
        
        EMObject read_attr(hid_t attr);
        int write_attr(hid_t loc,const char *name,EMObject obj);
        
        void init();
};

#endif