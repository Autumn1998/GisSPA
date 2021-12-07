#ifndef _emhdf2_h_
#define _emhdf2_h_


#define H5_USE_16_API

#include <hdf5.h>
#include <cstring>

#include <iostream>
#include <map>
#include "util_func.h"
#include "emdata_.h"
class emdataHeader;
class Euler;
class EMObject;
using std::string;

class EMhdf2 {
public:
	EMhdf2(const char* hdf_filename, const char* rw_mode="rw");
    ~EMhdf2();
    
    /** test whether this hdf file is new HDF5 or old HDF5 implementation
     * 
     * After you make change to this class, please check the HDF5 file created 
	 * by EMAN2 with the h5check program from:
	 * ftp:://ftp.hdfgroup.org/HDF5/special_tools/h5check/
	 * to verify the HDF5 file is compliant with the HDF5 File Format Specification.
     * 
     * @param hdf_filename file name
     * @return 0 for new HDF5 file, -1 for old-style HDF5 file
     * */
    static int init_test(const char * hdf_filename);
    
    emdataHeader * read_header(int image_index, Euler& euler);
    
    int read_data(float* data, int image_index = 0);
    
    int erase_header(int image_index);
    
    int write_header(const emdataHeader & header, int image_index, const Euler& euler);
    
    int write_data(const float* const data, int image_index = 0);
    
    /** get the number of images in this file*/
    int get_num_dataset();
    
private:
	hid_t file;
	hid_t group;
	hid_t accprop;
	hid_t simple_space;
	string filename;   
	bool initialized;
	const char * rwmode;
	
	int write_attr(hid_t loc, const char *name, EMObject obj);
	EMObject read_attr(hid_t attr);
	
	void init();
};

#endif	//_emhdf2_h_