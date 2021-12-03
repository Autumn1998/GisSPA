#ifndef _emhdf_h_
#define _emhdf_h_
#undef H5Gopen_vers
#define H5Gopen_vers 1
#include <hdf5.h>
#include <vector>
#include <string>
#include <cassert>
#include <stdlib.h>
#include <stdio.h>
#include "emdata_.h"

using namespace std;
class EMObject;

class EMhdf {
    public:
        enum DataType { INT, FLOAT, STRING };
        enum Nametype { ROOT_GROUP, CTFIT, NUMDATASET, COMPOUND_DATA_MAGIC };
        
    public:
        EMhdf(const char* hdf_filename, const char* rw_mode="rw");
        ~EMhdf();

        /** test whether this hdf file is new HDF5 or old HDF5 implementation
        * return 0 for new HDF5 file, -1 for old-style HDF5 file
        * */
        static int init_test(const char * hdf_filename);

        static int is_valid(const char* filename);

        int* read_dims(int dataset_id, int* p_ndim);
        int read_data(int dataset_id, float** p_data);
        int read_int_attr(int dataset_id, const char* attr_name);
        float read_float_attr(int dataset_id, const char* attr_name);
        const char* read_string_attr(int dataset_id, const char* attr_name);
        void read_array_attr(int dataset_id, const char* attr_name, void* value);
        int read_euler_attr(int dataset_id, const char* attr_name);
        int read_ctfit(int dataset_id, CTFParam* value);
        int read_mapinfo_attr(int dataset_id, const char* attr_name);

        void set_dataset(int dataset_id);
        static void create_enum_types();
        
        static const char* get_item_name(Nametype type);     
    
    private:
        hid_t file;
        hid_t group;
        hid_t cur_dataset;
        int cur_dataset_id;
        
        herr_t (*old_func)(void*);
        void *old_client_data;

        static hid_t euler_type;
        static hid_t mapinfo_type;

        vector<int> dataset_ids;
        
        float read_float_attr(const char* attr_name);
        void close_dataset(hid_t dataset);  
        const char* get_compound_name(int id, const char* name);
        
        void hdf_err_off();
        void hdf_err_on();

};



#endif