#include "emhdf_.h"
#include "emdata_.h"

hid_t EMhdf::euler_type = -1;
hid_t EMhdf::mapinfo_type = -1;

void EMhdf::create_enum_types()
{
    static int enum_types_created = 0;
    
    if (!enum_types_created) {

        Euler::EulerType e;
        euler_type = H5Tcreate(H5T_ENUM, sizeof(Euler::EulerType));
        H5Tenum_insert(euler_type, "EMAN", (e=Euler::EMAN, &e));
        H5Tenum_insert(euler_type, "IMAGIC", (e=Euler::IMAGIC, &e));
        H5Tenum_insert(euler_type, "SPIN", (e=Euler::SPIN, &e));
        H5Tenum_insert(euler_type, "QUAT", (e=Euler::QUAT, &e));
        H5Tenum_insert(euler_type, "MATRIX", (e=Euler::MATRIX, &e));
        H5Tenum_insert(euler_type, "SGIROT", (e=Euler::SGIROT, &e));
        H5Tenum_insert(euler_type, "SPIDER", (e=Euler::SPIDER, &e));
        H5Tenum_insert(euler_type, "MRC", (e=Euler::MRC, &e));

        MapInfoType m;
        mapinfo_type = H5Tcreate(H5T_ENUM, sizeof(MapInfoType));
        
        H5Tenum_insert(mapinfo_type, "NORMAL", (m=NORMAL, &m));
        H5Tenum_insert(mapinfo_type, "ICOS2F_FIRST_OCTANT", (m=ICOS2F_FIRST_OCTANT, &m));
        H5Tenum_insert(mapinfo_type, "ICOS2F_FULL", (m=ICOS2F_FULL, &m));
        H5Tenum_insert(mapinfo_type, "ICOS2F_HALF", (m=ICOS2F_HALF, &m));
        H5Tenum_insert(mapinfo_type, "ICOS3F_HALF", (m=ICOS3F_HALF, &m));
        H5Tenum_insert(mapinfo_type, "ICOS3F_FULL", (m=ICOS3F_FULL, &m));
        H5Tenum_insert(mapinfo_type, "ICOS5F_HALF", (m=ICOS5F_HALF, &m));
        H5Tenum_insert(mapinfo_type, "ICOS5F_FULL", (m=ICOS5F_FULL, &m));

        enum_types_created = 1;
    }
}

EMhdf::~EMhdf()
{
    close_dataset(cur_dataset);
    H5Gclose(group);
    H5Fclose(file);
}

int EMhdf::is_valid(const char* filename)
{
    int valid = 0;
    
    if (filename) {
        FILE *in=fopen(filename,"r");
        if (in) {
            char x[8];
            char y[8] = { (char)137,72,68,70,13,10,26,10 };
            int realRead = fread(x,8,1,in);
            if(realRead == 0) valid = 0;
            fclose(in);
            if (strncmp(x,y,8)==0) valid=1;
        }
    }
    return valid;
}

void EMhdf::close_dataset(hid_t dataset)
{
    hdf_err_off();
    if (dataset >= 0) {
    	H5Dclose(dataset);
    }
    hdf_err_on();
}

herr_t file_info(hid_t loc_id, const char *name, void *opdata)
{
    loc_id = loc_id;
    vector<int>* dataset_ids = static_cast<vector<int>*>(opdata);

    const char* magic = EMhdf::get_item_name(EMhdf::COMPOUND_DATA_MAGIC);
    if (strncmp(name, magic, strlen(magic)) != 0) {
        int id = atoi(name);
        dataset_ids->push_back(id);
    }
    
    return 0;
 }

 const char* EMhdf::get_item_name(Nametype type)
{
    switch (type) {
        case ROOT_GROUP:
        return "/";
        case CTFIT:
        return "ctfit";
        case NUMDATASET:
        return "num_dataset";
        case COMPOUND_DATA_MAGIC:
        return "compound.";
    }
    return "unknown";	
}

int EMhdf::init_test(const char * hdf_filename)
{
	H5Eset_auto1(0, 0);	// Turn off console error logging.
	
	hid_t fileid = H5Fopen(hdf_filename, H5F_ACC_RDWR, H5Pcreate(H5P_FILE_ACCESS));
	hid_t groupid = H5Gopen1(fileid, "/");
	hid_t attid = H5Aopen_name(groupid, "num_dataset");
	
	if (attid < 0) {
		H5Gclose(groupid);
		H5Fclose(fileid);
		return -1;
	}
	else {
		H5Aclose(attid);
		H5Gclose(groupid);
		H5Fclose(fileid);
		return 0;
	}
}

int* EMhdf::read_dims(int dataset_id, int* p_ndim)
{
    set_dataset(dataset_id);

    if (cur_dataset < 0) {
        char cur_dataset_name[32];
        sprintf(cur_dataset_name, "%d", dataset_id);
        cur_dataset = H5Dopen1(file, cur_dataset_name);
        
        if (cur_dataset < 0) {
            fprintf(stderr, "Error in reading data dimensions in hdf: %d\n", dataset_id);
            return 0;
        }
    }
    
    hid_t dataspace = H5Dget_space(cur_dataset);
    int rank = H5Sget_simple_extent_ndims(dataspace);
    hsize_t* dims = new hsize_t[rank];
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    
    int* dims1 = new int[rank];
    for (int i = 0; i < rank; i++) {
	    dims1[i] = dims[i];	
    }
    
    
    H5Sclose(dataspace);
    delete [] dims;
    dims = 0;

    (*p_ndim) = rank;
    return dims1;    
}

void EMhdf::hdf_err_off()
{
    H5Eget_auto1(&old_func, &old_client_data);
    H5Eset_auto1(0, 0);
}

void EMhdf::hdf_err_on()
{
    H5Eset_auto1(old_func, old_client_data);
}


void EMhdf::set_dataset(int dataset_id)
{
    int need_update = 0;
    
    if (cur_dataset_id >= 0) {
        if (dataset_id != cur_dataset_id) {
            need_update = 1;
        }
    }
    else {
        cur_dataset_id = dataset_id;
        need_update = 1;
    }

    if (need_update) {
        char cur_dataset_name[32];
        sprintf(cur_dataset_name, "%d", dataset_id);
        hdf_err_off();
        close_dataset(cur_dataset);
        cur_dataset = H5Dopen1(file, cur_dataset_name);
        hdf_err_on();
    }
}


EMhdf::EMhdf(const char* filename, const char* rw_mode)
{
    assert(filename);
    assert(rw_mode);
    
    if (strcmp(rw_mode, "r") == 0) {
	    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    }
    else if (strcmp(rw_mode, "rw") == 0) {
        hdf_err_off();
        file =  H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
        hdf_err_on();
        if (file < 0) {
            file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            H5Fclose(file);
            file =  H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
        }
    }
    else if (strcmp(rw_mode, "w") == 0) {
	    file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    else {
    	file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    }
    if(file<0) {
		fprintf(stderr,"EMhdf: cannot open file %s\n", filename);	
	}
    
    group = H5Gopen1(file, get_item_name(ROOT_GROUP));
    cur_dataset = -1;
    cur_dataset_id = -1;
    
    old_func = 0;
    old_client_data = 0;
    H5Giterate(file, get_item_name(ROOT_GROUP), NULL, file_info, &dataset_ids);
    create_enum_types();
}

int EMhdf::read_data(int dataset_id, float** p_data)
{
    set_dataset(dataset_id);

    float* data = *p_data;
    
    if (cur_dataset < 0) {
        char cur_dataset_name[32];
        sprintf(cur_dataset_name, "%d", dataset_id);
        cur_dataset = H5Dopen1(file, cur_dataset_name);

        if (cur_dataset < 0) {
            fprintf(stderr, "Error: hdf file has no data with id = %d\n", dataset_id);
            return 1;
        }
    }
    
    hid_t datatype = H5Dget_type(cur_dataset);
    H5T_class_t t_class = H5Tget_class(datatype);

    if (t_class != H5T_FLOAT) {
        fprintf(stderr, "Error: can only read/write hdf FLOAT data\n");
        H5Tclose(datatype);
        return 1;
    }
    
    int err = H5Dread(cur_dataset, H5T_NATIVE_FLOAT,  H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (err < 0) {
        fprintf(stderr, "Error: read hdf data failed with id = %d\n", dataset_id);
        return 1;
    }

    H5Tclose(datatype);    
    return 0;
}


int EMhdf::read_int_attr(int dataset_id, const char* attr_name)
{
    set_dataset(dataset_id);

    int value = 0;
    hid_t attr = H5Aopen_name(cur_dataset, attr_name);
    H5Aread(attr, H5T_NATIVE_INT, &value);
    H5Aclose(attr);
    return value;
}

 
float EMhdf::read_float_attr(int dataset_id, const char* attr_name)
{
    set_dataset(dataset_id);
    return read_float_attr(attr_name);
}



float EMhdf::read_float_attr(const char* attr_name)
{
    float value = 0;
    hid_t attr = H5Aopen_name(cur_dataset, attr_name);
    H5Aread(attr, H5T_NATIVE_FLOAT, &value);    
    H5Aclose(attr);
    return value;
}



const char* EMhdf::read_string_attr(int dataset_id, const char* attr_name)
{
    set_dataset(dataset_id);

    char* tmp_value = new char[MAXPATHLEN];
    hid_t attr = H5Aopen_name(cur_dataset, attr_name);
    hid_t datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size(datatype, MAXPATHLEN);
    H5Aread(attr, datatype, tmp_value);

    H5Tclose(datatype);
    H5Aclose(attr);

    char* value = new char[strlen(tmp_value) + 1];
    strcpy(value, tmp_value);
    delete [] tmp_value;
    tmp_value = 0;
    
    return value;
}

void EMhdf::read_array_attr(int dataset_id, const char* attr_name, void* value)
{
    set_dataset(dataset_id);

    hid_t attr = H5Aopen_name(cur_dataset, attr_name);
    hid_t datatype = H5Aget_type(attr);    
    H5Aread(attr, datatype, value);    
    H5Tclose(datatype);
    H5Aclose(attr);
}


int EMhdf::read_euler_attr(int dataset_id, const char* attr_name)
{
    set_dataset(dataset_id);

    Euler::EulerType val;
    hid_t attr = H5Aopen_name(cur_dataset, attr_name);
    H5Aread(attr, euler_type, &val);
    H5Aclose(attr);
    return static_cast<int>(val);
}
   

int EMhdf::read_ctfit(int dataset_id, CTFParam* value)
{
    assert(value);
    hid_t cur_dataset_orig = cur_dataset;

    const char* cur_dataset_name = get_compound_name(dataset_id, get_item_name(CTFIT));
    cur_dataset = H5Dopen1(file, cur_dataset_name);
    int err = 0;
    if (cur_dataset < 0) {
	    err = 1;
    }
    else {
        value->defocus = read_float_attr("defocus");
        value->bfactor = read_float_attr("bfactor");
        value->amplitude = read_float_attr("amplitude");
        value->ampcont = read_float_attr("ampcont");
        value->noise1 = read_float_attr("noise1");
        value->noise2 = read_float_attr("noise2");
        value->noise3 = read_float_attr("noise3");
        value->noise4 = read_float_attr("noise4");
        value->voltage = read_float_attr("voltage");
        value->cs = read_float_attr("cs");
        value->apix = read_float_attr("apix");
    }

    delete []  cur_dataset_name;
    cur_dataset_name = 0;
    H5Dclose(cur_dataset);
    cur_dataset = cur_dataset_orig;
    return err;
}

const char* EMhdf::get_compound_name(int id, const char* name)
{
    const char* magic = get_item_name(COMPOUND_DATA_MAGIC);
    char* compound_name = new char[strlen(name) + strlen(magic) + 32];
    sprintf(compound_name, "%s%d.%s", magic, id, name);
    return compound_name;
}

int EMhdf::read_mapinfo_attr(int dataset_id, const char* attr_name)
{
    set_dataset(dataset_id);

    MapInfoType val;
    hid_t attr = H5Aopen_name(cur_dataset, attr_name);
    H5Aread(attr, mapinfo_type, &val);
    H5Aclose(attr);
    return static_cast<int>(val);
}

