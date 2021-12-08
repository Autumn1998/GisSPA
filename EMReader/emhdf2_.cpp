#include "emhdf2_.h"

using std::map;

int EMhdf2::init_test(const char * hdf_filename)
{
	H5Eset_auto(0, 0);	// Turn off console error logging.
	
	hid_t fileid = H5Fopen(hdf_filename, H5F_ACC_RDWR, H5Pcreate(H5P_FILE_ACCESS));
	hid_t groupid = H5Gopen(fileid, "/");
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

void EMhdf2::init()
{
	if (initialized) {
		return;
	}
	
	H5Eset_auto(0, 0);	// Turn off console error logging.
    
    if (strcmp(rwmode, "r") == 0) {
		file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, accprop);
    }
    else if (strcmp(rwmode, "rw") == 0) {
    	file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, accprop);
    	if (file < 0) {
    		file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, accprop);
    	}
    }
    
    if(file<0) {	
		fprintf(stderr,"EMhdf: cannot open file %s\n", filename.c_str());	
	}
	
	group = H5Gopen(file,"/MDF/images");
	if (group < 0) {
		if (strcmp(rwmode, "r") == 0) {
			fprintf(stderr, "HDF5 file has no image data (no /TEM group)\n");
		}
		group = H5Gcreate(file, "/MDF", 64);		// create the group for Macromolecular data
		if (group<0) {
			fprintf(stderr, "Unable to add image group (/MDF) to HDF5 file\n");
		}
		H5Gclose(group);
		group=H5Gcreate(file,"/MDF/images",4096);		// create the group for images/volumes
		if (group<0) {
			fprintf(stderr, "Unable to add image group (/MDF/images) to HDF5 file\n");
		}
		write_attr(group,"imageid_max", EMObject(-1));
	}
    
    initialized = true;
}

EMhdf2::EMhdf2(const char* hdf_filename, const char* rw_mode) 
	: file(-1), group(-1), simple_space(0),
		filename(hdf_filename), initialized(false), rwmode(rw_mode)
{
	accprop = H5Pcreate(H5P_FILE_ACCESS);
    
    //STDIO file driver has 2G file size limit on 32 bit Linux system, use SEC2 to solve this problem.
    //H5Pset_fapl_stdio( accprop );
    H5Pset_fapl_sec2( accprop );
    
    hsize_t dims = 1;
	simple_space = H5Screate_simple(1, &dims, NULL);
    
    init();	
}

EMhdf2::~EMhdf2()
{
	H5Sclose(simple_space);
	H5Pclose(accprop);
	if (group >= 0) {
        H5Gclose(group);
    }
    if (file >= 0) {
		H5Fflush(file,H5F_SCOPE_GLOBAL);	// If there were no resource leaks, this wouldn't be necessary...
		H5Fclose(file);
    }
}

int EMhdf2::write_attr(hid_t loc,const char *name,EMObject obj) {
	hid_t type=0;
	hid_t spc=0;
	hsize_t dims=1;
	vector <float> fv;
	switch(obj.get_type()) {
	case EMObject::INT: type=H5Tcopy(H5T_STD_I32LE); spc=H5Scopy(simple_space); break;
	case EMObject::FLOAT: type=H5Tcopy(H5T_IEEE_F32LE); spc=H5Scopy(simple_space); break;
	case EMObject::DOUBLE: type=H5Tcopy(H5T_IEEE_F64LE); spc=H5Scopy(simple_space); break;
	case EMObject::STRING: 
		type=H5Tcopy(H5T_C_S1); 
		H5Tset_size(type,strlen((const char *)obj)+1);
		spc=H5Screate(H5S_SCALAR);
		break;
	case EMObject::FLOATARRAY:
		type=H5Tcopy(H5T_IEEE_F32LE);
		fv=obj;
		dims=fv.size();
		spc=H5Screate_simple(1,&dims,NULL);
		break;
	case EMObject::STRINGARRAY:
	case EMObject::EMDATA:
	case EMObject::XYDATA:
		return -1;
		break;
	case EMObject::UNKNOWN:
		break;
	}

    //we need this delete attribute call here, even we called erase_header()
    //at the biginning of write_header(), since the  "imageid_max" need be updated correctly.
	if( H5Adelete(loc,name) < 0 ) {
#ifdef DEBUGHDF
		LOGERR("Attribute %s deletion error in write_attr().\n", name);
#endif	
	}
	else {
#ifdef DEBUGHDF		
		printf("delete attribute %s successfully in write_attr().\n", name);
#endif	
	} 
	hid_t attr = H5Acreate(loc,name,type,spc,H5P_DEFAULT);

	unsigned int i;
	float f,*fa;
	double d;
	const char *s;
	switch(obj.get_type()) {
	case EMObject::INT:
		i=(int)obj;
		H5Awrite(attr,H5T_NATIVE_INT,&i);
		break;
	case EMObject::FLOAT:
		f=(float)obj;
		H5Awrite(attr,H5T_NATIVE_FLOAT,&f);
		break;
	case EMObject::DOUBLE:
		d=(double)obj;
		H5Awrite(attr,H5T_NATIVE_DOUBLE,&d);
		break;
	case EMObject::STRING: 
		s=(const char *)obj;
//		H5Awrite(attr,H5T_NATIVE_CHAR,s);
		H5Awrite(attr,type,s);
		break;
	case EMObject::FLOATARRAY:
		fa=(float *)malloc(fv.size()*sizeof(float));
		for (i=0; i<fv.size(); i++) fa[i]=fv[i];
		H5Awrite(attr,H5T_NATIVE_FLOAT,fa);
		free(fa);
		break;
//	case EMObject::STRINGARRAY:
//	case EMObject::EMDATA:
//	case EMObject::XYDATA:
//		return -1;
//		break;
	default:
		fprintf(stderr, "Unhandled HDF5 metadata '%s'", name);
		
	}

	herr_t ret1 = H5Tclose(type);
	herr_t ret2 = H5Sclose(spc);
	herr_t ret3 = H5Aclose(attr);
	if(ret1>=0 && ret2>=0 && ret3>=0) {
		return 0;
	}
	else {
		fprintf(stderr, "close error in write_attr()\n");
		return -1;
	}
}

EMObject EMhdf2::read_attr(hid_t attr) {
	hid_t type = H5Aget_type(attr);
	hid_t spc = H5Aget_space(attr);
	H5T_class_t cls = H5Tget_class(type);
	size_t sz = H5Tget_size(type);						// storage size, arrays handled in the 'space'
	hssize_t pts = H5Sget_simple_extent_npoints(spc);	// number of points > 1 if an array of floats
	EMObject ret(0);
	int i;
	float f,*fa;
	double d;
	char *s;
	vector <float> fv(pts);

	switch (cls) {
	case H5T_INTEGER:
		H5Aread(attr,H5T_NATIVE_INT,&i);
		ret=EMObject(i);
		break;
	case H5T_FLOAT:
		if (sz==4) {
			if (pts==1) {
				H5Aread(attr,H5T_NATIVE_FLOAT,&f);
				ret=EMObject(f);
			}
			else {
				fa=(float *)malloc(pts*sizeof(float));
				H5Aread(attr,H5T_NATIVE_FLOAT,fa);
				for (i=0; i<pts; i++) fv[i]=fa[i];
				free(fa);
				ret=EMObject(fv);
			}
		}
		else if (sz==8) {
			H5Aread(attr,H5T_NATIVE_DOUBLE,&d);
			ret=EMObject(d);
		}
		break;
	case H5T_STRING:
		s=(char *)malloc(sz+1);
		H5Aread(attr,type,s);
//		H5Aread(attr,H5T_NATIVE_CHAR,s);
		ret=EMObject(s);
		free(s);
		break;
	default:
		fprintf(stderr, "Unhandled HDF5 metadata %d", cls);
	}

	H5Sclose(spc);
	H5Tclose(type);

	return ret;
}


emdataHeader * EMhdf2::read_header(int image_index, Euler& euler)
{
	emdataHeader * head = new emdataHeader();
	
	int i;
	// Each image is in a group for later expansion. Open the group
	char ipath[50];
	sprintf(ipath,"/MDF/images/%d", image_index);
	hid_t igrp=H5Gopen1(file, ipath);
	
	int nattr=H5Aget_num_attrs(igrp);

	char name[128];
	map<string, EMObject>	attrs;
	for (i=0; i<nattr; i++) {
		hid_t attr=H5Aopen_idx(igrp, i);
		H5Aget_name(attr,127,name);
		if (strncmp(name,"EMAN.", 5)!=0) {
			H5Aclose(attr);
			continue;
		}
		EMObject val=read_attr(attr);
		char attr_name[128] = "";
		strcat(attr_name, name+5);
		string atname = string(attr_name);
		attrs[atname] = val; 
		
		H5Aclose(attr);
	}
	
	head->nx = attrs["nx"];
	head->ny = attrs["ny"];
	head->nz = attrs["nz"];
	head->min = attrs["minimum"];
	head->max = attrs["maximum"];
	head->mean = attrs["mean"];
	head->sigma = attrs["sigma"];
	head->sumsq = attrs["square_sum"];
	head->nimg = attrs["ptcl_repr"];
	
	float pixelx = attrs["spacing_row"];
    float pixely = attrs["spacing_col"];
    float pixelz = attrs["spacing_sec"];
    if (pixelx == pixely && pixelx == pixelz) {
    	head->pixel = pixelx;
    }
       
    head->xorigin = attrs["origin_row"];
    head->yorigin = attrs["origin_col"];
    head->zorigin = attrs["origin_sec"];
    
    strcpy( head->micrograph_id, (const char *)attrs["micrograph_id"] );
    
    head->particle_location_center_x = attrs["particle_location_center_x"];
    head->particle_location_center_y = attrs["particle_location_center_y"];
    
    if(attrs["orientation_convention"].get_type() == EMObject::INT) {    
    	head->orientation_convention = static_cast<Euler::EulerType>( (int)attrs["orientation_convention"]);
    }
    else if(attrs["orientation_convention"].get_type() == EMObject::STRING)  {
    	//this hdf file is created by EMAN2, the orientation convention is a string "EMAN", not enum value.
    }
    
    vector<float> orientation;
    orientation = attrs["orientation"];
    // Teepanis 10/03/06 orientation.size() should be
    // either 3 (Euler) or 4 (Quat and Spin).
    if(orientation.size() == 3 || orientation.size() == 4) {
		switch( head->orientation_convention ) {
	         case Euler::EMAN: case Euler::MRC: case Euler::IMAGIC:
	              euler = Euler(orientation[0], orientation[1], orientation[2], head->orientation_convention );
	              break;
	
	         case Euler::QUAT: case Euler::SPIN:
	              euler = Euler(orientation[0], orientation[1], orientation[2], orientation[3], head->orientation_convention );
	    }
    }
	
	head->center_x = attrs["center_x"];
    head->center_y = attrs["center_y"];
    head->alipeak = attrs["score"];
//	head->nimg = attrs["good"]; //boolean actually
	
	head->ctf[0] = attrs["defocus"];
	head->ctf[1] = attrs["bfactor"];
	head->ctf[2] = attrs["amplitude"];
	head->ctf[3] = attrs["ampcont"];
	head->ctf[4] = attrs["noise1"];
	head->ctf[5] = attrs["noise2"];
	head->ctf[6] = attrs["noise3"];
	head->ctf[7] = attrs["noise4"];
	head->ctf[8] = attrs["voltage"];
	head->ctf[9] = attrs["cs"];
	head->ctf[10] = attrs["apix"];

	head->map_type = static_cast<MapInfoType>((int)attrs["map_type"]); 

	H5Gclose(igrp);

	return head;
}


int EMhdf2::read_data(float* data, int image_index)
{
	char ipath[50];
	sprintf(ipath,"/MDF/images/%d/image", image_index);
	hid_t ds=H5Dopen(file, ipath);
	if (ds<0) {
		fprintf(stderr, "%s Image does not exist\n", filename.c_str());
	}
	hid_t spc=H5Dget_space(ds);
	
	H5Dread(ds,H5T_NATIVE_FLOAT,spc,spc,H5P_DEFAULT,data);
	
	H5Sclose(spc);
	H5Dclose(ds);
	return 0;
}

int EMhdf2::erase_header(int image_index)
{
	if(image_index < 0) return 0; //image_index<0 for appending image, no need for erasing 
	
	init();
	
	int i;
	// Each image is in a group for later expansion. Open the group
	char ipath[50];
	sprintf(ipath, "/MDF/images/%d", image_index);
	hid_t igrp = H5Gopen(file, ipath);
	
	int nattr = H5Aget_num_attrs(igrp);

	char name[128];
	for (i=0; i<nattr; i++) {
		hid_t attr = H5Aopen_idx(igrp, 0); //use 0 as index here, since the H5Adelete() will change the index
		H5Aget_name(attr, 127, name);
		H5Aclose(attr);
		if( H5Adelete(igrp, name) < 0 ) {
			fprintf(stderr, "attribute %s deletion error in erase_header().\n", name);
		}
	}

	H5Gclose(igrp);
	return 0;
}

int EMhdf2::write_header(const emdataHeader & header, int image_index, const Euler& euler)
{
	init();
	
	// If image_index<0 append, and make sure the max value in the file is correct
	// though this is normally handled by EMData.write_image()
	hid_t attr = H5Aopen_name(group,"imageid_max");
	int nimg = read_attr(attr);
	H5Aclose(attr);
	unsigned int i;
	if (image_index<0) image_index=nimg+1;
	if (image_index>nimg) {
		write_attr(group,(const char *)"imageid_max",EMObject(image_index));
	}
	
	// Each image is in a group for later expansion. Open the group, create if necessary
	char ipath[50];
	sprintf(ipath,"/MDF/images/%d", image_index);
	hid_t igrp = H5Gopen(file, ipath);
	
	if (igrp<0) {
		// Need to create a new image group
		igrp=H5Gcreate(file, ipath, 64);		// The image is a group, with attributes on the group
		if (igrp<0) {
			fprintf(stderr,"%s: Unable to add /MDF/images/# to HDF5 file", filename.c_str());
		}
		
		sprintf(ipath,"/MDF/images/%d/image",image_index);
		// Now create the actual image dataset
		hid_t space;
		hid_t ds;
		if (header.nz == 1) { 
			hsize_t dims[2]= { (hsize_t)header.ny, (hsize_t)header.nx };
			space=H5Screate_simple(2,dims,NULL);
		}
		else {
			hsize_t dims[3]= { (hsize_t)header.ny, (hsize_t)header.nx, (hsize_t)header.nx };
			space=H5Screate_simple(3,dims,NULL);
		}
		ds = H5Dcreate(file, ipath, H5T_NATIVE_FLOAT, space, H5P_DEFAULT);
		H5Dclose(ds);
		H5Sclose(space);
	}
	//if group already existed, and the overwiting image is in different size, 
	//need unlink the existing dataset first
	else {
		int nattr=H5Aget_num_attrs(igrp);
		
		char name[128];
		map<string, EMObject>	attrs;
		for (i=0; i<nattr; i++) {
			hid_t attr=H5Aopen_idx(igrp, i);
			H5Aget_name(attr,127,name);
			if (strncmp(name,"EMAN.", 5)!=0) {
				H5Aclose(attr);
				continue;
			}
			EMObject val=read_attr(attr);
			char attr_name[128] = "";
			strcat(attr_name, name+5);
			string atname = string(attr_name);
			attrs[atname] = val; 
			
			H5Aclose(attr);
		}
		
		erase_header(image_index);
		
		if( header.nx != (int)attrs["nx"] 
			|| header.ny != (int)attrs["ny"] 
			|| header.nz != (int)attrs["nz"] ) {
			sprintf(ipath,"/MDF/images/%d/image",image_index);
			H5Gunlink(igrp, ipath);	
				
			// Now create the actual image dataset	
			hid_t space;
			hid_t ds;
			
			if (header.nz == 1) { 
				hsize_t dims[2]= { (hsize_t)header.ny, (hsize_t)header.nx };
				space=H5Screate_simple(2,dims,NULL);
			}
			else {
				hsize_t dims[3]= { (hsize_t)header.ny, (hsize_t)header.nx, (hsize_t)header.nx };
				space=H5Screate_simple(3,dims,NULL);
			}
			ds=H5Dcreate(file,ipath, H5T_NATIVE_FLOAT, space, H5P_DEFAULT );
			H5Dclose(ds);
			H5Sclose(space);
		}
	}
		
	// Write the attributes to the group
	if(is_little_endian()) {
		write_attr(igrp, "EMAN.ImageEndian", "little");
	}
	else {
		write_attr(igrp, "EMAN.ImageEndian", "big");
	}
	
	write_attr(igrp, "EMAN.nx", header.nx); 
	write_attr(igrp, "EMAN.ny", header.ny);
	write_attr(igrp, "EMAN.nz", header.nz);
	write_attr(igrp, "EMAN.minimum", header.min);
	write_attr(igrp, "EMAN.maximum", header.max);
	write_attr(igrp, "EMAN.mean", header.mean);
	write_attr(igrp, "EMAN.sigma", header.sigma);
	write_attr(igrp, "EMAN.square_sum", header.sumsq);
	write_attr(igrp, "EMAN.pixel", header.pixel);
	write_attr(igrp, "EMAN.origin_row", header.xorigin);
	write_attr(igrp, "EMAN.origin_col", header.yorigin);
	write_attr(igrp, "EMAN.origin_sec", header.zorigin);
	write_attr(igrp, "EMAN.micrograph_id", header.micrograph_id);
	write_attr(igrp, "EMAN.particle_location_center_x", header.particle_location_center_x);
	write_attr(igrp, "EMAN.particle_location_center_y", header.particle_location_center_y);
	write_attr(igrp, "EMAN.center_x", header.center_x); 
	write_attr(igrp, "EMAN.center_y", header.center_y); 
	write_attr(igrp, "EMAN.score", header.alipeak); 
	write_attr(igrp, "EMAN.good", header.nimg); 
	write_attr(igrp, "EMAN.map_type", header.map_type);
	
	write_attr(igrp, "EMAN.defocus", header.ctf[0]);
	write_attr(igrp, "EMAN.bfactor", header.ctf[1]);
	write_attr(igrp, "EMAN.amplitude", header.ctf[2]);
	write_attr(igrp, "EMAN.ampcont", header.ctf[3]);
	write_attr(igrp, "EMAN.noise1", header.ctf[4]);
	write_attr(igrp, "EMAN.noise2", header.ctf[5]);
	write_attr(igrp, "EMAN.noise3", header.ctf[6]);
	write_attr(igrp, "EMAN.noise4", header.ctf[7]);
	write_attr(igrp, "EMAN.voltage", header.ctf[8]);
	write_attr(igrp, "EMAN.cs", header.ctf[9]);
	write_attr(igrp, "EMAN.apix", header.ctf[10]);
	
	write_attr(igrp, "EMAN.orientation_convention", header.orientation_convention);
	
	vector<float> orientation;
	orientation.push_back(euler.alt());
	orientation.push_back(euler.az());
	orientation.push_back(euler.phi());
	orientation.push_back(0);
	write_attr(igrp, "EMAN.orientation", orientation);
		
	H5Gclose(igrp);
	
	return 0;
}
    
int EMhdf2::write_data(const float* const data, int image_index)
{	
	if( image_index < 0 ) {
		hid_t attr = H5Aopen_name(group, "imageid_max");
		image_index = read_attr(attr);
		H5Aclose(attr);
	}
	
	char ipath[50];
	sprintf(ipath, "/MDF/images/%d/image", image_index);
	
	hid_t ds = H5Dopen(file, ipath);
	if (ds<0) {
		fprintf(stderr, "%s: Image dataset does not exist", filename.c_str());
	}
	
	hid_t spc = H5Dget_space(ds);
	H5Dwrite(ds, H5T_NATIVE_FLOAT, spc, spc, H5P_DEFAULT, data);
	
	H5Sclose(spc);
	H5Dclose(ds);
	return 0;
}

int EMhdf2::get_num_dataset()
{
	init();
	hid_t attr = H5Aopen_name(group,"imageid_max");
	int n = read_attr(attr);
	H5Aclose(attr);

	return n+1;
}

