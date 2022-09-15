#include "emdata_.h"

const int NVALIDFMT = 7;// was 10;
const char *VALIDFMT[] = { "gif","pgm","ppm","pnm","rgb","tga","bmp" };

Euler::~Euler() {}

Euler::Euler() { 
	alt_=0.; az_=0.; phi_=0.;
	init();
} // already rectified

	
Euler::Euler(float alt, float az, float phi){// EMAN by default
	setAngle(alt,az,phi);
	init();
}

Euler::Euler(float a1, float a2, float a3, int typ) { 
	setAngle(a1,a2,a3,typ);
	init();
} // end EMAN, IMAGIC, SPIDER constructor

Euler::Euler(float e1, float e2, float e3, float e0, int typ) {
	setAngle(e1,e2,e3,e0,typ); 
	init();
}           // end quaternion, and spin constructor	

Euler::Euler(float *m){
	setAngle(m);
	init();
}

void Euler::init() {
	strcpy(sym,"i");
	nsym=1;
	csym=0;
}

void Euler::setAngle(float a1, float a2, float a3, int typ){ 
	alt_ =a1 ; az_ = a2 ; phi_ = a3;
	if (typ == SPIDER) { az_ = a2-PI/2.-PI; phi_ = a3+PI/2.-PI;}    // then a1=theta, a2=phi, a3=gamma
    	if (typ !=EMAN && typ !=IMAGIC && typ != SPIDER) {	 
			printf("unknown declaration of type for Euler constructor /n");} 
	rectify();
}  
	
void  Euler::setAngle(float e1, float e2, float e3, float e0, int typ) {
	typ_ = typ ;
	if (typ !=SPIN && typ!=QUAT && typ != SGIROT) {	 
		printf("unknown declaration of type for Euler constructor /n");}
//   if not already, put in S4 language
	if (typ == SPIN || typ==SGIROT) { 
		float normal = sqrt(e1*e1+e2*e2+e3*e3);
		float sinQo2 = sin(e0/2.);		e0 = cos(e0/2.);
		e1 *= sinQo2/normal; e2 *= sinQo2/normal; e3 *= sinQo2/normal; }
	float w = e0, x = e1, y = e2, z = e3;
	m_[0] = 1-2*y*y-2*z*z;;
	m_[1] = 2*x*y + 2*w*z;
	m_[2] = 2*x*z - 2*w*y;
	m_[3] = 2*x*y - 2*w*z;
	m_[4] = 1-2*x*x-2*z*z;
	m_[5] = 2*y*z + 2*w*x;
	m_[6] = 2*x*z + 2*w*y;
	m_[7] = 2*y*z - 2*w*x;
	m_[8] = 1 - 2*x*x - 2*y*y;
	
	setAngle(m_);
	
	rectify();
	//if (typ==SGIROT) inverse();
}           // end quaternion, and spin setAngle

void Euler::setAngle(float *m){
	typ_ = MATRIX;
	alt_ = acos( m[8]    );
	az_  = atan2(m[6],-m[7]);
	phi_ = atan2(m[2], m[5]);
	if (m[8] > 0.999999) {
		alt_=0;
		az_ = atan2(m[1],m[4]);
		phi_ =0;}
	if (m[8] < -0.999999) {
		alt_ = PI;
		az_ = atan2(m[1],-m[4]);
		phi_ =0;
	}		 		 
	rectify();
}

void Euler::rectify(){
	if (!goodf(&phi_)) phi_=0;
	if (!goodf(&alt_)) alt_=0;
	if (!goodf(&az_)) az_=0;
	if (sqrt(alt_*alt_) <= .0001) { az_ += phi_; phi_=0;}
	
	phi_=fmod(phi_,(float)(2.0*PI));
	if (phi_<0) phi_+=2.0*PI;
	
	az_=fmod(az_,(float)(2.0*PI));
	if (az_<0) az_+=2.0*PI;
}


bool is_machine_big_endian()
{
	int one = 1;
	char* p_one = (char*) (&one);
	
	if(p_one[0] == 1 && p_one[1] == 0 && p_one[2] == 0 && p_one[3] == 0) {
	    return  false;
	}
	else {
	    return true;
	}
}

emdata::emdata() {
	//int i;
	
	header.nx=header.ny=header.nz=0;
	header.dx=header.dy=header.dz=0;
	header.min=header.max=header.mean=header.sigma=0;
	header.pixel=1.0f;
	header.xorigin=header.yorigin=header.zorigin=0.0;
	header.path[0]=0;
	header.pathnum=0;
	header.name[0]=0;
	header.rocount=0;
	header.micrograph_id[0]= '\0';
	header.particle_location_center_x = 0.0;
	header.particle_location_center_y = 0.0;
	header.center_x = -999.0;	// negative means no initialization
	header.center_y = -999.0;
	header.map_type = NORMAL;
	header.orientation_convention = Euler::EMAN;

	rdata=NULL;
	euler.setAngle(0,0,0);
	mrch=NULL;
	imagich=NULL;
	spiderH=NULL;
	fft=NULL;	
	rfp=NULL;
	
}

emdata::~emdata() {
	//if (rdata) Sfree(rdata,1);
	if (rdata) free(rdata);
	//if (mrch) Sfree(mrch,1);
	if (mrch) free(mrch);
	if (spiderH) free(spiderH);
	//if (imagich) Sfree(imagich,1);
	if (imagich) free(imagich);
	if (fft) delete fft;
	if (rfp) delete rfp;
}


int emdata::setSize(int x,int y, int z) {
	if (x<=0) x=1;
	if (y<=0) y=1;
	if (z<=0) z=1;
	if (y==1&&z!=1) { y=z; z=1; }
	if (x==1&&y!=1) { x=y; y=1; }
	header.nx=x;
	header.ny=y;
	header.nz=z;
	rdata=(float *)realloc(rdata,x*y*z*sizeof(float));
	return 0;
}

int emdata::readImage(const char *filespec, int n, int nodata)
{	
	if (!filespec) { printf("NULL readimage()\n"); return(-1); }

    if (EMhdf::is_valid(filespec)) {	//valid HDF file
		int ret = EMhdf::init_test(filespec);
        int err = 0;
        if( ret == 0 ) {	//old-style HDF format
			EMhdf hdf_file(filespec);
            err = readHDF(filespec, n, nodata);
            if (err) {
                printf("failed in readHDF()\n");
                setSize(10, 10, 1);
                zero();
            }
        }
        else {	//new-style HDF format
			err = readHDF2(filespec, n, nodata);
            if (err) {
				printf("failed in readHDF2()\n");
                setSize(10, 10, 1);
                zero();
            }
        }
        return err;
    }

	int r=readMRC(filespec,nodata,n);
	if (r) { setSize(10,10,1); zero(); }
	return r;

}

void emdata::zero(float sig) {
	if (sig==0) memset(rdata,0,header.nx*header.ny*header.nz*sizeof(float));
	else {
		int i;
		for (i=0; i<header.nx*header.ny*header.nz; i++) rdata[i]=grand(0,sig);
	}
}


int emdata::readHDF(const char *fsp, int image_index, int nodata)
{
    if (image_index < 0)image_index = 0;

    if (strrchr(fsp,'/') == NULL) setName(fsp);
    else setName(strrchr(fsp,'/') + 1);
    
    setPath(fsp);
    
    EMhdf hdf_file(fsp, "r");
    
    int ndim = 0;
    int* dims = hdf_file.read_dims(image_index, &ndim);

    if (ndim == 3) {
		setSize(dims[0], dims[1], dims[2]);
    }
    else if (ndim == 2) {
		setSize(dims[0], dims[1], 1);
    }
    else {
		setSize(10, 10, 1);
		zero();
		return -1;
    }


    if (!nodata) { 
		int err = hdf_file.read_data(image_index, &rdata);
	 	if (err) {
	   		setSize(10, 10, 1);
	   		zero();
	   		return -1;
	 	}
    }    

    float pixelx = hdf_file.read_float_attr(image_index, "spacing_row");
    float pixely = hdf_file.read_float_attr(image_index, "spacing_col");
    float pixelz = hdf_file.read_float_attr(image_index, "spacing_sec");
    if (pixelx == pixely && pixelx == pixelz)
      header.pixel = pixelx;

    else {
        setSize(10, 10, 1);
		zero();
		return -1;
    }

    header.xorigin = hdf_file.read_float_attr(image_index, "origin_row");
    header.yorigin = hdf_file.read_float_attr(image_index, "origin_col");
    header.zorigin = hdf_file.read_float_attr(image_index, "origin_sec");
    header.min = hdf_file.read_float_attr(image_index, "minimum");
    header.max = hdf_file.read_float_attr(image_index, "maximum");
    header.mean = hdf_file.read_float_attr(image_index, "average");
    header.sigma = hdf_file.read_float_attr(image_index, "std"); 
    header.nimg = hdf_file.read_int_attr(image_index, "ptcl_repr");

    const char *tmp = hdf_file.read_string_attr(image_index, "micrograph_id");
    setMicrographId(tmp); 
    delete[] tmp;
    tmp = 0;

    header.particle_location_center_x = hdf_file.read_float_attr(image_index, "particle_location_center_x");
    header.particle_location_center_y = hdf_file.read_float_attr(image_index, "particle_location_center_y");

    float orientation[4];
    header.orientation_convention = static_cast<Euler::EulerType>( hdf_file.read_euler_attr(image_index, "orientation_convention"));
    hdf_file.read_array_attr(image_index, "orientation", orientation);
    switch(header.orientation_convention) {
         case Euler::EMAN: case Euler::MRC: case Euler::IMAGIC:
              euler = Euler(orientation[0], orientation[1], orientation[2], header.orientation_convention );
              break;

         case Euler::QUAT: case Euler::SPIN:
              euler = Euler(orientation[0], orientation[1], orientation[2], orientation[3], header.orientation_convention );
    }

    header.center_x = hdf_file.read_float_attr(image_index, "center_x");
    header.center_y = hdf_file.read_float_attr(image_index, "center_y");
 
    header.alipeak = hdf_file.read_float_attr(image_index, "score");
    //nimg<=0 bad, >0 good.
    header.nimg = hdf_file.read_int_attr(image_index, "good"); //boolean actually
 
    CTFParam* ctfit = new CTFParam();
    if (hdf_file.read_ctfit(image_index, ctfit) != 1) {
      header.ctf[0] = ctfit->defocus;
      header.ctf[1] = ctfit->bfactor;
      header.ctf[2] = ctfit->amplitude;
      header.ctf[3] = ctfit->ampcont;
      header.ctf[4] = ctfit->noise1;
      header.ctf[5] = ctfit->noise2;
      header.ctf[6] = ctfit->noise3;
      header.ctf[7] = ctfit->noise4;
      header.ctf[8] = ctfit->voltage;
      header.ctf[9] = ctfit->cs;
      header.ctf[10] = ctfit->apix;
      // this is redundant but needed for legacy
	  char buf[160];
	  sprintf(buf,"!-%1.3g %1.3g %1.3g %1.3g %1.3g %1.3g %1.3g %1.3g %1.3g %1.3g %1.3g",
			header.ctf[0],header.ctf[1],header.ctf[2],header.ctf[3],header.ctf[4],header.ctf[5],header.ctf[6],header.ctf[7],header.ctf[8],header.ctf[9],header.ctf[10]);
	  for(int i=0;i<80;i++) header.name[i] = buf[i];
	} 
    delete ctfit;
    ctfit = 0;

    header.map_type = static_cast<MapInfoType>(hdf_file.read_mapinfo_attr(image_index, "map_type"));  
    
    header.pathnum = 0;
    header.dx = header.dy = header.dz = 0;
    return 0;
}


int emdata::readHDF2(const char *fsp, int image_index, int nodata)
{
    if (image_index < 0) {
		image_index = 0;
    }
    
    if (strrchr(fsp,'/') == NULL) {
		setName(fsp);
    }
    else {
		setName(strrchr(fsp,'/') + 1);
    }
    setPath(fsp);
    
    EMhdf2 * hdf_file = new EMhdf2(fsp, "r");
    emdataHeader * head = hdf_file->read_header(image_index, this->euler);
    
    if( head->nz > 1 ) {
    	setSize(head->nx, head->ny, head->nz);
    }
    else if( head->ny > 1 ) {
    	setSize(head->nx, head->ny, 1);
    }
    else if( head->nx > 1 ) {
    	setSize(head->nx, 1, 1);
    }
    else {
    	setSize(10, 10, 1);
    	zero();
    	return -1;
    }
    
    if( !nodata )  {
    	int err = hdf_file->read_data(rdata, image_index);
    	if( err ) {
    		setSize(10, 10, 1);
    		zero();
    		return -1;
    	}
    }
    header = *head;

    if(head) {
    	delete head;
    	head = 0;
    }
    if(hdf_file) {
    	delete hdf_file;
    	hdf_file = 0;
    }
	return 0; 
}

//Can not work when read complex
int emdata::readMRC(const char *fsp,int nodata, int n) {
	FILE *in;
	int i,j,l,pipe=0;
	int ord=0;	// 1 if the file order is MSB first
	int mord=0; // 1 if the machine order is MSB first
	float f,a,p;
	char s[800],*m2,w;
	unsigned char *cdata,*cdata2,c,t;
	short *sdata = 0;;
	unsigned short *usdata = 0;
	static const char *str[] = { "char", "short", "float", "short complex", "float complex" };
	static const char *str2[] = { "" , "byte reversed" , "" };
	static const char *str3[] = { "" , "compressed" };
	int nx = header.nx;
	int ny = header.ny;
	int nz = header.nz;
	
	if (!mrch) mrch=(mrcH *)malloc(sizeof(mrcH));

	if (fsp==NULL) { in=stdin; }
	else {
		if (strcmp(&fsp[strlen(fsp)-3],".gz")==0 || strcmp(&fsp[strlen(fsp)-2],".Z")==0) {
			sprintf(s,"zcat %s",fsp);
			in=popen(s,"rb");
			pipe=1;
		}
		if (strcmp(&fsp[strlen(fsp)-3],".bz")==0 || strcmp(&fsp[strlen(fsp)-4],".bz2")==0) {
			sprintf(s,"bzcat %s",fsp);
			in=popen(s,"rb");
			pipe=1;
		}
		else in=fopen(fsp,"rb");
		if (in==NULL) {
			printf("Cannot open MRC file\n");
			zero();
			return -1;
		}
		//if (strrchr(fsp,'/')==NULL) setName(fsp);
		//else setName(strrchr(fsp,'/')+1);
		setPath(fsp);
	}

	if (!fread(mrch,sizeof(mrcH),1,in)) 
	{ 
		setSize(10,10,1); 
		zero(); 
		return -1; 
	}
	
	mrch->labels[0][79]=0;
	setName(mrch->labels[0]);
	if (header.name[0]=='!' && header.name[1]=='-') {
		i=sscanf(header.name+2," %f %f %f %f %f %f %f %f %f %f %f",
			header.ctf,header.ctf+1,header.ctf+2,header.ctf+3,header.ctf+4,header.ctf+5,header.ctf+6,header.ctf+7,header.ctf+8,header.ctf+9,header.ctf+10); 
		if (i!=11) printf("Incomplete CTF info: %s\n",header.name);
	}
	if (header.name[0]=='!' && header.name[1]=='$') {
		sscanf(header.name+2," %f,%f",header.ctf,header.ctf+10);
	}

	// This swaps the byte order in the header on appropriate machines
	if (is_machine_big_endian()) mord=1;
	else mord=0;

	m2=(char *)mrch;
	if (m2[0]==0 && m2[1]==0) ord=1;
	if (mord ^ ord) {
		for (i=0; i<56*4; i+=4) {
			w=m2[i]; m2[i]=m2[i+3]; m2[i+3]=w;
			w=m2[i+1]; m2[i+1]=m2[i+2]; m2[i+2]=w;
		}
	}
	
	if (mrch->ny>0xffff || mrch->nz>0xffff || mrch->mode>4 || mrch->nlabl>11) {
		setSize(10,10,1); printf("Invalid MRC file\n"); zero(); return -1;
	}

	if (mrch->ny==1) { printf("This appears to be a 1d file. I only read 2d files.\n"); zero(); return -4; }
	//if (mrch->nz!=1) { printf("This appears to be a 3d file. I only read 2d files.\n"); return 0; }
	
	int startz;
	if (mrch->mode==MODE_float_COMPLEX || mrch->mode==MODE_short_COMPLEX) mrch->nx*=2;
	nx=mrch->nx;
	ny=mrch->ny;
	if (n>=0 && n< mrch->nz) {
		nz=1;
		startz=n;
	}
	else {
		//if(n!=-1) printf("WARNING from readMRC: section n=%d is out of valid range [0,%d] and thus ignored, whole map will be read in\n",n,mrch->nz-1); 	
		nz=mrch->nz;
		startz=0;
	}
	
	if (mrch->nz==1) sprintf(s,"Read a %s %s %dx%d 2D %s MRC File\n",
		str2[ord+mord],str3[pipe],nx,ny,str[mrch->mode]);
	else sprintf(s,"Read a %dx%dx%d 3D MRC File\n",nx,ny,nz);
	if (mrch->xlen==0) mrch->xlen=1.0;
	if (mrch->ylen==0) mrch->ylen=1.0;
	if (mrch->zlen==0) mrch->zlen=1.0;

	setSize(nx,ny,nz);

	fseek(in,sizeof(mrcH)+mrch->nsymbt,SEEK_SET);

	cdata=(unsigned char *)rdata;
	sdata=(short *)rdata;
	usdata=(unsigned short *)rdata;
	off_t cur_offset;

	switch(mrch->mode) 
	{
		case MODE_char:
			fseek(in,((off_t)startz)*nx*ny,SEEK_CUR);
			for (j=0; j<ny; j++) {
				if ((i=fread(&cdata[j*nx*nz],nx,nz,in))!=nz) {
					sprintf(s,"Incomplete data read %d/%d blocks\n",i,ny*nz);
				}
			}
			for (i=nx*ny*nz-1; i>=0; i--) rdata[i]=(float)cdata[i]/100.0-1.28;
			break;
		case MODE_short:
			fseek(in,((off_t)startz)*nx*ny*sizeof(short),SEEK_CUR);
			for (j=0; j<ny; j++) {
				if ((i=fread(&sdata[j*nx*nz],nx*sizeof(short),nz,in))!=nz) {
					sprintf(s,"Incomplete data read %d/%d blocks\n",i,ny*nz);
				}				
			}
			if ((ord+mord)&1) {
				for (i=nx*ny*nz-1; i>=0; i--) {
					t=cdata[i*2];
					cdata[i*2]=cdata[i*2+1];
					cdata[i*2+1]=t;
					rdata[i]=(float)(sdata[i]);
				} 
			}
			else {
				for (i=nx*ny*nz-1; i>=0; i--)  {
					rdata[i]=(float)(sdata[i]); 
				}
			}
			break;
		case MODE_ushort:
			fseek(in,((off_t)startz)*nx*ny*sizeof(unsigned short),SEEK_CUR);		
			for (j=0; j<ny; j++) {
				if ((i=fread(&cdata[j*nx*nz],nx*sizeof(unsigned short),nz,in))!=nz) {
					sprintf(s,"Incomplete data read %d/%d blocks\n",i,ny*nz);
				}				
			}
			if ((ord+mord)&1) {
				for (i=nx*ny*nz-1; i>=0; i--) {
					t=cdata[i*2];
					cdata[i*2]=cdata[i*2+1];
					cdata[i*2+1]=t;
					rdata[i]=(float)(usdata[i]);
				} 
			}
			else {
				for (i=nx*ny*nz-1; i>=0; i--)  {
					rdata[i]=(float)(usdata[i]); 
				}
			}
			
			break;
		case MODE_float:
			cur_offset = ((off_t)startz)*nx*ny*sizeof(float);
			fseek(in,cur_offset,SEEK_CUR);
			cdata=(unsigned char *)rdata;
			for (j=0; j<ny; j++) {
				//printf("fuck!\n");
				if ((i=fread(&rdata[j*nx*nz],nx*sizeof(float),nz,in))!=nz) {
					sprintf(s,"Incomplete data read %d/%d blocks\n",i,ny*nz);
				}
			}
			
			// swap order if machine and file have different order
			if (ord+mord==1) {
				for (i=nx*ny*nz*4-4; i>=0; i-=4) {
					c=cdata[i]; cdata[i]=cdata[i+3]; cdata[i+3]=c;
					c=cdata[i+1]; cdata[i+1]=cdata[i+2]; cdata[i+2]=c;
				}
			}
			break;
		case MODE_short_COMPLEX:
			fseek(in,((off_t)startz)*nx*ny*sizeof(short),SEEK_CUR);
			cdata=(unsigned char *)rdata;
			for (j=0; j<ny; j++) {
				if ((i=fread(&cdata[j*nx*nz*2],nx*sizeof(short),nz,in))!=nz) {
					sprintf(s,"Incomplete data read %d/%d blocks\n",i,ny*nz);
				}
			}
			if (ord) {
				for (i=nx*ny*nz-1; i>=0; i--) {
					rdata[i]=(float)((cdata[i*2]<<8)+cdata[i*2+1]); 
				}
			}
			else {
				for (i=nx*ny*nz-1; i>=0; i--) {
					rdata[i]=(float)((cdata[i*2+1]<<8)+cdata[i*2]); 
				}
			}
			break;
		case MODE_float_COMPLEX:
			fseek(in,((off_t)startz)*nx*ny*sizeof(float),SEEK_CUR);
			cdata=(unsigned char *)rdata;
			clearerr(in);
			for (j=0; j<ny; j++) {
				if ((i=fread(&rdata[j*nx*nz],nx*sizeof(float),nz,in))!=nz) {
					sprintf(s,"Incomplete data read %d/%d blocks\n",i,ny*nz);
				}
			}
			if (ord+mord==1) {
				for (i=nx*ny*nz*4-4; i>=0; i-=4) {
					c=cdata[i]; cdata[i]=cdata[i+3]; cdata[i+3]=c;
					c=cdata[i+1]; cdata[i+1]=cdata[i+2]; cdata[i+2]=c;
				}
			}
			if (ferror(in)) { zero(); return -1; }
			break;
		default: 
			printf("Bad MRC file\n"); 
			zero(); 
			return -1;
	}
	
  	if(fsp!=NULL)  if (pipe) pclose(in); else fclose(in);
	header.pixel=mrch->xlen/(float)(nx);	// somebody changed the writeMRC to use nx instead of nx-1, so this is update to use nx too for self consistence
	header.xorigin=mrch->xorigin;
	header.yorigin=mrch->yorigin;
	header.zorigin=mrch->zorigin;
	header.dx=header.dy=header.dz=0;
	setRAlign(0,0,0);
	return 0;
}

void emdata::ap2ri() {		// the opposite
	int i;
	float f,*data;
	
	for (i=0; i<header.nx*header.ny*header.nz; i+=2) {
		f=rdata[i]*sin(rdata[i+1]);
		rdata[i]=rdata[i]*cos(rdata[i+1]);
		rdata[i+1]=f;
	}
	
}


void emdata::rotate(float e)
{
	float *ddata=(float *)malloc(header.nx*header.ny*header.nz*sizeof(float));

	float cose=cos(e);
	float sine=sin(e);
	float x2c=header.nx/2.0-header.dx;
	float y2c=header.ny/2.0-header.dy;
	int i,j,k0,k1,k2,k3;
	float x,y,ii,jj,t,u;
	float p0,p1,p2,p3;
	for (j=0,y=-header.ny/2.0; j<header.ny; j++,y+=1.0) {
		for (i=0,x=-header.nx/2.0; i<header.nx; i++,x+=1.0) {

		float x2 = (cose*x+sine*y)+x2c;
		float y2 = (-sine*x+cose*y)+y2c;

		if (x2<0||x2>header.nx-1.0||y2<0||y2>header.ny-1.0) {
			ddata[i+j*header.nx]=0;
			continue;
		}

		ii=floor(x2);
		jj=floor(y2);
		k0=ii+jj*header.nx;
		k1=k0+1;
		k2=k0+header.nx+1;
		k3=k0+header.nx;

		if (ii==header.nx-1) { k1--; k2--; }
		if (jj==header.ny-1) { k2-=header.nx; k3-=header.nx; }

		// no interpolation
		//ddata[i+j*nx]=sdata[ii+jj*nx];

		// bilinear interpolation of raw data
		t=(x2-(float)ii);
		u=(y2-(float)jj);
		float tt=1.0-t;
		float uu=1.0-u;

		p0=rdata[k0]*tt*uu;
		p1=rdata[k1]*t*uu;
		p3=rdata[k3]*tt*u;
		p2=rdata[k2]*t*u;
		ddata[i+j*header.nx]=p0+p1+p2+p3;
		}
	}
	setXYZOrigin(	header.xorigin-header.dx*header.pixel,
				header.yorigin-header.dy*header.pixel,
				header.zorigin-header.dz*header.pixel
			);

	free(rdata);
	rdata=ddata;
	header.dx=header.dy=header.dz=0;
	setRAlign(0,0,0);
	
}


EMObject::operator int () const
{
	if (type == INT) {
		return n;
	}
	else if (type == FLOAT) {
		return (int) f;
	}
	else if (type == DOUBLE) {
		return (int) d;
	}
	else {
		if (type != UNKNOWN) {
			fprintf(stderr, "Cannot convert to int this data type\n");
		}
	}
	return 0;
}

EMObject::operator float () const
{
	if (type == FLOAT) {
		return f;
	}
	else if (type == INT) {
		return (float) n;
	}
	else if (type == DOUBLE) {
		return (float) d;
	}
	else {
		if (type != UNKNOWN) {
			fprintf(stderr, "Cannot convert to float from this data type\n");
		}
	}

	return 0;
}

EMObject::operator double () const
{
	if (type == DOUBLE) {
		return d;
	}
	else if (type == INT) {
		return (double) n;
	}
	else if (type == FLOAT) {
		return (double) f;
	}
	else {
		if (type != UNKNOWN) {
			fprintf(stderr, "Cannot convert to double from this data type\n");
		}
	}
	return 0;
}

EMObject::operator  const char *() const
{
	if (type != STRING) {
		if (type != UNKNOWN) {
			fprintf(stderr, "Cannot convert to string from this data type\n");
		}
		return "";
	}
	return str.c_str();
}


EMObject::operator vector < float > () const
{
	if (type != FLOATARRAY) {
		if (type != UNKNOWN) {
			fprintf(stderr, "Cannot convert to vector<float> from this data type\n");
		}
		return vector < float >();
	}
	return farray;
}

EMObject::operator vector<string> () const
{
	if (type != STRINGARRAY) {
		if (type != UNKNOWN) {
			fprintf(stderr, "Cannot convert to vector<string> from this data type\n");
		}
		return vector<string>();
	}
	return strarray;
}

    