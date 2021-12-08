#ifndef DATATYPE_H
#define DATATYPE_H

/* for MRC header */
#define MRC_LABEL_SIZE      80   
#define MRC_USER            25
#define MRC_NUM_LABELS      10

/*  The different modes supported by the MRC format. Values used by IMOD */
#define MODE_char           0
#define MODE_short          1
#define MODE_float          2
#define MODE_short_COMPLEX  3  
#define MODE_float_COMPLEX  4
#define MODE_ushort			6		//non-standard
#define MODE_uchar3			16		//unsigned char * 3, for rgb data, non-standard
/* All modes before MODE_short_COMPLEX must be real, all those after it must be complex. */


/*******************************************************************************
Protected Data Structures
*******************************************************************************/
struct mrcH
{
    int         nx;                 /* Number of columns -      */
                                    /* Fastest changing in map. */
    int         ny;

    int         nz;                 /* Number of sections -     */
                                    /* slowest changing in map. */

    int         mode;               /* See modes above. */

    int         nxstart;            /* No. of first column in map, default 0.*/
    int         nystart;            /* No. of first row in map, default 0.*/
    int         nzstart;            /* No. of first section in map,default 0.*/

    int         mx;                 /* Number of intervals along X. */
    int         my;                 /* Number of intervals along Y. */
    int         mz;                 /* Number of intervals along Z. */

    float       xlen;               /* Cell dimensions (Angstroms). */
    float       ylen;               /* Cell dimensions (Angstroms). */
    float       zlen;               /* Cell dimensions (Angstroms). */

    float       alpha;              /* Cell angles (Degrees). */
    float       beta;               /* Cell angles (Degrees). */
    float       gamma;              /* Cell angles (Degrees). */

    int         mapc;               /* Which axis corresponds to Columns.  */
                                    /* (1,2,3 for X,Y,Z.                   */
    int         mapr;               /* Which axis corresponds to Rows.     */
                                    /* (1,2,3 for X,Y,Z.                   */
    int         maps;               /* Which axis corresponds to Sections. */
                                    /* (1,2,3 for X,Y,Z.                   */

    float       amin;               /* Minimum density value. */
    float       amax;               /* Maximum density value. */
    float       amean;              /* Mean density value.    */

    int         ispg;               /* Space group number (0 for images). */

    int         nsymbt;             /* Number of chars used for storing   */
                                    /* symmetry operators. */


    int         user[MRC_USER];  
	
    float       xorigin;            /* X origin. */
    float       yorigin;            /* Y origin. */
    float       zorigin;            /* Y origin. */

    char        map[4];				/* constant string "MAP "  */
    int         machinestamp;	    /* machine stamp in CCP4 convention: big endian=0x11110000 little endian=0x4444000 */
                                    /* The specification is ambiguous, so we are writing 0x11111111 and 0x44444444 */

    float       rms;                /* 	rms deviation of map from mean density */

    int         nlabl;              /* Number of labels being used. */
                                    /* 10 text labels of 80 characters each. */
    char        labels[MRC_NUM_LABELS][MRC_LABEL_SIZE];
    
};


struct imagicH {
    int imnum;	// image number
    int count;	// total # images in file (1st record only)
    int error;
    int headrec;	// # header records/image
    int mday;
    int month;
    int year;
    int hour;
    int minute;
    int sec;
    int reals;	// image size in reals
    int pixels;	// image size in pixels
    int ny;
    int nx;
    char type[4];	// usu 'REAL', or 'INTG'
    int ixold;
    int iyold;
    float avdens;	// average density
    float sigma;	// deviation of density
    float varia;	// variance
    float oldav;	// old average dens
    float max;	// max dens
    float min;	// min dens
    int complex;	// ?
    float cellx;
    float celly;
    float cellz;
    float cella1;
    float cella2;
    char label[80];	//30
    int SPACE[8];	//50
    float MRC1[4];	//58
    int MRC2;	//62
    int SPACE2[7];
    int lbuf;		// effective buffer len = nx
    int inn;		// lines in buffer = 1
    int iblp;		// buffer lines/image = ny
    int ifb;		// 1st line in buf = 0
    int lbr;		// last buf line read = -1
    int lbw;		// last buf line written = 0
    int lastlr;		// last line called for read = -1
    int lastlw;		// last line called for write = 1
    int ncflag;		// decode to complex = 0
    int num;		// file number = 40 (?)
    int nhalf;		// leff/2
    int ibsd;		// record size for r/w (words) = nx*2
    int ihfl;		// file # = 8
    int lcbr;		// lin count read buf = -1
    int lcbw;		// lin count wr buf = 1
    int imstr;		// calc stat on rd = -1
    int imstw;		// calc stat on wr = -1
    int istart;		// begin line in buf = 1
    int iend;		// end line in buf = nx
    int leff;		// eff line len	= nx
    int linbuf;		// line len (16 bit) nx *2
    int ntotbuf;	// total buf in pgm = -1
    int SPACE3[5];
    int icstart;	// complex line start = 1
    int icend;		// complex line end = nx/2
    int rdonly;		// read only = 0

    int clsrep;		// EMAN specific, classes represented with 0x7a6b5c00 mask
    int emanmisc[6];
    float qual[50];	// quality info from EMAN classification
    int cls[50];	// number of best class
    int flags[50];	// eman flags
};

struct SPIDERH {
    float nslice;		//0    1 for images
    float ny;
    float u1;			//unused, actually used for SPIDER --Grant
    float u2;			
    float type;		//4   type of data 1=2d image,3=volume,-11=2d fft odd, -12=2d fft even, -21=3d fft odd, -22=3dfft even 
    float mmvalid;		// 1 if max and min are valid
    float max;
    float min;
    float mean;		//8
    float sigma;		// std dev, -1=unknown
    float u3;
    float nx;
    float headrec;		//12   # records in header
    float angvalid;		// 1 if tilt angles have been computed
    float phi,theta;
    float gamma;		//16
    float dx,dy,dz;		// translations
    float scale;		//20
    float headlen;		// in bytes
    float reclen;
    float istack;		// file contains a stack of images (2 in overall header, 0 in images)
    float inuse;		//24    indicates that the current image is present in the stack
    float maxim,imgnum,lastindx, u6,u7;
    float ang2valid;	// additional angles
    float phi1,theta1,psi1;
    float phi2,theta2,psi2;
    char u8[48];
    float xf[27];		// Jose Maria's transforms
    float u9[135];
    char date[11];
    char time[8];
    char name[160];
};
enum MapInfoType {
    NORMAL,
    ICOS2F_FIRST_OCTANT,
    ICOS2F_FULL,
    ICOS2F_HALF,
    ICOS3F_HALF,
    ICOS3F_FULL,
    ICOS5F_HALF,
    ICOS5F_FULL
};

struct gatanH {
    //short tag;
    //short tagl;
    short version;
    short un1;
    short un2;
    short nx;
    short ny;
    short len;	// 1=byte, 2=short, 4=float,long
    short type;	// 1=short, 2=float, 3=complex, 5=packed complex, 6=byte, 7=long
};


struct CTFParam {
    float defocus;
    float bfactor;
    float amplitude;
    float ampcont;
    float noise1;
    float noise2;
    float noise3;
    float noise4;
    float voltage;
    float cs;
    float apix;
};

enum ImageType {
    IMAGE_MRC = 0,
    IMAGE_SPIDER = 1,
    IMAGE_IMAGIC = 2,
    IMAGE_PGM = 4,
    IMAGE_LST = 5,
    IMAGE_PIF = 6,
    IMAGE_DM3 = 7,
    IMAGE_TIF = 8,
    IMAGE_VTK = 9,
    IMAGE_HDF = 10,
    IMAGE_PNG = 11,
    IMAGE_EM = 12,
    IMAGE_SAL,    
    IMAGE_ICOS,
    IMAGE_IMAGIC3D,
    IMAGE_EMIM,
    IMAGE_SINGLE_SPIDER,
    IMAGE_GATAN,
    IMAGE_AMIRA,
    IMAGE_XPLOR,
    IMAGE_SPIDER_SWAP,
    IMAGE_DF3,
    IMAGE_ANY
};

#endif

