# GisSPA

#### Introduction
GPU parallel for isspa, pick particle with CC.  
CPU version from CJ, BeiJing

#### structure
./deprated_code/:ignored this file   
./EMReader/: read .hdf and .mrc file(just read data). Extracted from EMAN1.9.   
./main.cu: main function   
./GPU_func.cu(h): function processed on GPU  
./Makefile: complier  
./hdf5: the source package of hdf5  
temp file:  
./Data: test data  
./Output: res for test data and obj files.  
./Config_example: configure file to run test data.  

#### Install

0.  install git lfs to get the test_Data
1.  install hdf5 1.8 (https://portal.hdfgroup.org/display/support/HDF5+1.8.21), the source is attached at ./hdf5
2.  Fill the "LIB_HDF5" and "INCLUDE_HDF5" (row 4 & 5) in Makefile with your install path to set hdf5 available.
3.  complie with Makefile. 
4.  write a config file
5.  ./main + config_file

#### Quick start  

-> inatll git lfs  

1. sudo apt-get update  
2. sudo apt-get install git-lfs  
3. git lfs install  
  
-> install hdf5     

1. cd ./hdf5 
2. tar zvxf hdf5-1.8.21.tar.gz     
3. cd hdf5-1.8.21   
4. ./configure --prefix="your hdf5 install path"   
5. make   
6. make install   
7. export LD_LIBRARY_PATH="your hdf5 install path"/lib:$LD_LIBRARY_PATH  
   at /etc/profile or ~/.bashrc    
8. source /etc/profile or ~/.bashrc   
  
-> run GisSPA   
   
1. vim Makfile, set LIB_HDF5="your hdf5 install path"/lib,  set INCLUDE_HDF5="your hdf5 install path"/include  
2. make clean  
3. make  
	If "can not create /Output/Objects/main.o", mkdir /Output/Objects  
4. ./main ./Data_TEST/config   

WARNING:The Data_TEST need at least 13GB memory at GPU!     
If temp_size = 10,10 and img_size = 10,10, check the git lfs.  
   
The result will be find at ./Output/test_Image_bin2_output.lst. We attacted the result at ./Data_TEST/test_Image_bin2_output.lst.   
  
If "error while loading shared libraries", check your LD_LIBRARY_PATH, and source. If it not work, do  
	cp "your hdf5 install path"/lib/libhdf5.so.10 .  
Make ./main  and libhdf5.so.10 in the same path.   

#### Parameters
This program detect targets with orientations and tanslations.

How to use this problem:  "./main config_file"  
use "./main -h" to get help msg

<br />--------------------------------config example:--------------------------------------<br />
\# Do leave spaces on either side of '='  
\# Parameter request  
input     = Data_TEST/test_Image.lst 
template  = Data_TEST/emd_9976_apix3p336_proj.hdf  
eulerfile = Data_TEST/proj_step3_c2.lst  
angpix    = 3.336  
phistep   = 2   
kk        = 3   
energy    = 300  
cs        = 2.7  
Highres   = 9  
Lowres    = 100  
diameter  = 180  
  
\# Parameter optional  
threshold = 8  
output    = Data_TEST/test_Image_output.lst 
first     = 0  
last      = 12  
window_size = 480  
phase_flip  = 1  
GPU_ID      = 0  
overlap     = 180  
<br />--------------------------------------------------------------------------------------<br />

HDF5 and CUDA lib are needed to run  
All parameters should be set on the config file.  
All parameters are listed as below which can be set at config file.  
(for .eg,  "input  =  /Data/inputfile" )Requested Parameters:  
input            = input contast-inverted images lstfile with ctf information  
template         = input 2D projections templates in .hdf format  
eulerfile        = euler file with euler values  
angpix           = input pixel size in angstroms  
phistep          = inplane rotation sampling  
kk               = overlapping density parameter, default is 3.  
energy           = accerlerating voltage in kV.  
cs               = spherical aberration in um.  
Highres          = high resolution cut   
Lowres           = low resolution cut  
diameter         = target diameter in pixels 
Optional Parameters:  
threshold        = cc threshold value, only output score beyond this value  
output           = output lstfile filaname  
first            = the first image id to process.  
last             = the last image id to process.  
window_size      = the window size which is splitted from raw IMG.  
GPU_ID           = ID of GPU device.   
phase_flip       = Whether do filtering on images, operation(1) or not(0, in case of input being filtered already).  
overlap          = size of overlap between diff window.   

attention:  
1. should be:  
template size: tx,ty  
img size:ix.iy  
a)max(tx,ty) < padding_size < max(ix,iy)  
b)padding_size%32 = 0  
2. as padding_size increasing, consumed memory decreases↓ first, then increases↑.  
3. for parameter "input". prefix will be auto-added. prefix is current dir  
4. all imgs should have same size.  
<br />   
注意：    
1.padding_size必须大于template的边长，小于raw img短边边长，并且为32的倍数，图像会被分割成为padding_size*padding_size1的子区域，overlap大小为padding_size的13%.   <br />   
2.padding_size过小会使并行度不足，过大会提高申请存储的时间和overlap的消耗。在显存足够的情况下，padding_size设置为320左右效果较好。但是仍需根据实际实验结果进行调试(数据量大的时候适当增大).    <br />
4.使用的时候,--input的文件所在的文件夹与其中包含的.MRC/.hdf文件相同,即程序在解析raw image的时候，文件所在的文件夹(前缀)与--input的输入相同.  <br />
5.所有img的size必须相同     <br />

#### Python script

relion2lst_only_particles.py (rewritten from Wen Jiang's script from JSPR)   
This script read particles.star file and generate file in .lst format and rewrite images in .hdf format.   
Please run ./relion2lst_only_particles.py -h for details.  

remove_repeat_particles_from_list.py   
<Detection file> <number of windows> <center thres> <euler thres> <output>
This python script is used to merge duplicate detections:  
<Dection file> = Data_TEST/test_Image_bin2_output.lst  
<number of windows> = last-first (eg. 12)  
<center thres> = Threshold of coordinate distance (eg. 4)  
<euler thres> = Threshold of euler distance (eg. 6)    
<outpot> = Merged list filename (eg. Data_TEST/test_Image_bin2_output_merged.lst)    
Detections that are within both the center and euler thresholds are considered to be duplicate detections   

convert_my-output_to_relion.py  
<Merged lstfile> <contast-inverted lstfile> <scale factor> <scaled window size> <unbinned pixel size> <out star file>  
<Merged lstfile> = detection file after removing duplicate detections (eg. Data_TEST/test_Image_output_merged.lst)  
<contast-inverted lstfile> = contast-inverted images lstfile with ctf information (eg. test_Image.lst)  
<scale factor> = scale factor of images used in localization (eg. 2)  
<scaled window size> = window size in localization (eg. 720)  
<unbinned pixel size> = pixel size in original micrograph (eg. 1.668)  
<out star file> = output file in .star format  

A complete workflow of the demo data can be found at folder Data_TEST/workflow.txt
##Micrographs are clipped into multiple images to save GPU memory. In case of sufficient GPU memory, this step can be ignored, but the contast of micrographs must be inverted.

#### contributor  

LT & CJ.

