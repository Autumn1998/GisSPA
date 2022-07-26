# GisSPA

#### Introduction
GPU parallel for isspa, pick particle with CC.  
CPU version from CJ, BeiJing

#### structure
./deprated_code/:ignored this file   
./EMReader/: read .hdf and .mrc file(just read data). Extracted from EMAN1.9.   
./main => main function   
./GPU_func.cu(h) => function processed on GPU  
./Makefile => complier  
./config => example for config file

#### Install

1.  install hdf5
2.  replace the "LINK_FLAG" and "COMP_FLAG" in Makefile to set hdf5 available.
3.  complie with Makefile. HDF5 denpendcy is needed to read ".hdf"
4.  ./main + args

#### Parameters
This program detect targets with orientations and tanslations.

How to use this problem:  "./main config_file"  
use "./main -h" to get help msg

<br />--------------------------------config example:--------------------------------------<br />
\# Do leave spaces on either side of '='
\# Parameter request
input     = Data/Data_2/test.lst
template  = Data/Data_2/emd_9976_apix3p336_proj.hdf
eulerfile = Data/Data_2/proj_step3_c2.lst
angpix    = 1.688
phistep   = 2 
kk        = 3 
energy    = 300
cs        = 2.7
Highres   = 8
Lowres    = 100
diameter  = 180

\# Parameter optional
threshold = 7
output    = Output/test_out.lst
first     = 0
last      = 1
window_size = 320
phase_flip  = 1
GPU_ID      = 1
overlap     = 24
<br />--------------------------------------------------------------------------------------<br />

HDF5 and CUDA lib are needed to run
All parameters should be set on the config file.
All parameters are listed as below which can be set at config file.
(for .eg,  "input  =  /Data/inputfile" )Requested Parameters:
input            = input filtered images lstfile with ctf information
template         = input 2D projections templates in .hdf format
eulerfile        = euler file with euler values
angpix           = input pixel size in angstroms
phistep          = inplane rotation sampling
kk               = overlapping density parameter, default is 3.
energy           = accerlerating voltage in kV.
cs               = spherical aberration in um.
Highres          = high resolution cut 
Lowres           = low resolution cut
diameter         = target diameter in pixel
Optional Parameters:
threshold        = cc threshold value, only output score beyond this value
output           = output lstfile filaname
first            = the first image id to process.
last             = the last image id to process.
window_size      = the window size which is splitted from raw IMG.
GPU_ID           = ID of GPU device.
phase_flip       = Whether use phase filp operation(1) or not(0).
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
1.padding_size必须大于template的边长，小于raw img短边边长，并且为32的倍数，图像会被分割成为padding_size*padding_size1的子区域，overlap大小为padding_size的13%  
2.padding_size过小会使并行度不足，过大会提高申请存储的时间和overlap的消耗。在显存足够的情况下，padding_size设置为320左右效果较好。但是仍需根据实际实验结果进行调试(数据量大的时候适当增大).  
4.使用的时候,--input的文件所在的文件夹与其中包含的.MRC/.hdf文件相同,即程序在解析raw image的时候，文件所在的文件夹(前缀)与--input的输入相同
5.所有img的size必须相同

#### contributor

LT & CJ.

