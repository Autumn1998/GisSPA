# GisSPA
 
  
Time for compling: short     
No need to install  

#### Introduction
GPU parallel for isSPA, pick particle with CC.  
CPU version from CJ, BeiJing

######################################   
This site is no longer updated, please check "https://www.x-mol.com/groups/Zhang_Xinzheng" for subsequent updates.  
Please find the test data on this website.  


#### structure
./deprated_code/:ignored this file   
./EMReader/: read .hdf and .mrc file(just read data). Extracted from EMAN1.9.   
./main.cu: main function   
./GPU_func.cu(h): function processed on GPU  
./Makefile: complier  
./hdf5: the source package of hdf5  
temp file:  
./Data_TEST: test data  
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
4. git lfs clone https://github.com/Autumn1998/GisSPA.git   
  
-> install hdf5     

0. cd GisSPA  
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
   
0. cd GisSPA  
1. vim Makfile, set LIB_HDF5="your hdf5 install path"/lib,  set INCLUDE_HDF5="your hdf5 install path"/include  
2. make clean  
3. make  
	If "can not create /Output/Objects/main.o", mkdir /Output/Objects  
4. ./main ./Data_TEST/config   

WARNING:The Data_TEST need at least 13GB memory at GPU!     
If error occurred and temp_size = 10,10 and img_size = 10,10, check the git lfs.  
If "error while loading shared libraries", check your LD_LIBRARY_PATH, and source. If it not work, do  
	cp "your hdf5 install path"/lib/libhdf5.so.10 .  
Make ./main  and libhdf5.so.10 in the same path.   

  

How to use this problem:  "./main config_file"  
use "./main -h" to get help msg


#### contributor  

LT & CJ.

