# isspa_GPU_forCJ

#### 介绍
GPU parallel for isspa, pick particle with CC.  
CPU version from CJ, BeiJing

#### 软件架构
./deprated_code/:ignored this file   
./EMReader/: read .hdf and .mrc file(just read data). Extracted from EMAN1.9.  
./Get_sigma_from_template/: ignored this file, to read Sigma from .hdf  
./main => main function   
./GPU_func.cu(h) => function processed on GPU  
./Makefile => complier  
./run.sh => example for uesd  

#### 安装教程

1.  install hdf5
2.  complie with Makefile. HDF5 denpendcy is needed to read .hdf
3.  ./main + args

#### 参数说明
This program detect targets with orientations and tanslations.  
--input            :input filtered images lstfile with ctf information  
--angpix           :input pixel size in angstroms  
--template         :input 2D projections templates in .hdf format  
--eulerfile        :euler file with euler values  
--phistep          :inplane rotation sampling  
--weight           :optimized weight with SSNR parmeters.  
--kk               :overlapping density parameter, default is 3.  
--first            :the first image to process.  
--last             :the last image to process.  
--energy           :accerlerating voltage in kV.  
--cs               :spherical aberration in um.  
--Highres          :high resolution cut   
--Lowres           :low resolution cut  
--diameter         :target diameter in pixel  
--threshold        :cc threshold value, only output LOCs beyond this value  
--output           :output lstfile filaname  

To run on GPU  
--padding_size     :size of padded template and raw img will be splitted into  
--device           :which GPU will be used   
注意：  
1.padding_size必须大于template的边长，小于raw img短边边长，并且为32的倍数，图像会被分割成为padding_size*padding_size1的子区域，overlap大小为padding_size的13%  
2.padding_size过小会使并行度不足，过大会提高申请存储的时间和overlap的消耗。在显存足够的情况下，padding_size设置为320左右效果较好。但是仍需根据实际实验结果进行调试(数据量大的时候适当增大).  
3.如果不想使程序打印信息，请在Makefile中删除-DDEBUG  
4.使用的时候,--input的文件所在的文件夹与其中包含的.MRC/.hdf文件相同,即程序在解析raw image的时候，文件所在的文件夹(前缀)与--input的输入相同

#### 参与贡献

LT & CJ.

