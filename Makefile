CC=nvcc 
LINK_FLAG=-lcufft -L/home/liutong/software/hdf5/lib -lhdf5 
COMP_FLAG=-I/home/liutong/software/hdf5/include -O3 #-DDEBUG #-g -G 
OBJ=main.o emdata_.o util_func.o emhdf_.o emhdf2_.o DataReader.o GPU_func.o
SRC=emdata_.h MyEmhdf.h MyUtil.h DataType.h
main: $(OBJ) 
	$(CC) $(LINK_FLAG) -o main $(OBJ)
main.o:main.cu 
	$(CC) $(COMP_FLAG) -c main.cu 
GPU_func.o:GPU_func.cu
	$(CC) $(COMP_FLAG) -c GPU_func.cu
emdata_.o:EMReader/emdata_.cpp
	$(CC) $(COMP_FLAG) -c EMReader/emdata_.cpp
util_func.o:EMReader/util_func.cpp
	$(CC) $(COMP_FLAG) -c EMReader/util_func.cpp
emhdf_.o:EMReader/emhdf_.cpp
	$(CC) $(COMP_FLAG) -c EMReader/emhdf_.cpp
emhdf2_.o:EMReader/emhdf2_.cpp
	$(CC) $(COMP_FLAG) -c EMReader/emhdf2_.cpp
DataReader.o:EMReader/DataReader.cpp
	$(CC) $(COMP_FLAG) -c EMReader/DataReader.cpp
.PHONY : clean
clean:
	rm $(OBJ)