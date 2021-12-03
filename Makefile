CC=nvcc 
LINK_FLAG=-lcufft -Xcompiler -fopenmp -pg -L/home/liutong/software/hdf5/lib -lhdf5 
COMP_FLAG=-I/home/liutong/software/hdf5/include -O3 -DH5Gopen_vers=1 -DDEBUG -g -G
OBJ=main.o emdata_.o util_func.o emhdf_.o emhdf2_.o DataReader.o GPU_func.o
SRC=emdata_.h MyEmhdf.h MyUtil.h DataType.h
main: $(OBJ) 
	$(CC) $(LINK_FLAG) -o main $(OBJ)
main.o:main.cu 
	$(CC) $(COMP_FLAG) -c main.cu 
emdata_.o:emdata_.cpp
	$(CC) $(COMP_FLAG) -c emdata_.cpp
util_func.o:util_func.cpp
	$(CC) $(COMP_FLAG) -c util_func.cpp
emhdf_.o:emhdf_.cpp
	$(CC) $(COMP_FLAG) -c emhdf_.cpp
emhdf2_.o:emhdf2_.cpp
	$(CC) $(COMP_FLAG) -c emhdf2_.cpp
DataReader.o:DataReader.cpp
	$(CC) $(COMP_FLAG) -c DataReader.cpp
GPU_func.o:GPU_func.cu
	$(CC) $(COMP_FLAG) -c GPU_func.cu
.PHONY : clean
clean:
	rm $(OBJ)