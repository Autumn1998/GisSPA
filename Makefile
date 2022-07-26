# replace -Lpath -Ipath to hdf5 path.
# path in Lpath: lib in hdf5
# path in Ipath: include in hdf5
LIB_HDF5 = /home/liutong/software/hdf5/lib
INCLUDE_HDF5 = /home/liutong/software/hdf5/include
LINK_FLAG=-lcufft -L$(LIB_HDF5) -lhdf5 
COMP_FLAG=-std=c++11 -I$(INCLUDE_HDF5) -O3 -DDEBUG #-g -G 

OBJ_FILE=./Output/Objects
CC=nvcc 

OBJ=$(OBJ_FILE)/main.o $(OBJ_FILE)/emdata_.o $(OBJ_FILE)/util_func.o $(OBJ_FILE)/emhdf_.o $(OBJ_FILE)/emhdf2_.o $(OBJ_FILE)/DataReader.o $(OBJ_FILE)/GPU_func.o
SRC=emdata_.h MyEmhdf.h MyUtil.h DataType.h
main: $(OBJ) 
	$(CC) $(LINK_FLAG) -o main $(OBJ)
$(OBJ_FILE)/main.o:main.cu 
	$(CC) $(COMP_FLAG) -c main.cu -o $(OBJ_FILE)/main.o
$(OBJ_FILE)/GPU_func.o:GPU_func.cu
	$(CC) $(COMP_FLAG) -c GPU_func.cu -o $(OBJ_FILE)/GPU_func.o
$(OBJ_FILE)/emdata_.o:EMReader/emdata_.cpp
	$(CC) $(COMP_FLAG) -c EMReader/emdata_.cpp -o $(OBJ_FILE)/emdata_.o
$(OBJ_FILE)/util_func.o:EMReader/util_func.cpp
	$(CC) $(COMP_FLAG) -c EMReader/util_func.cpp -o $(OBJ_FILE)/util_func.o
$(OBJ_FILE)/emhdf_.o:EMReader/emhdf_.cpp
	$(CC) $(COMP_FLAG) -c EMReader/emhdf_.cpp -o $(OBJ_FILE)/emhdf_.o
$(OBJ_FILE)/emhdf2_.o:EMReader/emhdf2_.cpp
	$(CC) $(COMP_FLAG) -c EMReader/emhdf2_.cpp -o $(OBJ_FILE)/emhdf2_.o
$(OBJ_FILE)/DataReader.o:EMReader/DataReader.cpp
	$(CC) $(COMP_FLAG) -c EMReader/DataReader.cpp -o $(OBJ_FILE)/DataReader.o
.PHONY : clean
clean:
	rm $(OBJ)