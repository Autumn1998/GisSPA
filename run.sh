make 

# Test Data 1
# Template:150*150*200     Raw Image:512*512
if [ $1 == '0' ]; then
./main  --input Data/Data_1/sample_bin2_filter_input.lst \
        --angpix 2.08 \
        --template Data/Data_1/sample_step5_p1.hdf \
        --eulerfile Data/Data_1/global_euler_step5_p1.lst \
        --phistep 2 \
        --weight Data/Data_1/wiener_reo.txt \
        --kk 3 \
        --first 0 \
        --last 1 \
        --energy 300 \
        --cs 2.7 \
        --Highres 8 \
        --Lowres 100 \
        --diameter 180 \
        --threshold 7 \
        --output Output/test_out.lst \
        --padding_size 256 \
        --device 0 \ 
        --phase_flip 1

# Test Data 2
# Template:182*182*2295    Raw Image:5760*4092
elif [ $1 == '1' ]; then
./main  --input Data/Data_2/test.lst \
        --angpix 1.688 \
        --template Data/Data_2/emd_9976_apix3p336_proj.hdf \
        --eulerfile  Data/Data_2/proj_step3_c2.lst \
        --phistep 2 \
        --weight Data/Data_2/snr_ncov_search.txt \
        --kk 3 \
        --first 0 \
        --last 1 \
        --energy 300 \
        --cs 2.7 \
        --Highres 8 \
        --Lowres 100 \
        --diameter 180 \
        --threshold 7 \
        --output  Output/test_out.lst \
        --padding_size 320 \
        --device 0 \ 
        --phase_flip 1

# Test  performance 
elif [ $1 == '2' ]; then
sudo nvprof ./main --input Data/Data_2/test.lst \
        --angpix 1.688 \
        --template Data/Data_2/emd_9976_apix3p336_proj.hdf \
        --eulerfile  Data/Data_2/proj_step3_c2.lst \
        --phistep 2 \
        --weight Data/Data_2/snr_ncov_search.txt \
        --kk 3 \
        --first 0 \
        --last 1 \
        --energy 300 \
        --cs 2.7 \
        --Highres 8 \
        --Lowres 100 \
        --diameter 180 \
        --threshold 7 \
        --output  Output/test_out.lst \
        --padding_size 320 \
        --device 0 \ 
        --phase_flip 1
fi
