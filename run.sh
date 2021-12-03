make 
./main --input ../sample_bin2_filter_input.lst --angpix 2.08 --template ../sample_step5_p1.hdf --eulerfile ../global_euler_step5_p1.lst --phistep 2 --weight ../wiener_reo.txt --kk 3 --first 0 --last 1 --energy 300 --cs 2.7 --Highres 8 --Lowres 100 --diameter 180 --threshold 7 --output test_out.lst --padding_size 256 --device 8
#./main --input ../test/test.lst --angpix 1.688 --template ../test/emd_9976_apix3p336_proj.hdf --eulerfile  ../test/proj_step3_c2.lst --phistep 2 --weight ../test/snr_ncov_search.txt --kk 3 --first 0 --last 1 --energy 300 --cs 2.7 --Highres 8 --Lowres 100 --diameter 180 --threshold 7 --output  test_out.lst --padding_size 387 --device 7

