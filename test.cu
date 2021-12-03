#include "emdata_.h"

int main()
{

    emdata *test = new emdata();
    test->readImage("/home/liutong/chengji/test/emd_9976_apix3p336_proj.hdf",0);
    printf("%d %d %d\n",test->header.nx,test->header.ny,test->header.nz);
    printf("%f %f %f\n",test->header.dx,test->header.dy,test->header.dz);
    float *data = test->getData();
    int cnt = 0;
    for(int i=0;i<10000;i++) {
        if(data[i] >0){
            printf("%d : %f \n",i,data[i]);
            cnt ++;
        }
        if(cnt == 10) break;
    }

    return 0;
}