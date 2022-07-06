#include "EMData.h"

int main(int argc,char **argv)
{
    int N = atoi(argv[1]);
    char path[1000];
    strcpy(path,argv[2]);
    char out[1000];
    strcpy(out,argv[3]);
    FILE *fp=fopen(out,"wb");
    EMData * reader = new EMData();
    for(int i=0;i<N;i++)
    {
        reader->readImage(path,i);
        fprintf(fp,"%lf\n",reader->Sigma());
    } 
    fclose(fp);
    return 0;
}