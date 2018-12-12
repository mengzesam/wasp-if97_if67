#include <stdio.h>
#include "deboorbispline.h"
#define MAX 160000

int main(){
    float x[MAX];
    float y[MAX];
    float z[MAX]; 
    int m=0;
    FILE* fp=fopen("fitdata.txt","r");
    while(fscanf(fp,"%f\t%f\t%f\n",&x[m],&y[m],&z[m])==3 && m<MAX){
        m++;
    }
    fclose(fp);
    free(fp);
    if(m<=0) return -1;
    int nx,ny,nc;
    float* tx;
    float* ty;
    float* c;
    int ier=csurfit(x,y,z,m,&tx,&nx,&ty,&ny,&c,&nc);
    if(ier>0 || ier<-2) return -1;
    printf("double x[%d]={\n\t",nx);
    for(int i=0;i<nx-1;i++){
        printf("%f,",tx[i]);
    };
    printf("%f\n};\n",tx[nx-1]);
    printf("double y[%d]={\n\t",ny);
    for(int i=0;i<ny-1;i++){
        printf("%f,",ty[i]);
    };
    printf("%f\n};\n",ty[ny-1]);
    printf("double c[%d]={\n\t",nc);
    for(int i=0;i<nc-1;i++){
        printf("%f,",c[i]);
    };
    printf("%f\n};\n",c[nc-1]);
    float xx=-1.5042;
    float yy=-1.4471;
    int errFlag;
    float zz,err;
    float maxZ,errMax=0;
    int maxI;
    for(int i=0;i<m;i++){
        zz=cbisp(tx,nx,ty,ny,c,nc,3,3,x[i],y[i],&errFlag);
        err=fabs(zz-z[i]);
        if(err>errMax) {
            errMax=err;
            maxZ=zz;
            maxI=i;
        }
        printf("zz:%15.10f,z:%15.10f,%15.10f\n",zz,z[i],err);
    }
    printf("max error:%15.10f,%15.10f,%15.10f,%15.10f,,%15.10f",
            x[maxI],y[maxI],maxZ,z[maxI],fabs(maxZ-z[maxI]));
    free(tx);
    free(ty);
    free(c);
    return 0;
}