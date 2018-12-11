#include <stdlib.h>
#include <math.h>
#include <stdio.h>
/*
* cbspl cbspl and cbisp function are bivariate b-spline functions, they base on fpbspl.f and fpbisp.f
*which were programmed by Prof. P. Dierckx and are re-written with cpp.
*surfit_ is c interface function of surfit.f
*fpbspl.f and fpbisp.f author:
        Prof. P. Dierckx
		Department of Computer Science
		K.U.Leuven
		Celestijnenlaan 200 A
		B 3001, Heverlee, Belgium fb are re-written with cpp
*please refer to http://www.netlib.org/dierckx/index.html for details
*Re-written by:
        mengzesen<mengzesam@126.com>
*/
extern void surfit_(int* iopt,int* m,float x[],float y[],float z[],float w[],float* xb,float* xe,
                     float* yb,float* ye,int* kx,int* ky,float* s,int* nxest,int* nyest,int* nmax,
                     float* eps,int* nx,float tx[],int* ny,float ty[],float c[],float* fp,
                     float wrk1[],int* lwrk1,float wrk2[],int* lwrk2,int iwrk[],int* kwrk,int* ier);

int csurfit(float x[],float y[],float z[],int m,float** tx,int* nx,
            float** ty,int* ny,float** c,int* nc){
    /*
    calling csurfit function such as:
    int ier=csurfit(x,y,z,m,&tx,&nx,&ty,&ny,&c,&nc);
    */
        int iopt=0;        
        int kx=3;
        int ky=3;
        if(m<(kx+1)*(ky+1)){
            //printf("m >= (kx+1)(ky+1) must hold");
            return -1;
        }        
        int nxest=kx+1+(int)sqrt(m/2);
        nxest=(nxest>=2*(kx+1))?nxest:2*(kx+1); 
        int nyest=nxest;
        nyest=(nyest>=2*(ky+1))?nyest:2*(ky+1); 
        int nmax=nxest>=nyest?nxest:nyest;
        float* txx=(float*)malloc(sizeof(float)*nmax);  
        float* tyy=(float*)malloc(sizeof(float)*nmax);  
        float* cc=(float*)malloc(sizeof(float)*((nxest-kx-1)*(nyest-ky-1))); 
        float s=m-sqrt(2*m);
        float eps=1E-7;
        float fp;
        float* w=(float*)malloc(sizeof(float)*m);        
        float xb,xe,yb,ye;
        xb=x[0];
        xe=x[0];
        yb=y[0];
        ye=y[0];
        w[0]=1.0;
        for(int i=1;i<m;i++){
            w[i]=1.0;
            if(x[i]<xb)
                xb=x[i];
            else if(x[i]>xe)
                xe=x[i];
            if(y[i]<yb)
                yb=y[i];
            else if(y[i]>ye)
                ye=y[i];            
        } 
        int u=nxest-kx-1;
        int v=nyest-ky-1;
        int km=(kx>kx?kx:ky)+1;
        int ne = nmax;
        int bx=kx*v+ky+1;
        int by=ky*u+kx+1;
        int b1=bx;
        int b2=bx+v-ky;
        if(bx>by){
            b1=by;
            b2=by+u-kx;
        }
        int lwrk1=u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+10;
        int lwrk2=u*v*(b2+1)+b2+10;       
        int kwrk=m+(nxest-2*kx-1)*(nyest-2*ky-1)+10;
        /* int pause;
        printf("%d:%d:%d\n",lwrk1,lwrk2,kwrk);
        scanf("%d",&pause); */
        float* wrk1=(float*)malloc(sizeof(float)*lwrk1);  
        float* wrk2=(float*)malloc(sizeof(float)*lwrk2); 
        int* iwrk=(int*)malloc(sizeof(int)*kwrk); 
        int ier;
        //scanf("%d",&pause);
        surfit_(&iopt,&m,x,y,z,w,&xb,&xe,
                &yb,&ye,&kx,&ky,&s,&nxest,&nyest,&nmax,
                &eps,nx,txx,ny,tyy,cc,&fp,
                wrk1,&lwrk1,wrk2,&lwrk2,iwrk,&kwrk,&ier);
        *nc=((*nx)-kx-1)*((*ny)-ky-1);
        *tx=(float*)malloc(sizeof(float)*(*nx));  
        *ty=(float*)malloc(sizeof(float)*(*ny));  
        *c=(float*)malloc(sizeof(float)*(*nc)); 
        for(int i=0;i<(*nx);i++)
            (*tx)[i]=txx[i];
        free(txx);
        for(int i=0;i<(*ny);i++)
            (*ty)[i]=tyy[i];
        free(tyy);
        for(int i=0;i<(*nc);i++)
            (*c)[i]=cc[i];
        free(cc);
        free(w);
        free(wrk1);
        free(wrk2);
        free(iwrk);
        //scanf("%d",&pause);
        return ier;
}
void cbspl(float t[],int n,int k,float x,int idx,float h[6]){
    /*
    *as per: FITPACK fpbspl.f,which was coded with fortran,this function translates it into cpp  
    *cbspl evaluates the (k+1) non-zero b-splines of degree k at t(idx) <= x < t(idx+1) using 
    *the stable recurrence relation of de boor and cox.
    *inputs:t,n,k,x,idx
    *output:h[6]
    *NOTE:k<=5
    */
    float one=1.0;
    h[0]=one;
    float hh[6];
    for(int j=0;j<k;j++){
        for(int i=0;i<=j;i++){
            hh[i]=h[i];
        }
        h[0]=0.0;
        for(int i=0;i<=j;i++){
            int idx1=idx+(i+1);
            int idx0=idx1-(j+1);
            float f=hh[i]/(t[idx1]-t[idx0]);
            h[i]+=f*(t[idx1]-x);
            h[i+1]=f*(x-t[idx0]);
        }
    }
}
float cbisp(float tx[],int nx,float ty[],int ny,float c[],int nc,int kx,int ky,
            float x,float y,int* errFlag){
    /*
    *as per:FITPACK fpbisp.f,which was code with fortran,this function translate it into cpp
    *inputs:tx[nx],ty[ny],c[(nx-kx-1)*(ny-ky-1)],kx,ky bispline knot points,coefficents 
    *and x,y direction degree
    * tx[kx]<tx[kx+1]...<tx[nx-2]
    * ty[ky]<ty[ky+1]...<ty[ny-2]
    *inputs: (x,y)to be evaled points
    *return: evaled value
    *
    *NOTE:kx,ky<=5
    */
    float z=0.0;
    *errFlag=0;
    float tb=tx[kx];
    float te=tx[nx-kx-1];
    float arg=x;
    int idx=kx;
    if(arg<tb){
        *errFlag=-1;
        return z;
    }else if(arg>te){
        *errFlag=-1;
        return z;
    }else if(abs(arg-te)<1E-7){//arg==te
        idx=(nx-kx-2);
    }else
        while(arg>=tx[idx+1] && idx<(nx-kx-2))
            idx++;
    float hx[6];
    int xidx=idx;
    cbspl(tx,nx,kx,arg,idx,hx);
    tb=ty[ky];
    te=ty[ny-ky-1];
    arg=y;
    idx=ky;
    if(arg<tb){
        *errFlag=-1;
        return z;
    }else if(arg>te){
        *errFlag=-1;
        return z;
    }else if(abs(arg-te)<1E-7){//arg==te
        idx=(ny-ky-2);
    }else
        while(arg>=ty[idx+1] && idx<(ny-ky-2))
            idx++;
    float hy[6];
    int yidx=idx;
    cbspl(ty,ny,ky,arg,idx,hy);
    idx=(xidx-kx)*(ky+1)+(yidx-ky);
    for(int i=0;i<kx+1;i++){
        int cidx=idx;
        for(int j=0;j<ky+1;j++){
            z+=c[cidx]*hx[i]*hy[j];
            cidx++;
        }
        idx+=ky+1;
    }
    return z;
}

int main(){
/*     float x[80]={
        -1.9867,-1.5042,-1.7742,-1.6004,-1.5655,-1.3779,-1.9129,-1.5747,-1.6036,-1.953,
        -1.3627,-1.2517,-1.0566,-1.6402,-1.3958,-1.671,-1.1381,-1.7332,-1.792,-1.0052,-0.8194,
        -0.9082,-0.339,-0.0594,-0.3295,-0.7246,-0.3325,-0.7718,-0.6555,-0.1243,-0.0757,-0.7134,
        -0.0422,-0.9322,-0.3798,-0.5879,-0.4919,-0.3558,-0.0638,-0.7318,0.7348,0.5143,0.5456,
        0.3855,0.6197,0.2048,0.1138,0.6618,0.887,0.3598,0.7982,0.8622,0.5027,0.2402,0.7488,
        0.1189,0.1752,0.2245,0.6075,0.4312,1.8552,1.3081,1.0085,1.1658,1.2353,1.6036,1.078,
        1.5563,1.2505,1.8409,1.9347,1.2037,1.5779,1.8671,1.1892,1.5763,1.1368,1.1464,1.4444,1.5092
    };
    float y[80]={
        -1.947,-1.4471,-1.0253,-1.7173,-1.8963,-0.3124,-0.2167,-0.4469,-0.3483,-0.4068,0.5334,0.5905,
        0.8297,0.0373,0.4089,1.2759,1.9973,1.1393,1.0699,1.0919,-1.7594,-1.13,-1.8467,-1.3774,-1.1278,
        -0.4482,-0.2274,-0.9464,-0.0687,-0.6415,0.0176,0.6721,0.6835,0.8158,0.4549,1.3671,1.1213,
        1.3006,1.6276,1.1287,-1.5365,-1.4583,-1.4447,-1.451,-1.3444,-0.8807,-0.7657,-0.6786,-0.1306,
        -0.035,0.1373,0.741,0.0037,0.8623,0.7138,1.2044,1.1834,1.946,1.8895,1.7096,-1.7396,-1.7975,
        -1.4331,-1.6434,-1.6872,-0.0301,-0.8749,-0.1535,-0.8546,-0.5623,0.9197,0.9823,0.7897,0.1472,
        0.7828,1.1683,1.0035,1.4909,1.5232,1.9105
    };
    float z[80]={
        0.0076,-0.0035,0.0144,0.0036,0.0157,0.1309,0.0402,0.0625,0.0553,0.0133,0.1119,0.1269,0.1708,
        0.0638,0.1077,0.0274,-0.0092,0.0169,0.0178,0.1106,0.0154,0.1125,0.0388,0.1521,0.2237,0.4845,
        0.8418,0.2268,0.6497,0.6507,1.0069,0.369,0.6235,0.2226,0.7087,0.1145,0.2194,0.163,0.079,0.161,
        0.0682,0.0947,0.0775,0.1061,0.117,0.4397,0.556,0.393,0.4664,0.8872,0.5129,0.2846,0.7766,0.4501,
        0.3513,0.2338,0.2242,0.0087,0.0106,0.0472,-0.0076,0.008,0.0401,0.0157,0.0259,0.0826,0.1429,
        0.076,0.0995,0.0319,0.0138,0.0934,0.0479,0.0276,0.1353,0.0123,0.1042,0.0203,0.0108,0.0022
    }; 
*/
    float x[110000];
    float y[110000];
    float z[110000];    
    int MAX=110000;
    int m=0;
    FILE* fp=fopen("fitdata.csv","r");
    int r=fscanf(fp,"%f,%f,%f\n",&x[m],&y[m],&z[m]);
    while(r==3 && m<MAX){
        m++;
        r=fscanf(fp,"%f,%f,%f\n",&x[m],&y[m],&z[m]);
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
    float zz;
    for(int i=0;i<m;i++){
        zz=cbisp(tx,nx,ty,ny,c,nc,3,3,x[i],y[i],&errFlag);
        printf("zz:%15.10f,z:%15.10f,%15.10f\n",zz,z[i],fabs(zz-z[i]));
    }
    return 0;
}