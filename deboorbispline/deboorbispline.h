#include <stdlib.h>
#include <math.h>
#include <stdio.h>
/*
* csurfit,cbspl and cbisp function are bivariate b-spline functions, they base on fpbspl.f and fpbisp.f
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
/*
    surfit_ is a  extern function from fortran,pls refer to surfit.f or combined fortran file bisurfit.f
*/
extern void surfit_(int* iopt,int* m,float x[],float y[],float z[],float w[],float* xb,float* xe,
                     float* yb,float* ye,int* kx,int* ky,float* s,int* nxest,int* nyest,int* nmax,
                     float* eps,int* nx,float tx[],int* ny,float ty[],float c[],float* fp,
                     float wrk1[],int* lwrk1,float wrk2[],int* lwrk2,int iwrk[],int* kwrk,int* ier);
/*
    calling csurfit function such as:
    int ier=csurfit(x,y,z,m,&tx,&nx,&ty,&ny,&c,&nc);
    calling process is respone for free tx,ty and c
*/
int csurfit(float x[],float y[],float z[],int m,float** tx,int* nx,
            float** ty,int* ny,float** c,int* nc);
    
/*
    *cbspl is as per: FITPACK fpbspl.f,which was coded with fortran,this function translates it into cpp  
    *cbspl evaluates the (k+1) non-zero b-splines of degree k at t(idx) <= x < t(idx+1) using 
    *the stable recurrence relation of de boor and cox.
    *inputs:t,n,k,x,idx
    *output:h[6]
    *NOTE:k<=5
*/
void cbspl(float t[],int n,int k,float x,int idx,float h[6]);

/*
    *cbisp is as per:FITPACK fpbisp.f,which was code with fortran,this function translate it into cpp
    *inputs:tx[nx],ty[ny],c[(nx-kx-1)*(ny-ky-1)],kx,ky bispline knot points,coefficents 
    *and x,y direction degree
    * tx[kx]<tx[kx+1]...<tx[nx-2]
    * ty[ky]<ty[ky+1]...<ty[ny-2]
    *inputs: (x,y)to be evaled points
    *return: evaled value
    *
    *NOTE:kx,ky<=5
*/   
float cbisp(float tx[],int nx,float ty[],int ny,float c[],int nc,int kx,int ky,
            float x,float y,int* errFlag);
    
   

