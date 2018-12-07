#include <math.h>
#include <iostream>
using namespace std;

/*
*cbspl and cbisp function are bivariate b-spline functions, they base on fpbspl.f and fpbisp.f
*which were programmed by Prof. P. Dierckx and are re-written with cpp.
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

void cbspl(double t[],int n,int k,double x,int idx,double h[6]){
/*
*as per: FITPACK fpbspl.f,which was coded with fortran,this function translates it into cpp  
*cbspl evaluates the (k+1) non-zero b-splines of degree k at t(idx) <= x < t(idx+1) using 
*the stable recurrence relation of de boor and cox.
*inputs:t,n,k,x,idx
*output:h[6]
*NOTE:k<=5
*/
    double one=1.0;
    h[0]=one;
    double hh[6];
    for(int j=0;j<k;j++){
        for(int i=0;i<=j;i++){
            hh[i]=h[i];
        }
        h[0]=0.0;
        for(int i=0;i<=j;i++){
            int idx1=idx+(i+1);
            int idx0=idx1-(j+1);
            double f=hh[i]/(t[idx1]-t[idx0]);
            h[i]+=f*(t[idx1]-x);
            h[i+1]=f*(x-t[idx0]);
        }
    }
}
double cbisp(double tx[],int nx,double ty[],int ny,double c[],int nc,int kx,int ky,
            double x,double y,int& errFlag){
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
    double z=0.0;
    errFlag=0;
    double tb=tx[kx];
    double te=tx[nx-kx-1];
    double arg=x;
    int idx=kx;
    if(arg<tb){
        errFlag=-1;
        return z;
    }else if(arg>te){
        errFlag=-1;
        return z;
    }else if(abs(arg-te)<1E-7){//arg==te
        idx=(nx-kx-2);
    }else
        while(arg>=tx[idx+1] && idx<(nx-kx-2))
            idx++;
    double hx[6];
    int xidx=idx;
    cbspl(tx,nx,kx,arg,idx,hx);
    tb=ty[ky];
    te=ty[ny-ky-1];
    arg=y;
    idx=ky;
    if(arg<tb){
        errFlag=-1;
        return z;
    }else if(arg>te){
        errFlag=-1;
        return z;
    }else if(abs(arg-te)<1E-7){//arg==te
        idx=(ny-ky-2);
    }else
        while(arg>=ty[idx+1] && idx<(ny-ky-2))
            idx++;
    double hy[6];
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
    /*bispline tck of TH2P for IF97 region2,p=60-80MPa in region2,*/
    double tx[8]={512.0181309,512.0181309,512.0181309,512.0181309,800.0,800.0,800.0,800.0};
    double ty[8]={2658.55376,2658.55376,2658.55376,2658.55376,
                3880.153938,3880.153938, 3880.153938, 3880.153938
                };
    double c[16]={59.98234252,28.23100873,48.73665872, -294.99913313,
                113.06313626,25.62157911,136.07160665,-2.88146108,
                264.11921208,-32.92174268,39.55665078,-6.28943117,
                985.40016598,278.94076817,150.03702759,60.01468942
                };
    int dx=3;
    int dy=3;
    int errFlag=0;
    double x,y,z,zz;

    x=800.0;	
    y=3814.829707;
    zz=74.879;	
    z=cbisp(tx,8,ty,8,c,16,dx,dy,x,y,errFlag);	
    cout<<z<<"\t"<<zz<<"\t"<<abs(zz-z)<<endl;

    return 0;
}