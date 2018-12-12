#include "deboorbispline.h"
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
int csurfit(float x[],float y[],float z[],int m,float** tx,int* nx,
            float** ty,int* ny,float** c,int* nc){
    /*
    calling csurfit function such as:
    int ier=csurfit(x,y,z,m,&tx,&nx,&ty,&ny,&c,&nc);
    calling process is respone for free tx,ty and c
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

//cbisp and cbspl of double type 
void cbspl_db(double t[],int n,int k,double x,int idx,double h[6]){
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
double cbisp_db(double tx[],int nx,double ty[],int ny,double c[],int nc,int kx,int ky,
            double x,double y,int* errFlag){
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
    *errFlag=0;
    double tb=tx[kx];
    double te=tx[nx-kx-1];
    double arg=x;
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
    double hx[6];
    int xidx=idx;
    cbspl_db(tx,nx,kx,arg,idx,hx);
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
    double hy[6];
    int yidx=idx;
    cbspl_db(ty,ny,ky,arg,idx,hy);
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