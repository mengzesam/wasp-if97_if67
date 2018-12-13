/*
    reference book:IAPWS R7-97(2012)
*/
#include <math.h>
class IF97Base{
protected://构造函数为protected，禁止非子类实例化
    IF97Base(){};
public: //static function,region4 Saturation-Pressure/temperature Equation
    static double P2T(double ps);
    static double T2P(double ts);
public:
	static double T2P_B23(double t);
	static double P2T_B23(double p);
//bivariate b-spline fitting function
private:
	static void bfit_cbspl(double t[],int n,int k,double x,int idx,double h[6]);
public:
	static double bfit_cbisp(double tx[],int nx,double ty[],int ny,double c[],int nc,int kx,int ky,
            double x,double y,int& errFlag);
protected:
    const static double ERR0;
    const static double ERR;
    const static double ERR2;
    const static double T0;
    const static double R;
    const static double Tc;
    const static double pc;
    const static double rhoc;
private:
    const static double ni[10];//region4 系数,page34，table34
};
const double IF97Base::ERR0=1E-5;
const double IF97Base::ERR=1E-7;
const double IF97Base::ERR2=1E-10;
const double IF97Base::T0=273.15;
const double IF97Base::R=0.461526;//page5 (1)
const double IF97Base::Tc=647.096;//page5 (2)
const double IF97Base::pc=22.064;//page5 (3)
const double IF97Base::rhoc=322;//page5 (4)
const double IF97Base::ni[10]={ //page34 table34
			 0.11670521452767E4,
			-0.72421316703206E6,
			-0.17073846940092E2,
			 0.12020824702470E5,
			-0.32325550322333E7,
			 0.14915108613530E2,
			-0.48232657361591E4,
			 0.40511340542057E6,
			-0.23855557567849E0,
			 0.65017534844798E3
 };
double IF97Base::P2T(double ps){
		double Pstar=1; //p* 1Mpa
		double Tstar=1; //T* 1K
		double beta=pow(ps/Pstar,0.25);
		double E=beta*beta+ni[2]*beta+ni[5];
		double F=ni[0]*beta*beta+ni[3]*beta+ni[6];
		double G=ni[1]*beta*beta+ni[4]*beta+ni[7];
		double D=2*G/(-F-sqrt(F*F-4*E*G));
		double Ts=(ni[9]+D-sqrt(pow(ni[9]+D,2)-4*(ni[8]+ni[9]*D)))/2;
		double ts=Ts*Tstar-273.15;
		return ts;
}
double IF97Base::T2P(double ts){
        double Pstar=1; //p* 1Mpa
		double Tstar=1; //T* 1K
		double theta=(ts+273.15)/Tstar+ni[8]/((ts+273.15)/Tstar-ni[9]);
		double A=theta*theta+ni[0]*theta+ni[1];
		double B=ni[2]*theta*theta+ni[3]*theta+ni[4];
		double C=ni[5]*theta*theta+ni[6]*theta+ni[7];
		double ps=pow(2*C/(-B+sqrt(B*B-4*A*C)),4);
		ps*=Pstar;
		return ps;
}
double IF97Base::T2P_B23(double t){
	double ni[3]={
		0.34805185628969E3,
		-0.11671859879975E1,
		0.10192970039326E-2
	};
	double theta=(t+T0)/1.0;		
	double pi=ni[0]+ni[1]*theta+ni[2]*theta*theta;
	return 1.0*pi;
}
double IF97Base::P2T_B23(double p){
	double ni[3]={
		0.10192970039326E-2,
		0.57254459862746E3,
		0.13918839778870E2
	};	
	double pi=p/1.0;
	double theta=ni[1]+sqrt((pi-ni[2])/ni[0]);
	return 1.0*theta-T0;
}

void IF97Base::bfit_cbspl(double t[],int n,int k,double x,int idx,double h[6]){
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

double IF97Base::bfit_cbisp(double tx[],int nx,double ty[],int ny,double c[],int nc,int kx,int ky,
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
    bfit_cbspl(tx,nx,kx,arg,idx,hx);
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
    bfit_cbspl(ty,ny,ky,arg,idx,hy);
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