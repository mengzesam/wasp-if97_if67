#include "IF97Region3.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <time.h>
using namespace std;

double IF97Region3::TV2P(double t,double v){
    double delta=1.0/(rhoc*v);
    double tau=Tc/(t+T0);
    double phi_delta=ni[39]/delta;
    for(int i =0;i<39;i++){
        phi_delta+=ni[i]*Ii[i]*pow(delta,Ii[i]-1)*pow(tau,Ji[i]);
    }
    return delta*phi_delta*R*(t+T0)/(v*1000.0);
    //in if97--table31 there is not  "/1000"
}
double IF97Region3::TV2H(double t,double v){
    double delta=1.0/(rhoc*v);
    double tau=Tc/(t+T0);
    double phi_delta=ni[39]/delta;
    double phi_tau=0;
    for(int i =0;i<39;i++){
        phi_delta+=ni[i]*Ii[i]*pow(delta,Ii[i]-1)*pow(tau,Ji[i]);
        phi_tau+=ni[i]*pow(delta,Ii[i])*Ji[i]*pow(tau,Ji[i]-1);
    }
	return (tau*phi_tau+delta*phi_delta)*R*(t+T0);
}
double IF97Region3::TV2U(double t,double v){
    double delta=1.0/(rhoc*v);
    double tau=Tc/(t+T0);
    double phi_tau=0;
    for(int i =0;i<39;i++){
        phi_tau+=ni[i]*pow(delta,Ii[i])*Ji[i]*pow(tau,Ji[i]-1);
    }
    return tau*phi_tau*R*(t+T0);
}
double IF97Region3::TV2S(double t,double v){
    double delta=1.0/(rhoc*v);
    double tau=Tc/(t+T0);
    double phi=ni[39]*log(delta);
    double phi_tau=0;
    for(int i =0;i<39;i++){
		phi+=ni[i]*pow(delta,Ii[i])*pow(tau,Ji[i]);
        phi_tau+=ni[i]*pow(delta,Ii[i])*Ji[i]*pow(tau,Ji[i]-1);
    }
	return (tau*phi_tau-phi)*R;
}
double IF97Region3::TV2Cp(double t,double v){
    double delta=1.0/(rhoc*v);
    double tau=Tc/(t+T0);
    double phi_delta=ni[39]/delta;
    double phi_delta2=-ni[39]/(delta*delta);
    double phi_tau2=0;
    double phi_deltatau=0;
    for(int i =0;i<39;i++){
        phi_delta+=ni[i]*Ii[i]*pow(delta,Ii[i]-1)*pow(tau,Ji[i]);
        phi_delta2+=ni[i]*Ii[i]*(Ii[i]-1)*pow(delta,Ii[i]-2)*pow(tau,Ji[i]);
		phi_tau2+=ni[i]*pow(delta,Ii[i])*Ji[i]*(Ji[i]-1)*pow(tau,Ji[i]-2);
		phi_deltatau+=ni[i]*Ii[i]*pow(delta,Ii[i]-1)*Ji[i]*pow(tau,Ji[i]-1);
    }
	double cp=-tau*tau*phi_tau2;
	cp+=pow(delta*phi_delta-delta*tau*phi_deltatau,2)
        /(2*delta*phi_delta+delta*delta*phi_delta2);
    cp*=R;
    return cp;
}
double IF97Region3::TV2Cv(double t,double v){
    double delta=1.0/(rhoc*v);
    double tau=Tc/(t+T0);
    double phi_tau2=0;
    for(int i =0;i<39;i++){
		phi_tau2+=ni[i]*pow(delta,Ii[i])*Ji[i]*(Ji[i]-1)*pow(tau,Ji[i]-2);
    }
	return -tau*tau*phi_tau2*R;
}
double IF97Region3::TV2W(double t,double v){
    double p,h;
    h=H3ab(25);                                        
    t=PH2T3a(20,1700);
    t=PH2T3a(20,1827.10060895);
    t=PH2T3b(20,1827.10060895);
    h=H3ab(20);
    p=H2P3sat(1827.10060895);
    t=PH2T3a(50,2000);
    t=PH2T3a(100,2100);
    t=PH2T3b(20,2500);
    t=PH2T3b(50,2400);
    t=PH2T3b(100,2700);
    v=PH2V3a(20,1700);
    v=PH2V3a(50,2000);
    v=PH2V3a(100,2100);
    v=PH2V3b(20,2500);
    v=PH2V3b(50,2400);
    v=PH2V3b(100,2700);
    t=PS2T3a(20,3.8);
    t=PS2T3a(50,3.6);
    t=PS2T3a(100,4.0);
    t=PS2T3b(20,5.0);
    t=PS2T3b(50,4.5);
    t=PS2T3b(100,5.0);
    v=PS2V3a(20,3.8);
    v=PS2V3a(50,3.6);
    v=PS2V3a(100,4.0);
    v=PS2V3b(20,5.0);
    v=PS2V3b(50,4.5);
    v=PS2V3b(100,5.0);
    p=S2P3sat(3.8);
    p=S2P3sat(4.2);
    p=S2P3sat(5.2);
    p=H2P3sat(1700);
    p=H2P3sat(2000);
    p=H2P3sat(2400);
    return 0.0;
}
/*P,T 2*/
double IF97Region3::PT2V(double p,double t,int& itera){
    double err=ERR;
    double v;
    if(p>40){ //&& p<=100
        if(t<=T3ab(p))
            v=PT2V3a(p,t);
        else
            v=PT2V3b(p,t);
    }else if(p>25){
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3ab(p))
            v=PT2V3d(p,t);
        else if(t<=T3ef(p))
            v=PT2V3e(p,t);
        else
            v=PT2V3f(p,t);
    }else if(p>23.5){
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3gh(p))
            v=PT2V3g(p,t);
        else if(t<=T3ef(p))
            v=PT2V3h(p,t);
        else if(t<=T3ij(p))
            v=PT2V3i(p,t);
        else if(t<=T3jk(p))
            v=PT2V3j(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>23){
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3gh(p))
            v=PT2V3l(p,t);
        else if(t<=T3ef(p))
            v=PT2V3h(p,t);
        else if(t<=T3ij(p))
            v=PT2V3i(p,t);
        else if(t<=T3jk(p))
            v=PT2V3j(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>22.5){
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3gh(p))
            v=PT2V3l(p,t);
        else if(t<=T3mn(p))
            v=PT2V3m(p,t);
        else if(t<=T3ef(p))
            v=PT2V3n(p,t);
        else if(t<=T3op(p))
            v=PT2V3o(p,t);
        else if(t<=T3ij(p))
            v=PT2V3p(p,t);
        else if(t<=T3jk(p))
            v=PT2V3j(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>22.11){
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3qu(p))
            v=PT2V3q(p,t);
        else if(t<=T3uv(p))
            v=PT2V3u(p,t);
        else if(t<=T3ef(p))
            v=PT2V3v(p,t);
        else if(t<=T3wx(p))
            v=PT2V3w(p,t);
        else if(t<=T3rx(p))
            v=PT2V3x(p,t);
        else if(t<=T3jk(p))
            v=PT2V3r(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>22.064){
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3qu(p))
            v=PT2V3q(p,t);
        else if(t<=T3uv(p))
            v=PT2V3u(p,t);
        else if(t<=T3ef(p))
            v=PT2V3y(p,t);
        else if(t<=T3wx(p))
            v=PT2V3z(p,t);
        else if(t<=T3rx(p))
            v=PT2V3x(p,t);
        else if(t<=T3jk(p))
            v=PT2V3r(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>21.93161551){//21.93161551==psat97(v=0.00264m3/kg)f
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3qu(p))
            v=PT2V3q(p,t);
        else if(t<=T3uv(p))
            v=PT2V3u(p,t);
        else if(t<=P2T(p))//<Tsat97(p)
            v=PT2V3y(p,t);         
        else if(t<=T3wx(p))
            v=PT2V3z(p,t);         
        else if(t<=T3rx(p))
            v=PT2V3x(p,t);            
        else if(t<=T3jk(p))
            v=PT2V3r(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>21.90096265){//21.90096265==psat97(v=0.00385m3/kg)f
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3qu(p))
            v=PT2V3q(p,t);
        else if(t<=P2T(p))//<Tsat97(p)
            v=PT2V3u(p,t);       
        else if(t<=T3wx(p))
            v=PT2V3z(p,t);         
        else if(t<=T3rx(p))
            v=PT2V3x(p,t);            
        else if(t<=T3jk(p))
            v=PT2V3r(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>21.04336731897525){//p>psat97(643.15K)
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3qu(p))
            v=PT2V3q(p,t);
        else if(t<=P2T(p))//<Tsat97(p)
            v=PT2V3u(p,t);         
        else if(t<=T3rx(p))
            v=PT2V3x(p,t);            
        else if(t<=T3jk(p))
            v=PT2V3r(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>20.5){
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=P2T(p))//t<=Tsat97(p)
            v=PT2V3s(p,t);
        else if(t<=T3jk(p))
            v=PT2V3r(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>19.00881189173929){//p>p3cd
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=P2T(p))//t<=Tsat97(p)
            v=PT2V3s(p,t);
        else
            v=PT2V3t(p,t);
    }else{//p>Psat97(623.15K)
        if(t<=P2T(p))//t<=Tsat97(p)
            v=PT2V3c(p,t);
        else
            v=PT2V3t(p,t);
    }
    double tau=Tc/(t+T0);
    double lefts[39];
    for(int i=0;i<39;i++){
        lefts[i]=ni[i]*Ii[i]*pow(tau,Ji[i]);
    }
    double delta=1.0/(rhoc*v);
    double phi_delta=ni[39]/delta;
    for(int i=0;i<39;i++){
        phi_delta+=lefts[i]*pow(delta,Ii[i]-1);
    }
    double pp=delta*phi_delta*R*(t+T0)/(v*1000.0);
    itera=0;
    if(abs(p-pp)>err){
	    double v0=v;
	    double p0=pp;
	    double v1=v0+0.00001;
	    delta=1.0/(rhoc*v1);
	    phi_delta=ni[39]/delta;
	    for(int i=0;i<39;i++)
		phi_delta+=lefts[i]*pow(delta,Ii[i]-1);
	    double p1=delta*phi_delta*R*(t+T0)/(v1*1000.0);
	    v=v1+(p-p1)/(p1-p0)*(v1-v0);
	    delta=1.0/(rhoc*v);
	    phi_delta=ni[39]/delta;
	    for(int i=0;i<39;i++)
		phi_delta+=lefts[i]*pow(delta,Ii[i]-1);
	    pp=delta*phi_delta*R*(t+T0)/(v*1000.0);
	    while(abs(p-pp)>err){
		    itera++;
		    v0=v1;
		    p0=p1;
		    v1=v;
		    p1=pp;  
		    v=v1+(p-p1)/(p1-p0)*(v1-v0);
            delta=1.0/(rhoc*v);
            phi_delta=ni[39]/delta;
            for(int i=0;i<39;i++)
            phi_delta+=lefts[i]*pow(delta,Ii[i]-1);
            pp=delta*phi_delta*R*(t+T0)/(v*1000.0);    
	    }
    }
    return v;    
}
double IF97Region3::PT2V(double p,double t){
    double err=ERR;
    double v;
    if(p>40){ //&& p<=100
        if(t<=T3ab(p))
            v=PT2V3a(p,t);
        else
            v=PT2V3b(p,t);
    }else if(p>25){
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3ab(p))
            v=PT2V3d(p,t);
        else if(t<=T3ef(p))
            v=PT2V3e(p,t);
        else
            v=PT2V3f(p,t);
    }else if(p>23.5){
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3gh(p))
            v=PT2V3g(p,t);
        else if(t<=T3ef(p))
            v=PT2V3h(p,t);
        else if(t<=T3ij(p))
            v=PT2V3i(p,t);
        else if(t<=T3jk(p))
            v=PT2V3j(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>23){
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3gh(p))
            v=PT2V3l(p,t);
        else if(t<=T3ef(p))
            v=PT2V3h(p,t);
        else if(t<=T3ij(p))
            v=PT2V3i(p,t);
        else if(t<=T3jk(p))
            v=PT2V3j(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>22.5){
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3gh(p))
            v=PT2V3l(p,t);
        else if(t<=T3mn(p))
            v=PT2V3m(p,t);
        else if(t<=T3ef(p))
            v=PT2V3n(p,t);
        else if(t<=T3op(p))
            v=PT2V3o(p,t);
        else if(t<=T3ij(p))
            v=PT2V3p(p,t);
        else if(t<=T3jk(p))
            v=PT2V3j(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>22.11){
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3qu(p))
            v=PT2V3q(p,t);
        else if(t<=T3uv(p))
            v=PT2V3u(p,t);
        else if(t<=T3ef(p))
            v=PT2V3v(p,t);
        else if(t<=T3wx(p))
            v=PT2V3w(p,t);
        else if(t<=T3rx(p))
            v=PT2V3x(p,t);
        else if(t<=T3jk(p))
            v=PT2V3r(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>22.064){
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3qu(p))
            v=PT2V3q(p,t);
        else if(t<=T3uv(p))
            v=PT2V3u(p,t);
        else if(t<=T3ef(p))
            v=PT2V3y(p,t);
        else if(t<=T3wx(p))
            v=PT2V3z(p,t);
        else if(t<=T3rx(p))
            v=PT2V3x(p,t);
        else if(t<=T3jk(p))
            v=PT2V3r(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>21.93161551){//21.93161551==psat97(v=0.00264m3/kg)f
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3qu(p))
            v=PT2V3q(p,t);
        else if(t<=T3uv(p))
            v=PT2V3u(p,t);
        else if(t<=P2T(p))//<Tsat97(p)
            v=PT2V3y(p,t);         
        else if(t<=T3wx(p))
            v=PT2V3z(p,t);         
        else if(t<=T3rx(p))
            v=PT2V3x(p,t);            
        else if(t<=T3jk(p))
            v=PT2V3r(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>21.90096265){//21.90096265==psat97(v=0.00385m3/kg)f
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3qu(p))
            v=PT2V3q(p,t);
        else if(t<=P2T(p))//<Tsat97(p)
            v=PT2V3u(p,t);       
        else if(t<=T3wx(p))
            v=PT2V3z(p,t);         
        else if(t<=T3rx(p))
            v=PT2V3x(p,t);            
        else if(t<=T3jk(p))
            v=PT2V3r(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>21.04336731897525){//p>psat97(643.15K)
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=T3qu(p))
            v=PT2V3q(p,t);
        else if(t<=P2T(p))//<Tsat97(p)
            v=PT2V3u(p,t);         
        else if(t<=T3rx(p))
            v=PT2V3x(p,t);            
        else if(t<=T3jk(p))
            v=PT2V3r(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>20.5){
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=P2T(p))//t<=Tsat97(p)
            v=PT2V3s(p,t);
        else if(t<=T3jk(p))
            v=PT2V3r(p,t);
        else
            v=PT2V3k(p,t);
    }else if(p>19.00881189173929){//p>p3cd
        if(t<=T3cd(p))
            v=PT2V3c(p,t);
        else if(t<=P2T(p))//t<=Tsat97(p)
            v=PT2V3s(p,t);
        else
            v=PT2V3t(p,t);
    }else{//p>Psat97(623.15K)
        if(t<=P2T(p))//t<=Tsat97(p)
            v=PT2V3c(p,t);
        else
            v=PT2V3t(p,t);
    }
    double tau=Tc/(t+T0);
    double lefts[39];
    for(int i=0;i<39;i++){
        lefts[i]=ni[i]*Ii[i]*pow(tau,Ji[i]);
    }
    double delta=1.0/(rhoc*v);
    double phi_delta=ni[39]/delta;
    for(int i=0;i<39;i++){
        phi_delta+=lefts[i]*pow(delta,Ii[i]-1);
    }
    double pp=delta*phi_delta*R*(t+T0)/(v*1000.0);
    if(abs(p-pp)>err){
	    double v0=v;
	    double p0=pp;
	    double v1=v0+0.00001;
	    delta=1.0/(rhoc*v1);
	    phi_delta=ni[39]/delta;
	    for(int i=0;i<39;i++)
		phi_delta+=lefts[i]*pow(delta,Ii[i]-1);
	    double p1=delta*phi_delta*R*(t+T0)/(v1*1000.0);
	    v=v1+(p-p1)/(p1-p0)*(v1-v0);
	    delta=1.0/(rhoc*v);
	    phi_delta=ni[39]/delta;
	    for(int i=0;i<39;i++)
		phi_delta+=lefts[i]*pow(delta,Ii[i]-1);
	    pp=delta*phi_delta*R*(t+T0)/(v*1000.0);
	    while(abs(p-pp)>err){
		    v0=v1;
		    p0=p1;
		    v1=v;
		    p1=pp;  
		    v=v1+(p-p1)/(p1-p0)*(v1-v0);
            delta=1.0/(rhoc*v);
            phi_delta=ni[39]/delta;
            for(int i=0;i<39;i++)
            phi_delta+=lefts[i]*pow(delta,Ii[i]-1);
            pp=delta*phi_delta*R*(t+T0)/(v*1000.0);    
	    }
    }
    return v;    
}
/*P,H 2*/
void IF97Region3::PH2TV(double p,double h,double& t,double& v,int& itera){
    double err=ERR;
    if(h<=H3ab(p)){
        t=PH2T3a(p,h);
        v=PH2V3a(p,h);
    }else{
        t=PH2T3b(p,h);
        v=PH2V3b(p,h);
    }
    double hh=TV2H(t,v);
    itera=0;
    //if(abs(hh-h)>err){ 不要判断，强制进行迭代，因为PH2T/V3a/b得到的t和v算出的hh有时候刚好在err范围内等于h，但不符合的情况
        double t0=t;
        double h0=hh;
        double t1=t0-0.02;
        t1=t1<350?t0+0.02:t1;
        v=PT2V(p,t1);
        double h1=TV2H(t1,v);
        t=t1+(h-h1)/(h0-h1)*(t0-t1);
        v=PT2V(p,t);
        hh=TV2H(t,v);
        while(abs(hh-h)>err){
            itera++;
            t0=t1;
            h0=h1;
            t1=t;
            h1=hh;
            t=t1+(h-h1)/(h0-h1)*(t0-t1);
            v=PT2V(p,t);
            hh=TV2H(t,v);
        }
    //}
}
void IF97Region3::PH2TV(double p,double h,double& t,double& v){
    double err=ERR;
    if(h<=H3ab(p)){
        t=PH2T3a(p,h);
        v=PH2V3a(p,h);
    }else{
        t=PH2T3b(p,h);
        v=PH2V3b(p,h);
    }
    double hh=TV2H(t,v);
    //if(abs(hh-h)>err){ 不要判断，强制进行迭代，因为PH2T/V3a/b得到的t和v算出的hh有时候刚好在err范围内等于h，但不符合的情况
        double t0=t;
        double h0=hh;
        double t1=t0-0.02;
        t1=t1<350?t0+0.02:t1;
        v=PT2V(p,t1);
        double h1=TV2H(t1,v);
        t=t1+(h-h1)/(h0-h1)*(t0-t1);
        v=PT2V(p,t);
        hh=TV2H(t,v);
        while(abs(hh-h)>err){
            t0=t1;
            h0=h1;
            t1=t;
            h1=hh;
            t=t1+(h-h1)/(h0-h1)*(t0-t1);
            v=PT2V(p,t);
            hh=TV2H(t,v);
        }
    //}
}
/*P，S 2*/
void IF97Region3::PS2TV(double p,double s,double& t,double& v,int& itera){
    double err=ERR;
    if(s<=S3ab){
        t=PS2T3a(p,s);
        v=PS2V3a(p,s);
    }else{
        t=PS2T3b(p,s);
        v=PS2V3b(p,s);
    }
    double ss=TV2S(t,v);
    itera=0;
    //if(abs(ss-s)>err){ 不要判断，强制进行迭代，因为PS2T/V3a/b得到的t和v算出的ss有时候刚好在err范围内等于s，但不符合的情况
        double t0=t;
        double s0=ss;
        double t1=t0-0.2;
        t1=t1<350?t0+0.2:t1;
        v=PT2V(p,t1);
        double s1=TV2S(t1,v);
        t=t1+(s-s1)/(s0-s1)*(t0-t1);
        v=PT2V(p,t);
        ss=TV2S(t,v);
        while(abs(ss-s)>err){
            itera++;
            t0=t1;
            s0=s1;
            t1=t;
            s1=ss;
            t=t1+(s-s1)/(s0-s1)*(t0-t1);
            v=PT2V(p,t);
            ss=TV2S(t,v);
        }
    //}
}
void IF97Region3::PS2TV(double p,double s,double& t,double& v){
    double err=ERR;
    if(s<=S3ab){
        t=PS2T3a(p,s);
        v=PS2V3a(p,s);
    }else{
        t=PS2T3b(p,s);
        v=PS2V3b(p,s);
    }
    double ss=TV2S(t,v);
    //if(abs(ss-s)>err){ 不要判断，强制进行迭代，因为PS2T/V3a/b得到的t和v算出的ss有时候刚好在err范围内等于s，但不符合的情况
        double t0=t;
        double s0=ss;
        double t1=t0-0.2;
        t1=t1<350?t0+0.2:t1;
        v=PT2V(p,t1);
        double s1=TV2S(t1,v);
        t=t1+(s-s1)/(s0-s1)*(t0-t1);
        v=PT2V(p,t);
        ss=TV2S(t,v);
        while(abs(ss-s)>err){
            t0=t1;
            s0=s1;
            t1=t;
            s1=ss;
            t=t1+(s-s1)/(s0-s1)*(t0-t1);
            v=PT2V(p,t);
            ss=TV2S(t,v);
        }
    //}
}
/*H,S 2*/
void IF97Region3::HS2TVP(double h,double s,double& t,double& v, double& p,int& itera){ 
    double err=ERR;
    if(s<=S3ab){
        p=HS2P3a(h,s);
    }else{
        p=HS2P3b(h,s);
    }
    PS2TV(p,s,t,v);
    double hh=TV2H(t,v);
    itera=0;
    if(abs(hh-h)>err){
        double t0=t;
        double h0=hh;
        double t1=(t-0.1)<=350?t+0.1:t-0.1;
        double ss=TV2S(t1,v);
        if(abs(ss-s)>err){
            double v0=v;
            double s0=ss;
            double v1=v+0.0001;
            double s1=TV2S(t1,v1);
            v=v1+(s-s1)/(s0-s1)*(v0-v1);
            ss=TV2S(t1,v);
            while(abs(ss-s)>err && itera<40){
                itera++;
                v0=v1;
                s0=s1;
                v1=v;
                s1=ss;
                v=v1+(s-s1)/(s0-s1)*(v0-v1);
                ss=TV2S(t1,v);
            }
        }
        double h1=TV2H(t1,v);
        t=t1+(h-h1)/(h0-h1)*(t0-t1);
        ss=TV2S(t,v);
        if(abs(ss-s)>err){
            double v0=v;
            double s0=ss;
            double v1=v+0.0001;
            double s1=TV2S(t,v1);
            v=v1+(s-s1)/(s0-s1)*(v0-v1);
            ss=TV2S(t,v);
            while(abs(ss-s)>err  && itera<40){
                itera++;
                v0=v1;
                s0=s1;
                v1=v;
                s1=ss;
                v=v1+(s-s1)/(s0-s1)*(v0-v1);
                ss=TV2S(t,v);
            }
        }
        hh=TV2H(t,v);
        while(abs(hh-h)>err && itera<40){//p=97.506048,t=374.1对应的h，s输入时，不能收敛需要优化
            itera++;
            t0=t1;
            h0=h1;
            t1=t;
            h1=hh;
            t=t1+(h-h1)/(h0-h1)*(t0-t1);
            ss=TV2S(t,v);
            if(abs(ss-s)>err){
                itera++;
                double v0=v;
                double s0=ss;
                double v1=v+0.0001;
                double s1=TV2S(t,v1);
                v=v1+(s-s1)/(s0-s1)*(v0-v1);
                ss=TV2S(t,v);
                while(abs(ss-s)>err  && itera<40){
                    v0=v1;
                    s0=s1;
                    v1=v;
                    s1=ss;
                    v=v1+(s-s1)/(s0-s1)*(v0-v1);
                    ss=TV2S(t,v);
                }
            }
            hh=TV2H(t,v);
        }
    }
    p=TV2P(t,v);
}
void IF97Region3::HS2TVP(double h,double s,double& t,double& v, double& p){ 
    int itera; 
    double err=ERR;
    if(s<=S3ab){
        p=HS2P3a(h,s);
    }else{
        p=HS2P3b(h,s);
    }
    PS2TV(p,s,t,v);
    double hh=TV2H(t,v);
    itera=0;
    if(abs(hh-h)>err){
        double p0=p;
        double h0=hh;
        double p1=(p+0.1)>100?p-0.1:p+0.1;
        PS2TV(p1,s,t,v);
        double h1=TV2H(t,v);
        p=p1+(h-h1)/(h0-h1)*(p0-p1);
        PS2TV(p,s,t,v);
        hh=TV2H(t,v);
        while(abs(hh-h)>err && itera<10){//p=97.506048,t=374.1对应的h，s输入时，不能收敛需要优化
            itera++;
            p0=p1;
            h0=h1;
            p1=p;
            h1=hh;
            p=p1+(h-h1)/(h0-h1)*(p0-p1);
            PS2TV(p,s,t,v);
            hh=TV2H(t,v);
        }
    }
    //p=TV2P(t,v);
}
/*H,S 2* 辅助函数*/
double IF97Region3::HS2P3a(double h,double s){
    double ni[33]={
      0.770889828326934E1,
     -0.260835009128688E2,
      0.267416218930389E3,
      0.172221089496844E2,
     -0.293542332145970E3,
      0.614135601882478E3,
     -0.610562757725674E5,
     -0.651272251118219E8,
      0.735919313521937E5,
     -0.116646505914191E11,
      0.355267086434461E2,
     -0.596144543825955E3,
     -0.475842430145708E3,
      0.696781965359503E2,
      0.335674250377312E3,
      0.250526809130882E5,
      0.146997380630766E6,
      0.538069315091534E20,
      0.143619827291346E22,
      0.364985866165994E20,
     -0.254741561156775E4,
      0.240120197096563E28,
     -0.393847464679496E30,
      0.147073407024852E25,
     -0.426391250432059E32,
      0.194509340621077E39,
      0.666212132114896E24,
      0.706777016552858E34,
      0.175563621975576E42,
      0.108408607429124E29,
      0.730872705175151E44,
      0.159145847398870E25,
      0.377121605943324E41
    };
    int Ii[33]={
        0,0,0,1,1,1,1,1,2,2,3,3,3,4,4,4,4,5,6,
        7,8,10,10,14,18,20,22,22,24,28,28,32,32
    };
    int Ji[33]={
        0,1,5,0,3,4,8,14,6,16,0,2,3,0,1,4,5,28,28,
        24,1,32,36,22,28,36,16,28,36,16,36,10,28
    };
    double eta=h/2300.0;
    double sigma=s/4.4;
    double pi=0;
    for(int i=0;i<33;i++){
        pi+=ni[i]*pow(eta-1.01,Ii[i])*pow(sigma-0.750,Ji[i]);
    }
    return 99.0*pi;
}
double IF97Region3::HS2P3b(double h,double s){
    double ni[35]={
      0.125244360717979E-12,
     -0.126599322553713E-1,
      0.506878030140626E1,
      0.317847171154202E2,
     -0.391041161399932E6,
     -0.975733406392044E-10,
     -0.186312419488279E2,
      0.510973543414101E3,
      0.373847005822362E6,
      0.299804024666572E-7,
      0.200544393820342E2,
     -0.498030487662829E-5,
     -0.102301806360030E2,
      0.552819126990325E2,
     -0.206211367510878E3,
     -0.794012232324823E4,
      0.782248472028153E1,
     -0.586544326902468E2,
      0.355073647696481E4,
     -0.115303107290162E-3,
     -0.175092403171802E1,
      0.257981687748160E3,
     -0.727048374179467E3,
      0.121644822609198E-3,
      0.393137871762692E-1,
      0.704181005909296E-2,
     -0.829108200698110E2,
     -0.265178818131250,
      0.137531682453991E2,
     -0.522394090753046E2,
      0.240556298941048E4,
     -0.227361631268929E5,
      0.890746343932567E5,
     -0.239234565822486E8,
      0.568795808129714E10
    };
    int Ii[35]={
        -12,-12,-12,-12,-12,-10,-10,-10,-10,-8,-8,-6,-6,-6,-6,
        -5,-4,-4,-4,-3,-3,-3,-3,-2,-2,-1,0,2,2,5,6,8,10,14,14
    };
    int Ji[35]={
        2,10,12,14,20,2,10,14,18,2,8,2,6,7,8,
        10,4,5,8,1,3,5,6,0,1,0,3,0,1,0,1,1,1,3,7
    };
    double eta=h/2800.0;
    double sigma=s/5.3;
    double reciproca_pi=0;
    for(int i=0;i<35;i++){
        reciproca_pi+=ni[i]*pow(eta-0.681,Ii[i])*pow(sigma-0.792,Ji[i]);
    }
    return 16.6/reciproca_pi;
}
/*
tab4
    double ni[35]={
      0.125244360717979E-12,
     -0.126599322553713E-1,
      0.506878030140626E1,
      0.317847171154202E2,
     -0.391041161399932E6,
     -0.975733406392044E-10,
     -0.186312419488279E2,
      0.510973543414101E3,
      0.373847005822362E6,
      0.299804024666572E-7,
      0.200544393820342E2,
     -0.498030487662829E-5,
     -0.102301806360030E2,
      0.552819126990325E2,
     -0.206211367510878E3,
     -0.794012232324823E4,
      0.782248472028153E1,
     -0.586544326902468E2,
      0.355073647696481E4,
     -0.115303107290162E-3,
     -0.175092403171802E1,
      0.257981687748160E3,
     -0.727048374179467E3,
      0.121644822609198E-3,
      0.393137871762692E-1,
      0.704181005909296E-2,
     -0.829108200698110E2,
     -0.265178818131250,
      0.137531682453991E2,
     -0.522394090753046E2,
      0.240556298941048E4,
     -0.227361631268929E5,
      0.890746343932567E5,
     -0.239234565822486E8,
      0.568795808129714E10
    };
    int Ii[35]={
        -12,-12,-12,-12,-12,-10,-10,-10,-10,-8,-8,-6,-6,-6,-6,
        -5,-4,-4,-4,-3,-3,-3,-3,-2,-2,-1,0,2,2,5,6,8,10,14,14
    };
    int Ji[35]={
        2,10,12,14,20,2,10,14,18,2,8,2,6,7,8,
        10,4,5,8,1,3,5,6,0,1,0,3,0,1,0,1,1,1,3,7
    };
tab9
    double ni[27]={
      0.332171191705237,
      0.611217706323496E-3,
     -0.882092478906822E1,
     -0.455628192543250,
     -0.263483840850452E-4,
     -0.223949661148062E2,
     -0.428398660164013E1,
     -0.616679338856916,
     -0.146823031104040E2,
      0.284523138727299E3,
     -0.113398503195444E3,
      0.115671380760859E4,
      0.395551267359325E3,
     -0.154891257229285E1,
      0.194486637751291E2,
     -0.357915139457043E1,
     -0.335369414148819E1,
     -0.664426796332460,
      0.323321885383934E5,
      0.331766744667084E4,
     -0.223501257931087E5,
      0.573953875852936E7,
      0.173226193407919E3,
     -0.363968822121321E-1,
      0.834596332878346E-6,
      0.503611916682674E1,
      0.655444787064505E2
    };
    int Ii[27]={
        0,0,1,1,2,2,3,3,4,4,4,5,5,7,8,12,
        12,14,14,16,20,20,22,24,28,32,32,
    };
    int Ji[27]={
        14,36,3,16,0,5,4,36,4,16,24,18,24,
        1,4,2,4,1,22,10,12,28,8,3,0,6,8
    };
tab10
    double ni[19]={
      0.822673364673336,
      0.181977213534479,
     -0.112000260313624E-1,
     -0.746778287048033E-3,
     -0.179046263257381,
      0.424220110836657E-1,
     -0.341355823438768,
     -0.209881740853565E1,
     -0.822477343323596E1,
     -0.499684082076008E1,
      0.191413958471069,
      0.581062241093136E-1,
     -0.165505498701029E4,
      0.158870443421201E4,
     -0.850623535172818E2,
     -0.317714386511207E5,
     -0.945890406632871E5,
     -0.139273847088690E-5,
      0.631052532240980
    };
    int Ii[19]={
        0,0,0,0,2,3,4,4,5,5,6,7,7,7,10,10,10,32,32
    };
    int Ji[19]={
        1,4,10,16,1,36,3,16,20,36,4,2,28,32,14,32,36,0,6
    };
tab16
    double ni[30]={
     -0.524581170928788E3,
     -0.926947218142218E7,
     -0.237385107491666E3,
      0.210770155812776E11,
     -0.239494562010986E2,
      0.221802480294197E3,
     -0.510472533393438E7,
      0.124981396109147E7,
      0.200008436996201E10,
     -0.815158509791035E3,
     -0.157612685637523E3,
     -0.114200422332791E11,
      0.662364680776872E16,
     -0.227622818296144E19,
     -0.171048081348406E32,
      0.660788766938091E16,
      0.166320055886021E23,
     -0.218003784381501E30,
     -0.787276140295618E30,
      0.151062329700346E32,
      0.795732170300541E7,
      0.131957647355347E16,
     -0.325097068299140E24,
     -0.418600611419248E26,
      0.297478906557467E35,
     -0.953588761745473E20,
      0.166957699620939E25,
     -0.175407764869978E33,
      0.347581490626396E35,
     -0.710971318427851E39
    };
    int Ii[30]={
        1,1,2,2,4,4,7,8,8,10,12,12,18,20,24,28,
        28,28,28,28,32,32,32,32,32,36,36,36,36,36
    };
    int Ji[30]={
        8,24,4,32,1,2,7,5,12,1,0,7,10,12,32,8,
        12,20,22,24,2,7,12,14,24,10,12,20,22,28
    };
tab17
    double ni[16]={
      0.104351280732769E1,
     -0.227807912708513E1,
      0.180535256723202E1,
      0.420440834792042,
     -0.105721244834660E6,
      0.436911607493884E25,
     -0.328032702839753E12,
     -0.678686760804270E16,
      0.743957464645363E4,
     -0.356896445355761E20,
      0.167590585186801E32,
     -0.355028625419105E38,
      0.396611982166538E12,
     -0.414716268484468E41,
      0.359080103867382E19,
     -0.116994334851995E41
    };
    int Ii[16]={
        0,0,0,1,1,5,6,7,8,8,12,16,22,22,24,36
    };
    int Ji[16]={
        0,3,4,0,12,36,12,16,2,20,32,36,2,32,7,20
    };
tab23
    double ni[6]={
      0.913965547600543,
     -0.430944856041991E-4,
      0.603235694765419E2,
      0.117518273082168E-17,
      0.220000904781292,
     -0.690815545851641E2
    };
    int Ii[6]={
        0,1,1,3,5,6
    };
    int Ji[6]={
        0,-2,2,-12,-4,-3
    };
tab25
    double ni[25]={
      0.629096260829810E-3,
     -0.823453502583165E-3,
      0.515446951519474E-7,
     -0.117565945784945E1,
      0.348519684726192E1,
     -0.507837382408313E-11,
     -0.284637670005479E1,
     -0.236092263939673E1,
      0.601492324973779E1,
      0.148039650824546E1,
      0.360075182221907E-3,
     -0.126700045009952E-1,
     -0.122184332521413E7,
      0.149276502463272,
      0.698733471798484,
     -0.252207040114321E-1,
      0.147151930985213E-1,
     -0.108618917681849E1,
     -0.936875039816322E-3,
      0.819877897570217E2,
     -0.182041861521835E3,
      0.261907376402688E-5,
     -0.291626417025961E5,
      0.140660774926165E-4,
      0.783237062349385E7
    };
    int Ii[25]={
        -12,-10,-8,-4,-3,-2,-2,-2,-2,0,1,
        1,1,3,3,5,6,6,8,8,8,12,12,14,14
    };
    int Ji[25]={
        10,8,3,4,3,-6,2,3,4,0,-3,-2,10,-2,
        -1,-5,-6,-3,-8,-2,-1,-12,-1,-12,1
    };
tab28
    double ni[36]={
      0.179882673606601,
     -0.267507455199603,
      0.116276722612600E1,
      0.147545428713616,
     -0.512871635973248,
      0.421333567697984,
      0.563749522189870,
      0.429274443819153,
     -0.335704552142140E1,
      0.108890916499278E2,
     -0.248483390456012,
      0.304153221906390,
     -0.494819763939905,
      0.107551674933261E1,
      0.733888415457688E-1,
      0.140170545411085E-1,
     -0.106110975998808,
      0.168324361811875E-1,
      0.125028363714877E1,
      0.101316840309509E4,
     -0.151791558000712E1,
      0.524277865990866E2,
      0.230495545563912E5,
      0.249459806365456E-1,
      0.210796467412137E7,
      0.366836848613065E9,
     -0.144814105365163E9,
     -0.179276373003590E-2,
      0.489955602100459E10,
      0.471262212070518E3,
     -0.829294390198652E11,
     -0.171545662263191E4,
      0.355777682973575E7,
      0.586062760258436E12,
     -0.129887635078195E8,
      0.317247449371057E11
    };
    int Ii[36]={
        0,0,0,1,1,1,1,2,2,2,3,3,3,3,4,4,5,5,5,5,6,
        6,6,8,10,10,12,14,14,16,16,18,18,18,20,28
    };
    int Ji[36]={
        0,3,12,0,1,2,5,0,5,8,0,2,3,4,0,1,1,2,4,16,6,
        8,22,1,20,36,24,1,28,12,32,14,22,36,24,36
    };

*/
/*P,H or P,S 2* 辅助函数*/
double IF97Region3::H3ab(double p){//s=sc=4.41202148223476
    double pi=p/1.0;
    double ni[4]={
        0.201464004206875E4,
        0.374696550136983E1,
       -0.219921901054187E-1,
        0.875131686009950E-4
    };
    double eta=ni[0]+ni[1]*pi+ni[2]*pi*pi+ni[3]*pi*pi*pi;
    return 1.0*eta;
}
double IF97Region3::PH2T3a(double p,double h){
    double ni[31]={
     -0.133645667811215E-6,
      0.455912656802978E-5,
     -0.146294640700979E-4,
      0.639341312970080E-2,
      0.372783927268847E3,
     -0.718654377460447E4,
      0.573494752103400E6,
     -0.267569329111439E7,
     -0.334066283302614E-4,
     -0.245479214069597E-1,
      0.478087847764996E2,
      0.764664131818904E-5,
      0.128350627676972E-2,
      0.171219081377331E-1,
     -0.851007304583213E1,
     -0.136513461629781E-1,
     -0.384460997596657E-5,
      0.337423807911655E-2,
     -0.551624873066791,
      0.729202277107470,
     -0.992522757376041E-2,
     -0.119308831407288,
      0.793929190615421,
      0.454270731799386,
      0.209998591259910,
     -0.642109823904738E-2,
     -0.235155868604540E-1,
      0.252233108341612E-2,
     -0.764885133368119E-2,
      0.136176427574291E-1,
     -0.133027883575669E-1
    };
    int Ii[31]={
        -12,-12,-12,-12,-12,-12,-12,-12,-10,-10,-10,-8,-8,
        -8,-8,-5,-3,-2,-2,-2,-1,-1,0,0,1,3,3,4,4,10,12
    };    
    int Ji[31]={
        0,1,2,6,14,16,20,22,1,5,12,0,2,4,
        10,2,0,1,3,4,0,2,0,1,1,0,1,0,3,4,5
    };
    double pi=p/100.0;
    double eta=h/2300.0;
    double theta=0;
    for(int i=0;i<31;i++){
        theta+=ni[i]*pow(pi+0.240,Ii[i])*pow(eta-0.615,Ji[i]);
    }
    return 760.0*theta-273.15;
}
double IF97Region3::PH2T3b(double p,double h){
    double ni[33]={
      0.323254573644920E-4,
     -0.127575556587181E-3,
     -0.475851877356068E-3,
      0.156183014181602E-2,
      0.105724860113781,
     -0.858514221132534E2,
      0.724140095480911E3,
      0.296475810273257E-2,
     -0.592721983365988E-2,
     -0.126305422818666E-1,
     -0.115716196364853,
      0.849000969739595E2,
     -0.108602260086615E-1,
      0.154304475328851E-1,
      0.750455441524466E-1,
      0.252520973612982E-1,
     -0.602507901232996E-1,
     -0.307622221350501E1,
     -0.574011959864879E-1,
      0.503471360939849E1,
     -0.925081888584834,
      0.391733882917546E1,
     -0.773146007130190E2,
      0.949308762098587E4,
     -0.141043719679409E7,
      0.849166230819026E7,
      0.861095729446704,
      0.323346442811720,
      0.873281936020439,
     -0.436653048526683,
      0.286596714529479,
     -0.131778331276228,
      0.676682064330275E-2
    };
    int Ii[33]={
        -12,-12,-10,-10,-10,-10,-10,-8,-8,-8,-8,-8,-6,-6,
        -6,-4,-4,-3,-2,-2,-1,-1,-1,-1,-1,-1,0,0,1,3,5,6,8
    };
    int Ji[33]={
        0,1,0,1,5,10,12,0,1,2,4,10,0,1,2,0,1,
        5,0,4,2,4,6,10,14,16,0,2,1,1,1,1,1
    };
    double pi=p/100.0;
    double eta=h/2800.0;
    double theta=0;
    for(int i=0;i<33;i++){
        theta+=ni[i]*pow(pi+0.298,Ii[i])*pow(eta-0.720,Ji[i]);
    }
    return 860.0*theta-273.15;
}
double IF97Region3::PH2V3a(double p,double h){
    double ni[32]={
      0.529944062966028E-2,
     -0.170099690234461,
      0.111323814312927E2,
     -0.217898123145125E4,
     -0.506061827980875E-3,
      0.556495239685324,
     -0.943672726094016E1,
     -0.297856807561527,
      0.939353943717186E2,
      0.192944939465981E-1,
      0.421740664704763,
     -0.368914126282330E7,
     -0.737566847600639E-2,
     -0.354753242424366,
     -0.199768169338727E1,
      0.115456297059049E1,
      0.568366875815960E4,
      0.808169540124668E-2,
      0.172416341519307,
      0.104270175292927E1,
     -0.297691372792847,
      0.560394465163593,
      0.275234661176914,
     -0.148347894866012,
     -0.651142513478515E-1,
     -0.292468715386302E1,
      0.664876096952665E-1,
      0.352335014263844E1,
     -0.146340792313332E-1,
     -0.224503486668184E1,
      0.110533464706142E1,
     -0.408757344495612E-1
    };
    int Ii[32]={
        -12,-12,-12,-12,-10,-10,-10,-8,-8,-6,-6,-6,-4,
        -4,-3,-2,-2,-1,-1,-1,-1,0,0,1,1,1,2,2,3,4,5,8
    };
    int Ji[32]={
        6,8,12,18,4,7,10,5,12,3,4,22,2,3,7,
        3,16,0,1,2,3,0,1,0,1,2,0,2,0,2,2,2
    };
    double pi=p/100.0;
    double eta=h/2100.0;
    double omega=0;
    for(int i=0;i<32;i++){
        omega+=ni[i]*pow(pi+0.128,Ii[i])*pow(eta-0.727,Ji[i]);
    }
    return 0.0028*omega;
}
double IF97Region3::PH2V3b(double p,double h){

    double ni[30]={
     -0.225196934336318E-8,
      0.140674363313486E-7,
      0.233784085280560E-5,
     -0.331833715229001E-4,
      0.107956778514318E-2,
     -0.271382067378863,
      0.107202262490333E1,
     -0.853821329075382,
     -0.215214194340526E-4,
      0.769656088222730E-3,
     -0.431136580433864E-2,
      0.453342167309331,
     -0.507749535873652,
     -0.100475154528389E3,
     -0.219201924648793,
     -0.321087965668917E1,
      0.607567815637771E3,
      0.557686450685932E-3,
      0.187499040029550,
      0.905368030448107E-2,
      0.285417173048685,
      0.329924030996098E-1,
      0.239897419685483,
      0.482754995951394E1,
     -0.118035753702231E2,
      0.169490044091791,
     -0.179967222507787E-1,
      0.371810116332674E-1,
     -0.536288335065096E-1,
      0.160697101092520E1
    };
    int Ii[30]={
        -12,-12,-8,-8,-8,-8,-8,-8,-6,-6,-6,-6,-6,
        -6,-4,-4,-4,-3,-3,-2,-2,-1,-1,-1,-1,0,1,1,2,2
    };
    int Ji[30]={
        0,1,0,1,3,6,7,8,0,1,2,5,6,10,3,
        6,10,0,2,1,2,0,1,4,5,0,0,1,2,6
    };
    double pi=p/100.0;
    double eta=h/2800.0;
    double omega=0;
    for(int i=0;i<30;i++){
        omega+=ni[i]*pow(pi+0.0661,Ii[i])*pow(eta-0.720,Ji[i]);
    }
    return 0.0088*omega;    
}
double IF97Region3::PS2T3a(double p,double s){
    double ni[33]={
      0.150042008263875E10,
     -0.159397258480424E12,
      0.502181140217975E-3,
     -0.672057767855466E2,
      0.145058545404456E4,
     -0.823889534888890E4,
     -0.154852214233853,
      0.112305046746695E2,
     -0.297000213482822E2,
      0.438565132635495E11,
      0.137837838635464E-2,
     -0.297478527157462E1,
      0.971777947349413E13,
     -0.571527767052398E-4,
      0.288307949778420E5,
     -0.744428289262703E14,
      0.128017324848921E2,
     -0.368275545889071E3,
      0.664768904779177E16,
      0.449359251958880E-1,
     -0.422897836099655E1,
     -0.240614376434179,
     -0.474341365254924E1,
      0.724093999126110,
      0.923874349695897,
      0.399043655281015E1,
      0.384066651868009E-1,
     -0.359344365571848E-2,
     -0.735196448821653,
      0.188367048396131,
      0.141064266818704E-3,
     -0.257418501496337E-2,
      0.123220024851555E-2
    };
    int Ii[33]={
        -12,-12,-10,-10,-10,-10,-8,-8,-8,-8,-6,-6,-6,-5,
        -5,-5,-4,-4,-4,-2,-2,-1,-1,0,0,0,1,2,2,3,8,8,10
    };
    int Ji[33]={
        28,32,4,10,12,14,5,7,8,28,2,6,32,0,14,
        32,6,10,36,1,4,1,6,0,1,4,0,0,3,2,0,1,2
    };
    double pi=p/100.0;
    double sigma=s/4.4;
    double theta=0;
    for(int i=0;i<33;i++){
        theta+=ni[i]*pow(pi+0.240,Ii[i])*pow(sigma-0.703,Ji[i]);
    }
    return 760.0*theta-273.15;
}
double IF97Region3::PS2T3b(double p,double s){
    double ni[28]={
      0.527111701601660,
     -0.401317830052742E2,
      0.153020073134484E3,
     -0.224799398218827E4,
     -0.193993484669048,
     -0.140467557893768E1,
      0.426799878114024E2,
      0.752810643416743,
      0.226657238616417E2,
     -0.622873556909932E3,
     -0.660823667935396,
      0.841267087271658,
     -0.253717501764397E2,
      0.485708963532948E3,
      0.880531517490555E3,
      0.265015592794626E7,
     -0.359287150025783,
     -0.656991567673753E3,
      0.241768149185367E1,
      0.856873461222588,
      0.655143675313458,
     -0.213535213206406,
      0.562974957606348E-2,
     -0.316955725450471E15,
     -0.699997000152457E-3,
      0.119845803210767E-1,
      0.193848122022095E-4,
     -0.215095749182309E-4
    };
    int Ii[28]={
        -12,-12,-12,-12,-8,-8,-8,-6,-6,-6,-5,-5,
        -5,-5,-5,-4,-3,-3,-2,0,2,3,4,5,6,8,12,14
    };
    int Ji[28]={
        1,3,4,7,0,1,3,0,2,4,0,1,2,4,6,
        12,1,6,2,0,1,1,0,24,0,3,1,2
    };
    double pi=p/100.0;
    double sigma=s/5.3;
    double theta=0;
    for(int i=0;i<28;i++){
        theta+=ni[i]*pow(pi+0.760,Ii[i])*pow(sigma-0.818,Ji[i]);
    }
    return 860.0*theta-273.15;
}
double IF97Region3::PS2V3a(double p,double s){
    double ni[28]={
      0.795544074093975E2,
     -0.238261242984590E4,
      0.176813100617787E5,
     -0.110524727080379E-2,
     -0.153213833655326E2,
      0.297544599376982E3,
     -0.350315206871242E8,
      0.277513761062119,
     -0.523964271036888,
     -0.148011182995403E6,
      0.160014899374266E7,
      0.170802322663427E13,
      0.246866996006494E-3,
      0.165326084797980E1,
     -0.118008384666987,
      0.253798642355900E1,
      0.965127704669424,
     -0.282172420532826E2,
      0.203224612353823,
      0.110648186063513E1,
      0.526127948451280,
      0.277000018736321,
      0.108153340501132E1,
     -0.744127885357893E-1,
      0.164094443541384E-1,
     -0.680468275301065E-1,
      0.257988576101640E-1,
     -0.145749861944416E-3
    };
    int Ii[28]={
        -12,-12,-12,-10,-10,-10,-10,-8,-8,-8,-8,
        -6,-5,-4,-3,-3,-2,-2,-1,-1,0,0,0,1,2,4,5,6
    };
    int Ji[28]={
        10,12,14,4,8,10,20,5,6,14,16,28,
        1,5,2,4,3,8,1,2,0,1,3,0,0,2,2,0
    };
    double pi=p/100.0;
    double sigma=s/4.4;
    double omega=0;
    for(int i=0;i<28;i++){
        omega+=ni[i]*pow(pi+0.187,Ii[i])*pow(sigma-0.755,Ji[i]);
    }
    return 0.0028*omega;
}
double IF97Region3::PS2V3b(double p,double s){
    double ni[31]={
      0.591599780322238E-4,
     -0.185465997137856E-2,
      0.104190510480013E-1,
      0.598647302038590E-2,
     -0.771391189901699,
      0.172549765557036E1,
     -0.467076079846526E-3,
      0.134533823384439E-1,
     -0.808094336805495E-1,
      0.508139374365767,
      0.128584643361683E-2,
     -0.163899353915435E1,
      0.586938199318063E1,
     -0.292466667918613E1,
     -0.614076301499537E-2,
      0.576199014049172E1,
     -0.121613320606788E2,
      0.167637540957944E1,
     -0.744135838773463E1,
      0.378168091437659E-1,
      0.401432203027688E1,
      0.160279837479185E2,
      0.317848779347728E1,
     -0.358362310304853E1,
     -0.115995260446827E7,
      0.199256573577909,
     -0.122270624794624,
     -0.191449143716586E2,
     -0.150448002905284E-1,
      0.146407900162154E2,
     -0.327477787188230E1
    };
    int Ii[31]={
        -12,-12,-12,-12,-12,-12,-10,-10,-10,-10,-8,-5,-5,
        -5,-4,-4,-4,-4,-3,-2,-2,-2,-2,-2,-2,0,0,0,1,1,2
    };
    int Ji[31]={
        0,1,2,3,5,6,0,1,2,4,0,1,2,3,0,1,
        2,3,1,0,1,2,3,4,12,0,1,2,0,2,2        
    };
    double pi=p/100.0;
    double sigma=s/5.3;
    double omega=0;
    for(int i=0;i<31;i++){
        omega+=ni[i]*pow(pi+0.298,Ii[i])*pow(sigma-0.816,Ji[i]);
    }
    return 0.0088*omega;
}
double IF97Region3::H2P3sat(double h){
    double ni[14]={
      0.600073641753024,
     -0.936203654849857E1,
      0.246590798594147E2,
     -0.107014222858224E3,
     -0.915821315805768E14,
     -0.862332011700662E4,
     -0.235837344740032E2,
      0.252304969384128E18,
     -0.389718771997719E19,
     -0.333775713645296E23,
      0.356499469636328E11,
     -0.148547544720641E27,
      0.330611514838798E19,
      0.813641294467829E38
    };
    int Ii[14]={
        0,1,1,1,1,5,7,8,14,20,22,24,28,36
    };
    int Ji[14]={   
        0,1,3,4,36,3,0,24,16,16,3,18,8,24
    };
    double eta=h/2600.0;
    double p3sat=0;
    for(int i=0;i<14;i++){
        p3sat+=ni[i]*pow(eta-1.02,Ii[i])*pow(eta-0.608,Ji[i]);
    }
    return 22.0*p3sat;
}
double IF97Region3::S2P3sat(double s){
    double ni[10]={
      0.639767553612785,
     -0.129727445396014E2,
     -0.224595125848403E16,
      0.177466741801846E7,
      0.717079349571538E10,
     -0.378829107169011E18,
     -0.955586736431328E35,
      0.187269814676188E24,
      0.119254746466473E12,
      0.110649277244882E37
    };
    int Ii[10]={
        0,1,1,4,12,12,16,24,28,32
    };
    int Ji[10]={
        0,1,32,7,4,14,36,10,0,18
    };
    double sigma=s/5.2;
    double p3sat=0;
    for(int i=0;i<10;i++){
        p3sat+=ni[i]*pow(sigma-1.03,Ii[i])*pow(sigma-0.699,Ji[i]);
    }
    return 22.0*p3sat;
}
/*P,T 2 辅助函数*/
double IF97Region3::T3ab(double p){
    double ni[5]={
        0.154793642129415E4,
       -0.187661219490113E3,
        0.213144632222113E2,
       -0.191887498864292E4,
        0.918419702359447E3
    };
    int Ii[5]={
        0,1,2,-1,-2
    };
    double lnpi=log(p/1.0);
    double theta=0;
    for(int i=0;i<5;i++)
        theta+=ni[i]*pow(lnpi,Ii[i]);
    return 1.0*theta-T0;
}
double IF97Region3::T3cd(double p){
    double ni[4]={
        0.585276966696349E3,
        0.278233532206915E1,
       -0.127283549295878E-1,
        0.159090746562729E-3
    };
    int Ii[4]={
        0,1,2,3
    };
    double pi=p/1.0;
    double theta=0;
    for(int i=0;i<4;i++)
        theta+=ni[i]*pow(pi,Ii[i]);
    return 1.0*theta-T0;
}
double IF97Region3::T3ef(double p){
    double pi=p/1.0;
    double theta=3.727888004*(pi-22.064)+647.096;
    return 1.0*theta-T0;
}
double IF97Region3::T3gh(double p){
    double ni[5]={
       -0.249284240900418E5,
        0.428143584791546E4,
       -0.269029173140130E3,
        0.751608051114157E1,
       -0.787105249910383E-1
    };
    int Ii[5]={
        0,1,2,3,4
    };
    double pi=p/1.0;
    double theta=0;
    for(int i=0;i<5;i++)
        theta+=ni[i]*pow(pi,Ii[i]);
    return 1.0*theta-T0;
}
double IF97Region3::T3ij(double p){
    double ni[5]={
        0.584814781649163E3,
       -0.616179320924617E0,
        0.260763050899562E0,
       -0.587071076864459E-2,
        0.515308185433082E-4
    };
    int Ii[5]={
        0,1,2,3,4
    };
    double pi=p/1.0;
    double theta=0;
    for(int i=0;i<5;i++)
        theta+=ni[i]*pow(pi,Ii[i]);
    return 1.0*theta-T0;
}
double IF97Region3::T3jk(double p){
    double ni[5]={
        0.617229772068439E3,
       -0.770600270141675E1,
        0.697072596851896E0,
       -0.157391839848015E-1,
        0.137897492684194E-3      
    };
    int Ii[5]={
        0,1,2,3,4
    };
    double pi=p/1.0;
    double theta=0;
    for(int i=0;i<5;i++)
        theta+=ni[i]*pow(pi,Ii[i]);
    return 1.0*theta-T0;
}
double IF97Region3::T3mn(double p){
    double ni[4]={
        0.535339483742384E3,
        0.761978122720128E1,
       -0.158365725441648E0,
        0.192871054508108E-2
    };
    int Ii[4]={
        0,1,2,3
    };
    double pi=p/1.0;
    double theta=0;
    for(int i=0;i<4;i++)
        theta+=ni[i]*pow(pi,Ii[i]);
    return 1.0*theta-T0;
}
double IF97Region3::T3op(double p){
    double ni[5]={
        0.969461372400213E3,
       -0.332500170441278E3,
        0.642859598466067E2,
        0.773845935768222E3,
       -0.152313732937084E4
    };
    int Ii[5]={
        0,1,2,-1,-2
    };
    double lnpi=log(p/1.0);
    double theta=0;
    for(int i=0;i<5;i++)
        theta+=ni[i]*pow(lnpi,Ii[i]);
    return 1.0*theta-T0;
}
double IF97Region3::T3qu(double p){
    double ni[4]={
        0.565603648239126E3,
        0.529062258221222E1,
       -0.102020639611016E0,
        0.122240301070145E-2
    };
    int Ii[4]={
        0,1,2,3
    };
    double pi=p/1.0;
    double theta=0;
    for(int i=0;i<4;i++)
        theta+=ni[i]*pow(pi,Ii[i]);
    return 1.0*theta-T0;
}
double IF97Region3::T3rx(double p){
    double ni[4]={
        0.584561202520006E3,
       -0.102961025163669E1,
        0.243293362700452,
       -0.294905044740799E-2
    };
    int Ii[4]={
        0,1,2,3
    };
    double pi=p/1.0;
    double theta=0;
    for(int i=0;i<4;i++)
        theta+=ni[i]*pow(pi,Ii[i]);
    return 1.0*theta-T0;
}
double IF97Region3::T3uv(double p){
    double ni[4]={
        0.528199646263062E3,	
        0.890579602135307E1,  
       -0.222814134903755,
        0.286791682263697E-2
    };
    int Ii[4]={
        0,1,2,3
    };
    double pi=p/1.0;
    double theta=0;
    for(int i=0;i<4;i++)
        theta+=ni[i]*pow(pi,Ii[i]);
    return 1.0*theta-T0;
}
double IF97Region3::T3wx(double p){
    double ni[5]={
        0.728052609145380E1,
        0.973505869861952E2,
        0.147370491183191E2,
        0.329196213998375E3,
        0.873371668682417E3
    };
    int Ii[5]={
        0,1,2,-1,-2
    };
    double lnpi=log(p/1.0);
    double theta=0;
    for(int i=0;i<5;i++)
        theta+=ni[i]*pow(lnpi,Ii[i]);
    return 1.0*theta-T0;
}
double IF97Region3::PT2V3a(double p,double t){
    double pi=p/100.0;
    double theta=(t+T0)/760.0;
    double ni[30]={
        0.110879558823853E-2,
        0.572616740810616E3,
       -0.767051948380852E5,
       -0.253321069529674E-1,
        0.628008049345689E4,
        0.234105654131876E6,
        0.216867826045856E0,
       -0.156237904341963E3,
       -0.269893956176613E5,
       -0.180407100085505E-3,
        0.116732227668261E-2,
        0.266987040856040E2,
        0.282776617243286E5,
       -0.242431520029523E4,
        0.435217323022733E-3,
       -0.122494831387441E-1,
        0.179357604019989E1,
        0.442729521058314E2,
       -0.593223489018342E-2,
        0.453186261685774E0,
        0.135825703129140E1,
        0.408748415856745E-1,
        0.474686397863312E0,
        0.118646814997915E1,
        0.546987265727549E0,
        0.195266770452643E0,
       -0.502268790869663E-1,
       -0.369645308193377E0,
        0.633828037528420E-2,
        0.797441793901017E-1
    };
    int Ii[30]={
        -12,-12,-12,-10,-10,-10,-8,-8,-8,-6,-5,-5,-5,-4,-3,
        -3,-3,-3,-2,-2,-2,-1,-1,-1,0,0,1,1,2,2
    };
    int Ji[30]={
        5,10,12,5,10,12,5,8,10,1,1,5,10,8,0,
        1,3,6,0,2,3,0,1,2,0,1,0,2,0,2
    };
    double omega=0;
    for(int i=0;i<30;i++){
        omega+=ni[i]*pow(pi-0.085,Ii[i])*pow(theta-0.817,Ji[i]);
    }
    return 0.0024*omega;
}
double IF97Region3::PT2V3b(double p,double t){
    double pi=p/100.0;
    double theta=(t+T0)/860.0;
    double ni[32]={
       -0.827670470003621E-1,
        0.416887126010565E2,
        0.483651982197059E-1,
       -0.291032084950276E5,
       -0.111422582236948E3,
       -0.202300083904014E-1,
        0.294002509338515E3,
        0.140244997609658E3,
       -0.344384158811459E3,
        0.361182452612149E3,
       -0.140699677420738E4,
       -0.202023902676481E-2,
        0.171346792457471E3,
       -0.425597804058632E1,
        0.691346085000334E-5,
        0.151140509678925E-2,
       -0.416375290166236E-1,
       -0.413754957011042E2,
       -0.506673295721637E2,
       -0.572212965569023E-3,
        0.608817368401785E1,
        0.239600660256161E2,
        0.122261479925384E-1,
        0.216356057692938E1,
        0.398198903368642E0,
       -0.116892827834085E0,
       -0.102845919373532E0,
       -0.492676637589284E0,
        0.655540456406790E-1,
       -0.240462535078530E0,
       -0.269798180310075E-1,
        0.128369435967012E0
    };
    int Ii[32]={
        -12,-12,-10,-10,-8,-6,-6,-6,-5,-5,-5,-4,-4,-4,
        -3,-3,-3,-3,-3,-2,-2,-2,-1,-1,0,0,1,1,2,3,4,4
   };
    int Ji[32]={
        10,12,8,14,8,5,6,8,5,8,10,2,4,5,0,
        1,2,3,5,0,2,5,0,2,0,1,0,2,0,2,0,1 
    };
    double omega=0;
    for(int i=0;i<32;i++){
        omega+=ni[i]*pow(pi-0.280,Ii[i])*pow(theta-0.779,Ji[i]);
    }
    return 0.0041*omega;
}
double IF97Region3::PT2V3c(double p,double t){
    double pi=p/40.0;
    double theta=(t+T0)/690.0;
    double ni[35]={
      0.311967788763030E1,
      0.276713458847564E5,
      0.322583103403269E8,
     -0.342416065095363E3,
     -0.899732529907377E6,
     -0.793892049821251E8,
      0.953193003217388E2,
      0.229784742345072E4,
      0.175336675322499E6,
      0.791214365222792E7,
      0.319933345844209E-4,
     -0.659508863555767E2,
     -0.833426563212851E6,
      0.645734680583292E-1,
     -0.382031020570813E7,
      0.406398848470079E-4,
      0.310327498492008E2,
     -0.892996718483724E-3,
      0.234604891591616E3,
      0.377515668966951E4,
      0.158646812591361E-1,
      0.707906336241843,
      0.126016225146570E2,
      0.736143655772152,
      0.676544268999101,
     -0.178100588189137E2,
     -0.156531975531713,
      0.117707430048158E2,
      0.840143653860447E-1,
     -0.186442467471949,
     -0.440170203949645E2,
      0.123290423502494E7,
     -0.240650039730845E-1,
     -0.107077716660869E7,
      0.438319858566475E-1
    };
    int Ii[35]={
        -12,-12,-12,-10,-10,-10,-8,-8,-8,-6,-5,-5,-5,-4,-4,
        -3,-3,-2,-2,-2,-1,-1,-1,0,0,0,1,1,2,2,2,2,3,3,8
    };
    int Ji[35]={
        6,8,10,6,8,10,5,6,7,8,1,4,7,2,8,0,3,0,
        4,5,0,1,2,0,1,2,0,2,0,1,3,7,0,7,1
    };
    double omega=0;
    for(int i=0;i<35;i++){
        omega+=ni[i]*pow(pi-0.259,Ii[i])*pow(theta-0.903,Ji[i]);
    }
    return 0.0022*omega;
}
double IF97Region3::PT2V3d(double p,double t){
    double pi=p/40.0;
    double theta=(t+T0)/690.0;
    double ni[38]={
     -0.452484847171645E-9,
      0.315210389538801E-4,
     -0.214991352047545E-2,
      0.508058874808345E3,
     -0.127123036845932E8,
      0.115371133120497E13,
     -0.197805728776273E-15,
      0.241554806033972E-10,
     -0.156481703640525E-5,
      0.277211346836625E-2,
     -0.203578994462286E2,
      0.144369489909053E7,
     -0.411254217946539E11,
      0.623449786243773E-5,
     -0.221774281146038E2,
     -0.689315087933158E5,
     -0.195419525060713E8,
      0.316373510564015E4,
      0.224040754426988E7,
     -0.436701347922356E-5,
     -0.404213852833996E-3,
     -0.348153203414663E3,
     -0.385294213555289E6,
      0.135203700099403E-6,
      0.134648383271089E-3,
      0.125031835351736E6,
      0.968123678455841E-1,
      0.225660517512438E3,
     -0.190102435341872E-3,
     -0.299628410819229E-1,
      0.500833915372121E-2,
      0.387842482998411,
     -0.138535367777182E4,
      0.870745245971773,
      0.171946252068742E1,
     -0.326650121426383E-1,
      0.498044171727877E4,
      0.551478022765087E-2
    };
    int Ii[38]={
        -12,-12,-12,-12,-12,-12,-10,-10,-10,-10,-10,-10,-10,-8,-8,-8,
        -8,-6,-6,-5,-5,-5,-5,-4,-4,-4,-3,-3,-2,-2,-1,-1,-1,0,0,1,1,3
    };
    int Ji[38]={
        4,6,7,10,12,16,0,2,4,6,8,10,14,3,7,8,10,6,8,
        1,2,5,7,0,1,7,2,4,0,1,0,1,5,0,2,0,6,0
    };
    double omega=0;
    for(int i=0;i<38;i++){
        omega+=ni[i]*pow(pi-0.559,Ii[i])*pow(theta-0.939,Ji[i]);
    }
    return 0.0029*pow(omega,4);
}
double IF97Region3::PT2V3e(double p,double t){
    double pi=p/40.0;
    double theta=(t+T0)/710.0;
    double ni[29]={
      0.715815808404721E9,
     -0.114328360753449E12,
      0.376531002015720E-11,
     -0.903983668691157E-4,
      0.665695908836252E6,
      0.535364174960127E10,
      0.794977402335603E11,
      0.922230563421437E2,
     -0.142586073991215E6,
     -0.111796381424162E7,
      0.896121629640760E4,
     -0.669989239070491E4,
      0.451242538486834E-2,
     -0.339731325977713E2,
     -0.120523111552278E1,
      0.475992667717124E5,
     -0.266627750390341E6,
     -0.153314954386524E-3,
      0.305638404828265,
      0.123654999499486E3,
     -0.104390794213011E4,
     -0.157496516174308E-1,
      0.685331118940253,
      0.178373462873903E1,
     -0.544674124878910,
      0.204529931318843E4,
     -0.228342359328752E5,
      0.413197481515899,
     -0.341931835910405E2
    };
    int Ii[29]={
        -12,-12,-10,-10,-10,-10,-10,-8,-8,-8,-6,-5,-4,-4,
        -3,-3,-3,-2,-2,-2,-2,-1,0,0,1,1,1,2,2,
    };
    int Ji[29]={
        14,16,3,6,10,14,16,7,8,10,6,6,2,4,2,
        6,7,0,1,3,4,0,0,1,0,4,6,0,2
    };
    double omega=0;
    for(int i=0;i<29;i++){
        omega+=ni[i]*pow(pi-0.587,Ii[i])*pow(theta-0.918,Ji[i]);
    }
    return 0.0032*omega;
}
double IF97Region3::PT2V3f(double p,double t){
    double pi=p/40.0;
    double theta=(t+T0)/730.0;
    double ni[42]={
     -0.251756547792325E-7,
      0.601307193668763E-5,
     -0.100615977450049E-2,
      0.999969140252192,
      0.214107759236486E1,
     -0.165175571959086E2,
     -0.141987303638727E-2,
      0.269251915156554E1,
      0.349741815858722E2,
     -0.300208695771783E2,
     -0.131546288252539E1,
     -0.839091277286169E1,
      0.181545608337015E-9,
     -0.591099206478909E-3,
      0.152115067087106E1,
      0.252956470663225E-4,
      0.100726265203786E-14,
     -0.149774533860650E1,
     -0.793940970562969E-9,
     -0.150290891264717E-3,
      0.151205531275133E1,
      0.470942606221652E-5,
      0.195049710391712E-12,
     -0.911627886266077E-8,
      0.604374640201265E-3,
     -0.225132933900136E-15,
      0.610916973582981E-11,
     -0.303063908043404E-6,
     -0.137796070798409E-4,
     -0.919296736666106E-3,
      0.639288223132545E-9,
      0.753259479898699E-6,
     -0.400321478682929E-12,
      0.756140294351614E-8,
     -0.912082054034891E-11,
     -0.237612381140539E-7,
      0.269586010591874E-4,
     -0.732828135157839E-10,
      0.241995578306660E-9,
     -0.405735532730322E-3,
      0.189424143498011E-9,
     -0.486632965074563E-9
    };
    int Ii[42]={
        0,0,0,0,0,0,1,1,1,1,2,2,3,3,3,4,5,5,6,7,7,10,12,12,
        12,14,14,14,14,14,16,16,18,18,20,20,20,22,24,24,28,32
    };
    int Ji[42]={
        -3,-2,-1,0,1,2,-1,1,2,3,0,1,-5,-2,0,-3,-8,1,-6,-4,1,-6,-10,
        -8,-4,-12,-10,-8,-6,-4,-10,-8,-12,-10,-12,-10,-6,-12,-12,-4,-12,-12
    };
    double omega=0;
    for(int i=0;i<42;i++){
        float I=0.5*Ii[i];
        omega+=ni[i]*pow(pi-0.587,I)*pow(theta-0.891,Ji[i]);
    }
    return 0.0064*pow(omega,4);
}
double IF97Region3::PT2V3g(double p,double t){
    double pi=p/25.0;
    double theta=(t+T0)/660.0;
    double ni[38]={
      0.412209020652996E-4,
     -0.114987238280587E7,
      0.948180885032080E10,
     -0.195788865718971E18,
      0.496250704871300E25,
     -0.105549884548496E29,
     -0.758642165988278E12,
     -0.922172769596101E23,
      0.725379072059348E30,
     -0.617718249205859E2,
      0.107555033344858E5,
     -0.379545802336487E8,
      0.228646846221831E12,
     -0.499741093010619E7,
     -0.280214310054101E31,
      0.104915406769586E7,
      0.613754229168619E28,
      0.802056715528378E32,
     -0.298617819828065E8,
     -0.910782540134681E2,
      0.135033227281565E6,
     -0.712949383408211E19,
     -0.104578785289542E37,
      0.304331584444093E2,
      0.593250797959445E10,
     -0.364174062110798E28,
      0.921791403532461,
     -0.337693609657471,
     -0.724644143758508E2,
     -0.110480239272601,
      0.536516031875059E1,
     -0.291441872156205E4,
      0.616338176535305E40,
     -0.120889175861180E39,
      0.818396024524612E23,
      0.940781944835829E9,
     -0.367279669545448E5,
     -0.837513931798655E16
    };
    int Ii[38]={
        -12,-12,-12,-12,-12,-12,-10,-10,-10,-8,-8,-8,-8,-6,-6,
        -5,-5,-4,-3,-2,-2,-2,-2,-1,-1,-1,0,0,0,1,1,1,3,5,6,8,10,10
    };
    int Ji[38]={
        7,12,14,18,22,24,14,20,24,7,8,10,12,8,22,7,20,22,
        7,3,5,14,24,2,8,18,0,1,2,0,1,3,24,22,12,3,0,6
    };
    double omega=0;
    for(int i=0;i<38;i++){
        omega+=ni[i]*pow(pi-0.872,Ii[i])*pow(theta-0.971,Ji[i]);
    }
    return 0.0027*pow(omega,4);
}
double IF97Region3::PT2V3h(double p,double t){
    double pi=p/25.0;
    double theta=(t+T0)/660.0;
    double ni[29]={
      0.561379678887577E-1,
      0.774135421587083E10,
      0.111482975877938E-8,
     -0.143987128208183E-2,
      0.193696558764920E4,
     -0.605971823585005E9,
      0.171951568124337E14,
     -0.185461154985145E17,
      0.387851168078010E-16,
     -0.395464327846105E-13,
     -0.170875935679023E3,
     -0.212010620701220E4,
      0.177683337348191E8,
      0.110177443629575E2,
     -0.234396091693313E6,
     -0.656174421999594E7,
      0.156362212977396E-4,
     -0.212946257021400E1,
      0.135249306374858E2,
      0.177189164145813,
      0.139499167345464E4,
     -0.703670932036388E-2,
     -0.152011044389648,
      0.981916922991113E-4,
      0.147199658618076E-2,
      0.202618487925578E2,
      0.899345518944240,
     -0.211346402240858,
      0.249971752957491E2
    };
    int Ii[29]={
        -12,-12,-10,-10,-10,-10,-10,-10,-8,-8,-8,-8,-8,-6,-6,
        -6,-5,-5,-5,-4,-4,-3,-3,-2,-1,-1,0,1,1
    };
    int Ji[29]={
        8,12,4,6,8,10,14,16,0,1,6,7,8,4,
        6,8,2,3,4,2,4,1,2,0,0,2,0,0,2
    };
    double omega=0;
    for(int i=0;i<29;i++){
        omega+=ni[i]*pow(pi-0.898,Ii[i])*pow(theta-0.983,Ji[i]);
    }
    return 0.0032*pow(omega,4);
}
double IF97Region3::PT2V3i(double p,double t){
    double pi=p/25.0;
    double theta=(t+T0)/660.0;
    double ni[42]={
      0.106905684359136E1,
     -0.148620857922333E1,
      0.259862256980408E15,
     -0.446352055678749E-11,
     -0.566620757170032E-6,
     -0.235302885736849E-2,
     -0.269226321968839,
      0.922024992944392E1,
      0.357633505503772E-11,
     -0.173942565562222E2,
      0.700681785556229E-5,
     -0.267050351075768E-3,
     -0.231779669675624E1,
     -0.753533046979752E-12,
      0.481337131452891E1,
     -0.223286270422356E22,
     -0.118746004987383E-4,
      0.646412934136496E-2,
     -0.410588536330937E-9,
      0.422739537057241E20,
      0.313698180473812E-12,
      0.164395334345040E-23,
     -0.339823323754373E-5,
     -0.135268639905021E-1,
     -0.723252514211625E-14,
      0.184386437538366E-8,
     -0.463959533752385E-1,
     -0.992263100376750E14,
      0.688169154439335E-16,
     -0.222620998452197E-10,
     -0.540843018624083E-7,
      0.345570606200257E-2,
      0.422275800304086E11,
     -0.126974478770487E-14,
      0.927237985153679E-9,
      0.612670812016489E-13,
     -0.722693924063497E-11,
     -0.383669502636822E-3,
      0.374684572410204E-3,
     -0.931976897511086E5,
     -0.247690616026922E-1,
      0.658110546759474E2
    };
    int Ii[42]={
        0,0,0,1,1,1,1,2,3,3,4,4,4,5,5,5,7,7,8,8,10,12,12,12,
        14,14,14,14,18,18,18,18,18,20,20,22,24,24,32,32,36,36
    };
    int Ji[42]={
        0,1,10,-4,-2,-1,0,0,-5,0,-3,-2,-1,-6,-1,12,-4,-3,-6,10,-8,-12,-6,
        -4,-10,-8,-4,5,-12,-10,-8,-6,2,-12,-10,-12,-12,-8,-10,-5,-10,-8
    };
    double omega=0;
    for(int i=0;i<42;i++){
        float I=0.5*Ii[i];
        omega+=ni[i]*pow(pi-0.910,I)*pow(theta-0.984,Ji[i]);
    }
    return 0.0041*pow(omega,4);
}
double IF97Region3::PT2V3j(double p,double t){
    double pi=p/25.0;
    double theta=(t+T0)/670.0;
    double ni[29]={
     -0.111371317395540E-3,
      0.100342892423685E1,
      0.530615581928979E1,
      0.179058760078792E-5,
     -0.728541958464774E-3,
     -0.187576133371704E2,
      0.199060874071849E-2,
      0.243574755377290E2,
     -0.177040785499444E-3,
     -0.259680385227130E-2,
     -0.198704578406823E3,
      0.738627790224287E-4,
     -0.236264692844138E-2,
     -0.161023121314333E1,
      0.622322971786473E4,
     -0.960754116701669E-8,
     -0.510572269720488E-10,
      0.767373781404211E-2,
      0.663855469485254E-14,
     -0.717590735526745E-9,
      0.146564542926508E-4,
      0.309029474277013E-11,
     -0.464216300971708E-15,
     -0.390499637961161E-13,
     -0.236716126781431E-9,
      0.454652854268717E-11,
     -0.422271787482497E-2,
      0.283911742354706E-10,
      0.270929002720228E1
    };
    int Ii[29]={
        0,0,0,1,1,1,2,2,3,4,4,5,5,5,6,10,12
        ,12,14,14,14,16,18,20,20,24,24,28,28
    };
    int Ji[29]={
        -1,0,1,-2,-1,1,-1,1,-2,-2,2,-3,-2,0,3,-6,
        -8,-3,-10,-8,-5,-10,-12,-12,-10,-12,-6,-12,-5
    };
    double omega=0;
    for(int i=0;i<29;i++){
        float I=0.5*Ii[i];
        omega+=ni[i]*pow(pi-0.875,I)*pow(theta-0.964,Ji[i]);
    }
    return 0.0054*pow(omega,4);
}
double IF97Region3::PT2V3k(double p,double t){
    double pi=p/25.0;
    double theta=(t+T0)/680.0;
    double ni[34]={
     -0.401215699576099E9,
      0.484501478318406E11,
      0.394721471363678E-14,
      0.372629967374147E5,
     -0.369794374168666E-29,
     -0.380436407012452E-14,
      0.475361629970233E-6,
     -0.879148916140706E-3,
      0.844317863844331,
      0.122433162656600E2,
     -0.104529634830279E3,
      0.589702771277429E3,
     -0.291026851164444E14,
      0.170343072841850E-5,
     -0.277617606975748E-3,
     -0.344709605486686E1,
      0.221333862447095E2,
     -0.194646110037079E3,
      0.808354639772825E-15,
     -0.180845209145470E-10,
     -0.696664158132412E-5,
     -0.181057560300994E-2,
      0.255830298579027E1,
      0.328913873658481E4,
     -0.173270241249904E-18,
     -0.661876792558034E-6,
     -0.395688923421250E-2,
      0.604203299819132E-17,
     -0.400879935920517E-13,
      0.160751107464958E-8,
      0.383719409025556E-4,
     -0.649565446702457E-14,
     -0.149095328506000E-11,
      0.541449377329581E-8
    };
    int Ii[34]={
        -2,-2,-1,-1,0,0,0,0,0,0,0,0,0,1,1,1,1,
        1,2,2,2,2,2,2,5,5,5,6,6,6,6,8,10,12
    };
    int Ji[34]={
        10,12,-5,6,-12,-6,-2,-1,0,1,2,3,14,-3,-2,0,1,2,
        -8,-6,-3,-2,0,4,-12,-6,-3,-12,-10,-8,-5,-12,-12,-10
    };
    double omega=0;
    for(int i=0;i<34;i++){
        omega+=ni[i]*pow(pi-0.802,Ii[i])*pow(theta-0.935,Ji[i]);
    }
    return 0.0077*omega;
}
double IF97Region3::PT2V3l(double p,double t){
    double pi=p/24.0;
    double theta=(t+T0)/650.0;
    double ni[43]={
      0.260702058647537E10,
     -0.188277213604704E15,
      0.554923870289667E19,
     -0.758966946387758E23,
      0.413865186848908E27,
     -0.815038000738060E12,
     -0.381458260489955E33,
     -0.123239564600519E-1,
      0.226095631437174E8,
     -0.495017809506720E12,
      0.529482996422863E16,
     -0.444359478746295E23,
      0.521635864527315E35,
     -0.487095672740742E55,
     -0.714430209937547E6,
      0.127868634615495,
     -0.100752127917598E2,
      0.777451437960990E7,
     -0.108105480796471E25,
     -0.357578581169659E-5,
     -0.212857169423484E1,
      0.270706111085238E30,
     -0.695953622348829E33,
      0.110609027472280,
      0.721559163361354E2,
     -0.306367307532219E15,
      0.265839618885530E-4,
      0.253392392889754E-1,
     -0.214443041836579E3,
      0.937846601489667,
      0.223184043101700E1,
      0.338401222509191E2,
      0.494237237179718E21,
     -0.198068404154428,
     -0.141415349881140E31,
     -0.993862421613651E2,
      0.125070534142731E3,
     -0.996473529004439E3,
      0.473137909872765E5,
      0.116662121219322E33,
     -0.315874976271533E16,
     -0.445703369196945E33,
      0.642794932373694E33
    };
    int Ii[43]={
        -12,-12,-12,-12,-12,-10,-10,-8,-8,-8,-8,-8,-8,-8,-6,-5,-5,-4,-4,
        -3,-3,-3,-3,-2,-2,-2,-1,-1,-1,0,0,0,0,1,1,2,4,5,5,6,10,10,14,
    };
    int Ji[43]={
        14,16,18,20,22,14,24,6,10,12,14,18,24,36,8,4,5,7,16,1,3,
        18,20,2,3,10,0,1,3,0,1,2,12,0,16,1,0,0,1,14,4,12,10,
    };
    double omega=0;
    for(int i=0;i<43;i++){
        omega+=ni[i]*pow(pi-0.908,Ii[i])*pow(theta-0.989,Ji[i]);
    }
    return 0.0026*pow(omega,4);
}
double IF97Region3::PT2V3m(double p,double t){
    double pi=p/23.0;
    double theta=(t+T0)/650.0;
    double ni[40]={
        0.811384363481847E0,
       -0.568199310990094E4,
       -0.178657198172556E11,
        0.795537657613427E32,
       -0.814568209346872E5,
       -0.659774567602874E8,
       -0.152861148659302E11,
       -0.560165667510446E12,
        0.458384828593949E6,
       -0.385754000383848E14,
        0.453735800004273E8,
        0.939454935735563E12,
        0.266572856432938E28,
       -0.547578313899097E10,
        0.200725701112386E15,
        0.185007245563239E13,
        0.185135446828337E9,
       -0.170451090076385E12,
        0.157890366037614E15,
       -0.202530509748774E16,
        0.368193926183570E60,
        0.170215539458936E18,
        0.639234909918741E42,
       -0.821698160721956E15,
       -0.795260241872306E24,
        0.233415869478510E18,
       -0.600079934586803E23,
        0.594584382273384E25,
        0.189461279349492E40,
       -0.810093428842645E46,
        0.188813911076809E22,
        0.111052244098768E36,
        0.291133958602503E46,
       -0.329421923951460E22,
       -0.137570282536696E26,
        0.181508996303902E28,
       -0.346865122768353E30,
       -0.211961148774260E38,
       -0.128617899887675E49,
        0.479817895699239E65
    };
    int Ii[40]={
        0,3,8,20,1,3,4,5,1,6,2,4,14,2,5,3,0,1,1,1,28,
        2,16,0,5,0,3,4,12,16,1,8,14,0,2,3,4,8,14,24
   };
    int Ji[40]={
        0,0,0,2,5,5,5,5,6,6,7,8,8,10,10,12,14,14,18,20,20,22,
        22,24,24,28,28,28,28,28,32,32,32,36,36,36,36,36,36,36
     };
    double omega=0;
    for(int i=0;i<40;i++){
        float J=0.25*Ji[i];
        omega+=ni[i]*pow(pi-1.00,Ii[i])*pow(theta-0.997,J);
    }
    return 0.0028*omega;
}
double IF97Region3::PT2V3n(double p,double t){
    double pi=p/23.0;
    double theta=(t+T0)/650.0;
    double ni[39]={
      0.280967799943151E-38,
      0.614869006573609E-30,
      0.582238667048942E-27,
      0.390628369238462E-22,
      0.821445758255119E-20,
      0.402137961842776E-14,
      0.651718171878301E-12,
     -0.211773355803058E-7,
      0.264953354380072E-2,
     -0.135031446451331E-31,
     -0.607246643970893E-23,
     -0.402352115234494E-18,
     -0.744938506925544E-16,
      0.189917206526237E-12,
      0.364975183508473E-5,
      0.177274872361946E-25,
     -0.334952758812999E-18,
     -0.421537726098389E-8,
     -0.391048167929649E-1,
      0.541276911564176E-13,
      0.705412100773699E-11,
      0.258585887897486E-8,
     -0.493111362030162E-10,
     -0.158649699894543E-5,
     -0.525037427886100,
      0.220019901729615E-2,
     -0.643064132636925E-2,
      0.629154149015048E2,
      0.135147318617061E3,
      0.240560808321713E-6,
     -0.890763306701305E-3,
     -0.440209599407714E4,
     -0.302807107747776E3,
      0.159158748314599E4,
      0.232534272709876E6,
     -0.792681207132600E6,
     -0.869871364662769E11,
      0.354542769185671E12,
      0.400849240129329E15
    };
    int Ii[39]={
        0,3,4,6,7,10,12,14,18,0,3,5,6,8,12,0,3,7,
        12,2,3,4,2,4,7,4,3,5,6,0,0,3,1,0,1,0,1,0,1
  };
    int Ji[39]={
        -12,-12,-12,-12,-12,-12,-12,-12,-12,-10,-10,-10,-10,-10,-10,-8,
        -8,-8,-8,-6,-6,-6,-5,-5,-5,-4,-3,-3,-3,-2,-1,-1,0,1,1,2,4,5,6 
    };
    double omega=0;
    for(int i=0;i<39;i++){
        omega+=ni[i]*pow(pi-0.976,Ii[i])*pow(theta-0.997,Ji[i]);
    }
    return 0.0031*exp(omega);
}
double IF97Region3::PT2V3o(double p,double t){
    double pi=p/23.0;
    double theta=(t+T0)/650.0;
    double ni[24]={
      0.128746023979718E-34,
     -0.735234770382342E-11,
      0.289078692149150E-2,
      0.244482731907223,
      0.141733492030985E-23,
     -0.354533853059476E-28,
     -0.594539202901431E-17,
     -0.585188401782779E-8,
      0.201377325411803E-5,
      0.138647388209306E1,
     -0.173959365084772E-4,
      0.137680878349369E-2,
      0.814897605805513E-14,
      0.425596631351839E-25,
     -0.387449113787755E-17,
      0.139814747930240E-12,
     -0.171849638951521E-2,
      0.641890529513296E-21,
      0.118960578072018E-10,
     -0.155282762571611E-17,
      0.233907907347507E-7,
     -0.174093247766213E-12,
      0.377682649089149E-8,
     -0.516720236575302E-10
    };
    int Ii[24]={
        0,0,0,2,3,4,4,4,4,4,5,5,6,7,
        8,8,8,10,10,14,14,20,20,24
    };
    int Ji[24]={
        -12,-4,-1,-1,-10,-12,-8,-5,-4,-1,-4,-3,-8,
        -12,-10,-8,-4,-12,-8,-12,-8,-12,-10,-12
    };
    double omega=0;
    for(int i=0;i<24;i++){
        float I=0.5*Ii[i];
        omega+=ni[i]*pow(pi-0.974,I)*pow(theta-0.996,Ji[i]);
    }
    return 0.0034*omega;
}
double IF97Region3::PT2V3p(double p,double t){
    double pi=p/23.0;
    double theta=(t+T0)/650.0;
    double ni[27]={
     -0.982825342010366E-4,
      0.105145700850612E1,
      0.116033094095084E3,
      0.324664750281543E4,
     -0.123592348610137E4,
     -0.561403450013495E-1,
      0.856677401640869E-7,
      0.236313425393924E3,
      0.972503292350109E-2,
     -0.103001994531927E1,
     -0.149653706199162E-8,
     -0.215743778861592E-4,
     -0.834452198291445E1,
      0.586602660564988,
      0.343480022104968E-25,
      0.816256095947021E-5,
      0.294985697916798E-2,
      0.711730466276584E-16,
      0.400954763806941E-9,
      0.107766027032853E2,
     -0.409449599138182E-6,
     -0.729121307758902E-5,
      0.677107970938909E-8,
      0.602745973022975E-7,
     -0.382323011855257E-10,
      0.179946628317437E-2,
     -0.345042834640005E-3       
    };
    int Ii[27]={
        0,0,0,0,1,2,3,3,4,6,7,7,8,10,12,12,
        12,14,14,14,16,18,20,22,24,24,36
    };
    int Ji[27]={
        -1,0,1,2,1,-1,-3,0,-2,-2,-5,-4,-2,-3,-12,
        -6,-5,-10,-8,-3,-8,-8,-10,-10,-12,-8,-12
    };
    double omega=0;
    for(int i=0;i<27;i++){
        float I=0.5*Ii[i];
        omega+=ni[i]*pow(pi-0.972,I)*pow(theta-0.997,Ji[i]);
    }
    return 0.0041*omega;
}
double IF97Region3::PT2V3q(double p,double t){
    double pi=p/23.0;
    double theta=(t+T0)/650.0;
    double ni[24]={
     -0.820433843259950E5,
      0.473271518461586E11,
     -0.805950021005413E-1,
      0.328600025435980E2,
     -0.356617029982490E4,
     -0.172985781433335E10,
      0.351769232729192E8,
     -0.775489259985144E6,
      0.710346691966018E-4,
      0.993499883820274E5,
     -0.642094171904570,
     -0.612842816820083E4,
      0.232808472983776E3,
     -0.142808220416837E-4,
     -0.643596060678456E-2,
     -0.428577227475614E1,
      0.225689939161918E4,
      0.100355651721510E-2,
      0.333491455143516,
      0.109697576888873E1,
      0.961917379376452,
     -0.838165632204598E-1,
      0.247795908411492E1,
     -0.319114969006533E4
    };
    int Ii[24]={
        -12,-12,-10,-10,-10,-10,-8,-6,-5,-5,
        -4,-4,-3,-2,-2,-2,-2,-1,-1,-1,0,1,1,1
    };
    int Ji[24]={
        10,12,6,7,8,10,8,6,2,5,3,
        4,3,0,1,2,4,0,1,2,0,0,1,3
    };
    double omega=0;
    for(int i=0;i<24;i++){
        omega+=ni[i]*pow(pi-0.848,Ii[i])*pow(theta-0.983,Ji[i]);
    }
    return 0.0022*pow(omega,4);
}
double IF97Region3::PT2V3r(double p,double t){
    double pi=p/23.0;
    double theta=(t+T0)/650.0;
    double ni[27]={
      0.144165955660863E-2,
     -0.701438599628258E13,
     -0.830946716459219E-16,
      0.261975135368109,
      0.393097214706245E3,
     -0.104334030654021E5,
      0.490112654154211E9,
     -0.147104222772069E-3,
      0.103602748043408E1,
      0.305308890065089E1,
     -0.399745276971264E7,
      0.569233719593750E-11,
     -0.464923504407778E-1,
     -0.535400396512906E-17,
      0.399988795693162E-12,
     -0.536479560201811E-6,
      0.159536722411202E-1,
      0.270303248860217E-14,
      0.244247453858506E-7,
     -0.983430636716454E-5,
      0.663513144224454E-1,
     -0.993456957845006E1,
      0.546491323528491E3,
     -0.143365406393758E5,
      0.150764974125511E6,
     -0.337209709340105E-9,
      0.377501980025469E-8
    };
    int Ii[27]={
        -8,-8,-3,-3,-3,-3,-3,0,0,0,0,3,3,8,
        8,8,8,10,10,10,10,10,10,10,10,12,14
    };
    int Ji[27]={
        6,14,-3,3,4,5,8,-1,0,1,5,-6,-2,-12,-10,
        -8,-5,-12,-10,-8,-6,-5,-4,-3,-2,-12,-12
    };
    double omega=0;
    for(int i=0;i<27;i++){
        omega+=ni[i]*pow(pi-0.874,Ii[i])*pow(theta-0.982,Ji[i]);
    }
    return 0.0054*omega;
}
double IF97Region3::PT2V3s(double p,double t){
    double pi=p/21.0;
    double theta=(t+T0)/640.0;
    double ni[29]={
     -0.532466612140254E23,
      0.100415480000824E32,
     -0.191540001821367E30,
      0.105618377808847E17,
      0.202281884477061E59,
      0.884585472596134E8,
      0.166540181638363E23,
     -0.313563197669111E6,
     -0.185662327545324E54,
     -0.624942093918942E-1,
     -0.504160724132590E10,
      0.187514491833092E5,
      0.121399979993217E-2,
      0.188317043049455E1,
     -0.167073503962060E4,
      0.965961650599775,
      0.294885696802488E1,
     -0.653915627346115E5,
      0.604012200163444E50,
     -0.198339358557937,
     -0.175984090163501E58,
      0.356314881403987E1,
     -0.575991255144384E3,
      0.456213415338071E5,
     -0.109174044987829E8,
      0.437796099975134E34,
     -0.616552611135792E46,
      0.193568768917797E10,
      0.950898170425042E54
    };
    int Ii[29]={
        -12,-12,-10,-8,-6,-5,-5,-4,-4,-3,-3,-2,
        -1,-1,-1,0,0,0,0,1,1,3,3,3,4,4,4,5,14

    };
    int Ji[29]={
        20,24,22,14,36,8,16,6,32,3,8,4,1,2,
        3,0,1,4,28,0,32,0,1,2,3,18,24,4,24
    };
    double omega=0;
    for(int i=0;i<29;i++){
        omega+=ni[i]*pow(pi-0.886,Ii[i])*pow(theta-0.990,Ji[i]);
    }
    return 0.0022*pow(omega,4);
}
double IF97Region3::PT2V3t(double p,double t){
    double pi=p/20.0;
    double theta=(t+T0)/650.0;
    double ni[33]={
      0.155287249586268E1,
      0.664235115009031E1,
     -0.289366236727210E4,
     -0.385923202309848E13,
     -0.291002915783761E1,
     -0.829088246858083E12,
      0.176814899675218E1,
     -0.534686695713469E9,
      0.160464608687834E18,
      0.196435366560186E6,
      0.156637427541729E13,
     -0.178154560260006E1,
     -0.229746237623692E16,
      0.385659001648006E8,
      0.110554446790543E10,
     -0.677073830687349E14,
     -0.327910592086523E31,
     -0.341552040860644E51,
     -0.527251339709047E21,
      0.245375640937055E24,
     -0.168776617209269E27,
      0.358958955867578E29,
     -0.656475280339411E36,
      0.355286045512301E39,
      0.569021454413270E58,
     -0.700584546433113E48,
     -0.705772623326374E65,
      0.166861176200148E53,
     -0.300475129680486E61,
     -0.668481295196808E51,
      0.428432338620678E69,
     -0.444227367758304E72,
     -0.281396013562745E77
    };
    int Ii[33]={
        0,0,0,0,1,1,2,2,2,3,3,4,4,7,7,7,7,7,10,
        10,10,10,10,18,20,22,22,24,28,32,32,32,36
    };
    int Ji[33]={
        0,1,4,12,0,10,0,6,14,3,8,0,10,3,4,7,20,36,10,
        12,14,16,22,18,32,22,36,24,28,22,32,36,36
    };
    double omega=0;
    for(int i=0;i<33;i++){
        omega+=ni[i]*pow(pi-0.803,Ii[i])*pow(theta-1.02,Ji[i]);
    }
    return 0.0088*omega;
}
double IF97Region3::PT2V3u(double p,double t){
    double pi=p/23.0;
    double theta=(t+T0)/650.0;
    double ni[38]={
      0.122088349258355E18,
      0.104216468608488E10,
     -0.882666931564652E16,
      0.259929510849499E20,
      0.222612779142211E15,
     -0.878473585050085E18,
     -0.314432577551552E22,
     -0.216934916996285E13,
      0.159079648196849E21,
     -0.339567617303423E3,
      0.884387651337836E13,
     -0.843405926846418E21,
      0.114178193518022E2,
     -0.122708229235641E-3,
     -0.106201671767107E3,
      0.903443213959313E25,
     -0.693996270370852E28,
      0.648916718965575E-8,
      0.718957567127851E4,
      0.105581745346187E-2,
     -0.651903203602581E15,
     -0.160116813274676E25,
     -0.510254294237837E-8,
     -0.152355388953402,
      0.677143292290144E12,
      0.276378438378930E15,
      0.116862983141686E-1,
     -0.301426947980171E14,
      0.169719813884840E-7,
      0.104674840020929E27,
     -0.108016904560140E5,
     -0.990623601934295E-12,
      0.536116483602738E7,
      0.226145963747881E22,
     -0.488731565776210E-9,
      0.151001548880670E-4,
     -0.227700464643920E5,
     -0.781754507698846E28
    };
    int Ii[38]={
        -12,-10,-10,-10,-8,-8,-8,-6,-6,-5,-5,-5,-3,-1,-1,-1,
        -1,0,0,1,2,2,3,5,5,5,6,6,8,8,10,12,12,12,14,14,14,14
    };
    int Ji[38]={
        14,10,12,14,10,12,14,8,12,4,8,12,2,-1,1,12,14,-3,1,
        -2,5,10,-5,-4,2,3,-5,2,-8,8,-4,-12,-4,4,-12,-10,-6,6
    };
    double omega=0;
    for(int i=0;i<38;i++){
        omega+=ni[i]*pow(pi-0.902,Ii[i])*pow(theta-0.988,Ji[i]);
    }
    return 0.0026*omega;
}
double IF97Region3::PT2V3v(double p,double t){
    double pi=p/23.0;
    double theta=(t+T0)/650.0;
    double ni[39]={
     -0.415652812061591E-54,
      0.177441742924043E-60,
     -0.357078668203377E-54,
      0.359252213604114E-25,
     -0.259123736380269E2,
      0.594619766193460E5,
     -0.624184007103158E11,
      0.313080299915944E17,
      0.105006446192036E-8,
     -0.192824336984852E-5,
      0.654144373749937E6,
      0.513117462865044E13,
     -0.697595750347391E19,
     -0.103977184454767E29,
      0.119563135540666E-47,
     -0.436677034051655E-41,
      0.926990036530639E-29,
      0.587793105620748E21,
      0.280375725094731E-17,
     -0.192359972440634E23,
      0.742705723302738E27,
     -0.517429682450605E2,
      0.820612048645469E7,
     -0.188214882341448E-8,
      0.184587261114837E-1,
     -0.135830407782663E-5,
     -0.723681885626348E17,
     -0.223449194054124E27,
     -0.111526741826431E-34,
      0.276032601145151E-28,
      0.134856491567853E15,
      0.652440293345860E-9,
      0.510655119774360E17,
     -0.468138358908732E32,
     -0.760667491183279E16,
     -0.417247986986821E-18,
      0.312545677756104E14,
     -0.100375333864186E15,
      0.247761392329058E27    
    };
    int Ii[39]={
        -10,-8,-6,-6,-6,-6,-6,-6,-5,-5,-5,-5,-5,-5,-4,-4,-4,-4,
        -3,-3,-3,-2,-2,-1,-1,0,0,0,1,1,3,4,4,4,5,8,10,12,14
    };
    int Ji[39]={
        -8,-12,-12,-3,5,6,8,10,1,2,6,8,10,14,-12,-10,-6,10,-3,
        10,12,2,4,-2,0,-2,6,10,-12,-10,3,-6,3,10,2,-12,-2,-3,1
    };
    double omega=0;
    for(int i=0;i<39;i++){
        omega+=ni[i]*pow(pi-0.960,Ii[i])*pow(theta-0.995,Ji[i]);
    }
    return 0.0031*omega;
}  
double IF97Region3::PT2V3w(double p,double t){
    double pi=p/23.0;
    double theta=(t+T0)/650.0;
    double ni[35]={
     -0.586219133817016E-7,
     -0.894460355005526E11,
      0.531168037519774E-30,
      0.109892402329239,
     -0.575368389425212E-1,
      0.228276853990249E5,
     -0.158548609655002E19,
      0.329865748576503E-27,
     -0.634987981190669E-24,
      0.615762068640611E-8,
     -0.961109240985747E8,
     -0.406274286652625E-44,
     -0.471103725498077E-12,
      0.725937724828145,
      0.187768525763682E-38,
     -0.103308436323771E4,
     -0.662552816342168E-1,
      0.579514041765710E3,
      0.237416732616644E-26,
      0.271700235739893E-14,
     -0.907886213483600E2,
     -0.171242509570207E-36,
      0.156792067854621E3,
      0.923261357901470,
     -0.597865988422577E1,
      0.321988767636389E7,
     -0.399441390042203E-29,
      0.493429086046981E-7,
      0.812036983370565E-19,
     -0.207610284654137E-11,
     -0.340821291419719E-6,
      0.542000573372233E-17,
     -0.856711586510214E-12,
      0.266170454405981E-13,
      0.858133791857099E-5
    };
    int Ii[35]={
        -12,-12,-10,-10,-8,-8,-8,-6,-6,-6,-6,-5,-4,-4,-3,
        -3,-2,-2,-1,-1,-1,0,0,1,2,2,3,3,5,5,5,8,8,10,10
    };
    int Ji[35]={
        8,14,-1,8,6,8,14,-4,-3,2,8,-10,-1,3,-10,3,1,2,-8,
        -4,1,-12,1,-1,-1,2,-12,-5,-10,-8,-6,-12,-10,-12,-8
    };
    double omega=0;
    for(int i=0;i<35;i++){
        omega+=ni[i]*pow(pi-0.959,Ii[i])*pow(theta-0.995,Ji[i]);
    }
    return 0.0039*pow(omega,4);
}  
double IF97Region3::PT2V3x(double p,double t){
    double pi=p/23.0;
    double theta=(t+T0)/650.0;
    double ni[36]={
      0.377373741298151E19,
     -0.507100883722913E13,
     -0.103363225598860E16,
      0.184790814320773E-5,
     -0.924729378390945E-3,
     -0.425999562292738E24,
     -0.462307771873973E-12,
      0.107319065855767E22,
      0.648662492280682E11,
      0.244200600688281E1,
     -0.851535733484258E10,
      0.169894481433592E22,
      0.215780222509020E-26,
     -0.320850551367334,
     -0.382642448458610E17,
     -0.275386077674421E-28,
     -0.563199253391666E6,
     -0.326068646279314E21,
      0.397949001553184E14,
      0.100824008584757E-6,
      0.162234569738433E5,
     -0.432355225319745E11,
     -0.592874245598610E12,
      0.133061647281106E1,
      0.157338197797544E7,
      0.258189614270853E14,
      0.262413209706358E25,
     -0.920011937431142E-1,
      0.220213765905426E-2,
     -0.110433759109547E2,
      0.847004870612087E7,
     -0.592910695762536E9,
     -0.183027173269660E-4,
      0.181339603516302,
     -0.119228759669889E4,
      0.430867658061468E7
    };
    int Ii[36]={
        -8,-6,-5,-4,-4,-4,-3,-3,-1,0,0,0,1,1,2,3,3,3,
        4,5,5,5,6,8,8,8,8,10,12,12,12,12,14,14,14,14
    };
    int Ji[36]={
        14,10,10,1,2,14,-2,12,5,0,4,10,-10,-1,6,-12,0,8,3,
        -6,-2,1,1,-6,-3,1,8,-8,-10,-8,-5,-4,-12,-10,-8,-6
    };
    double omega=0;
    for(int i=0;i<36;i++){
        omega+=ni[i]*pow(pi-0.910,Ii[i])*pow(theta-0.988,Ji[i]);
    }
    return 0.0049*omega;
}
double IF97Region3::PT2V3y(double p,double t){
    double pi=p/22.0;
    double theta=(t+T0)/650.0;
    double ni[20]={
     -0.525597995024633E-9,
      0.583441305228407E4,
     -0.134778968457925E17,
      0.118973500934212E26,
     -0.159096490904708E27,
     -0.315839902302021E-6,
      0.496212197158239E3,
      0.327777227273171E19,
     -0.527114657850696E22,
      0.210017506281863E-16,
      0.705106224399834E21,
     -0.266713136106469E31,
     -0.145370512554562E-7,
      0.149333917053130E28,
     -0.149795620287641E8,
     -0.381881906271100E16,
      0.724660165585797E-4,
     -0.937808169550193E14,
      0.514411468376383E10,
     -0.828198594040141E5
    };
    int Ii[20]={
        0,0,0,0,1,2,2,2,2,3,3,
        3,4,4,5,5,8,8,10,12
    };
    int Ji[20]={
        -3,1,5,8,8,-4,-1,4,5,-8,4,
        8,-6,6,-2,1,-8,-2,-5,-8
    };
    double omega=0;
    for(int i=0;i<20;i++){
        omega+=ni[i]*pow(pi-0.996,Ii[i])*pow(theta-0.994,Ji[i]);
    }
    return 0.0031*pow(omega,4);
}
double IF97Region3::PT2V3z(double p,double t){
    double pi=p/22.0;
    double theta=(t+T0)/650.0;
    double ni[23]={
      0.244007892290650E-10,
     -0.463057430331242E7,
      0.728803274777712E10,
      0.327776302858856E16,
     -0.110598170118409E10,
     -0.323899915729957E13,
      0.923814007023245E16,
      0.842250080413712E-12,
      0.663221436245506E12,
     -0.167170186672139E15,
      0.253749358701391E4,
     -0.819731559610523E-20,
      0.328380587890663E12,
     -0.625004791171543E8,
      0.803197957462023E21,
     -0.204397011338353E-10,
     -0.378391047055938E4,
      0.972876545938620E-2,
      0.154355721681459E2,
     -0.373962862928643E4,
     -0.682859011374572E11,
     -0.248488015614543E-3,
      0.394536049497068E7
    };
    int Ii[23]={
        -8,-6,-5,-5,-4,-4,-4,-3,-3,-3,
        -2,-1,0,1,2,3,3,6,6,6,6,8,8
    };
    int Ji[23]={
        3,6,6,8,5,6,8,-2,5,6,2,-6,3,
        1,6,-6,-2,-6,-5,-4,-1,-8,-4
    };
    double omega=0;
    for(int i=0;i<23;i++){
        omega+=ni[i]*pow(pi-0.993,Ii[i])*pow(theta-0.994,Ji[i]);
    }
    return 0.0038*pow(omega,4);
}
void IF97Region3::testPT(){
    double p,t,v,vv;
    p=50;	t=630-273.15;	v=0.0014708531;	    vv=PT2V3a(p,t);	cout<<setprecision(12)<<"3a"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=80;	t=670-273.15;	v=0.001503831359;	vv=PT2V3a(p,t);	cout<<setprecision(12)<<"3a"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=50;	t=710-273.15;	v=0.002204728587;	vv=PT2V3b(p,t);	cout<<setprecision(12)<<"3b"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=80;	t=750-273.15;	v=0.00197369294;	vv=PT2V3b(p,t);	cout<<setprecision(12)<<"3b"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=20;	t=630-273.15;	v=0.001761696406;	vv=PT2V3c(p,t);	cout<<setprecision(12)<<"3c"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=30;	t=650-273.15;	v=0.001819560617;	vv=PT2V3c(p,t);	cout<<setprecision(12)<<"3c"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=26;	t=656-273.15;	v=0.00224558772;	vv=PT2V3d(p,t);	cout<<setprecision(12)<<"3d"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=30;	t=670-273.15;	v=0.002506897702;	vv=PT2V3d(p,t);	cout<<setprecision(12)<<"3d"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=26;	t=661-273.15;	v=0.002970225962;	vv=PT2V3e(p,t);	cout<<setprecision(12)<<"3e"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=30;	t=675-273.15;	v=0.003004627086;	vv=PT2V3e(p,t);	cout<<setprecision(12)<<"3e"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=26;	t=671-273.15;	v=0.005019029401;	vv=PT2V3f(p,t);	cout<<setprecision(12)<<"3f"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=30;	t=690-273.15;	v=0.004656470142;	vv=PT2V3f(p,t);	cout<<setprecision(12)<<"3f"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=23.6;	t=649-273.15;	v=0.002163198378;	vv=PT2V3g(p,t);	cout<<setprecision(12)<<"3g"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=24;	t=650-273.15;	v=0.002166044161;	vv=PT2V3g(p,t);	cout<<setprecision(12)<<"3g"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=23.6;	t=652-273.15;	v=0.002651081407;	vv=PT2V3h(p,t);	cout<<setprecision(12)<<"3h"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=24;	t=654-273.15;	v=0.002967802335;	vv=PT2V3h(p,t);	cout<<setprecision(12)<<"3h"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=23.6;	t=653-273.15;	v=0.003273916816;	vv=PT2V3i(p,t);	cout<<setprecision(12)<<"3i"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=24;	t=655-273.15;	v=0.003550329864;	vv=PT2V3i(p,t);	cout<<setprecision(12)<<"3i"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=23.5;	t=655-273.15;	v=0.004545001142;	vv=PT2V3j(p,t);	cout<<setprecision(12)<<"3j"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=24;	t=660-273.15;	v=0.005100267704;	vv=PT2V3j(p,t);	cout<<setprecision(12)<<"3j"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=23;	t=660-273.15;	v=0.006109525997;	vv=PT2V3k(p,t);	cout<<setprecision(12)<<"3k"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=24;	t=670-273.15;	v=0.006427325645;	vv=PT2V3k(p,t);	cout<<setprecision(12)<<"3k"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.6;	t=646-273.15;	v=0.002117860851;	vv=PT2V3l(p,t);	cout<<setprecision(12)<<"3l"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=23;	t=646-273.15;	v=0.002062374674;	vv=PT2V3l(p,t);	cout<<setprecision(12)<<"3l"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.6;	t=648.6-273.15;	v=0.00253306378;	vv=PT2V3m(p,t);	cout<<setprecision(12)<<"3m"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.8;	t=649.3-273.15;	v=0.002572971781;	vv=PT2V3m(p,t);	cout<<setprecision(12)<<"3m"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.6;	t=649-273.15;	v=0.002923432711;	vv=PT2V3n(p,t);	cout<<setprecision(12)<<"3n"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.8;	t=649.7-273.15;	v=0.002913311494;	vv=PT2V3n(p,t);	cout<<setprecision(12)<<"3n"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.6;	t=649.1-273.15;	v=0.003131208996;	vv=PT2V3o(p,t);	cout<<setprecision(12)<<"3o"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.8;	t=649.9-273.15;	v=0.003221160278;	vv=PT2V3o(p,t);	cout<<setprecision(12)<<"3o"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.6;	t=649.4-273.15;	v=0.003715596186;	vv=PT2V3p(p,t);	cout<<setprecision(12)<<"3p"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.8;	t=650.2-273.15;	v=0.00366475479;	vv=PT2V3p(p,t);	cout<<setprecision(12)<<"3p"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=21.1;	t=640-273.15;	v=0.001970999272;	vv=PT2V3q(p,t);	cout<<setprecision(12)<<"3q"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=21.8;	t=643-273.15;	v=0.002043919161;	vv=PT2V3q(p,t);	cout<<setprecision(12)<<"3q"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=21.1;	t=644-273.15;	v=0.005251009921;	vv=PT2V3r(p,t);	cout<<setprecision(12)<<"3r"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=21.8;	t=648-273.15;	v=0.005256844741;	vv=PT2V3r(p,t);	cout<<setprecision(12)<<"3r"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=19.1;	t=635-273.15;	v=0.001932829079;	vv=PT2V3s(p,t);	cout<<setprecision(12)<<"3s"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=20;	t=638-273.15;	v=0.001985387227;	vv=PT2V3s(p,t);	cout<<setprecision(12)<<"3s"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=17;	t=626-273.15;	v=0.008483262001;	vv=PT2V3t(p,t);	cout<<setprecision(12)<<"3t"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=20;	t=640-273.15;	v=0.006227528101;	vv=PT2V3t(p,t);	cout<<setprecision(12)<<"3t"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=21.5;	t=644.6-273.15;	v=2.268366647E-3;	vv=PT2V3u(p,t);	cout<<setprecision(10)<<"3u"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.0;	t=646.1-273.15;	v=2.296350553E-3;	vv=PT2V3u(p,t);	cout<<setprecision(10)<<"3u"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.5;	t=648.6-273.15;	v=2.832373260E-3;	vv=PT2V3v(p,t);	cout<<setprecision(10)<<"3v"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.3;	t=647.9-273.15;	v=2.811424405E-3;	vv=PT2V3v(p,t);	cout<<setprecision(10)<<"3v"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.15;	t=647.5-273.15;	v=3.694032281E-3;	vv=PT2V3w(p,t);	cout<<setprecision(10)<<"3w"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.3;	t=648.1-273.15;	v=3.622226305E-3;	vv=PT2V3w(p,t);	cout<<setprecision(10)<<"3w"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.11;	t=648.0-273.15;	v=4.528072649E-3;	vv=PT2V3x(p,t);	cout<<setprecision(10)<<"3x"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.3;	t=649.0-273.15;	v=4.556905799E-3;	vv=PT2V3x(p,t);	cout<<setprecision(10)<<"3x"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.0;	t=646.84-273.15;	v=2.698354719E-3;	vv=PT2V3y(p,t);	cout<<setprecision(10)<<"3y"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.064;	t=647.05-273.15;	v=2.717655648E-3;	vv=PT2V3y(p,t);	cout<<setprecision(10)<<"3y"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.0;	t=646.89-273.15;	v=3.798732962E-3;	vv=PT2V3z(p,t);	cout<<setprecision(10)<<"3z"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
    p=22.064;	t=647.15-273.15;	v=3.701940010E-3;	vv=PT2V3z(p,t);	cout<<setprecision(10)<<"3z"<<"\t"<<v<<"\t"<<vv<<"\t"<<abs(100*(vv-v)/v)<<endl;
}
/*end P,T 2 辅助函数*/
void verifyPT(){
    double p,t,v,vv,pp,tt;
    int i=0,max=0;
    vv=IF97Region3::PT2V(21.9547685800,373,i);
    vv=IF97Region3::PT2V(22.2930643,376.85,i);
    vv=IF97Region3::PT2V(25.5837018,376.85,i);
    vv=IF97Region3::PT2V(78.3095639,476.85,i); 
    double dp=(100-IF97Region3::T2P_B23(350))/500;
    for(p=IF97Region3::T2P_B23(350);p<=100;p+=dp){
        double maxt=IF97Region3::P2T_B23(p);
        double dt=0.1;
        for(t=350;t<=maxt;t+=dt){
            v=IF97Region3::PT2V(p,t,i);
            pp=IF97Region3::TV2P(t,v);
            if(i>max) max=i;
            cout<<setprecision(10)<<i<<'\t'<<v<<"\t"<<t<<"\t"<<p<<"\t"<<pp<<"\t";
            cout<<setprecision(10)<<abs(100*(pp-p)/p)<<"\t"<<endl;
        }
    }
    cout<<max<<endl;
}
void verifyPH(){
    double p,t,v,vv,h,hh,pp,tt;
    int i=0,max=0;
    IF97Region3::PH2TV(25.5837018,1863.43019,t,v,i);
    IF97Region3::PH2TV(0.222930643E2,0.237512401E4,t,v,i);
    IF97Region3::PH2TV(0.783095639E2,0.225868845E4,t,v,i);
    v=IF97Region3::PT2V(83.945184,451.5,i);
    p=IF97Region3::TV2P(425.5,v);
    h=IF97Region3::TV2H(425.5,v);
    pp=IF97Region3::TV2P(425.4964925,0.001874206903);
    IF97Region3::PH2TV(56.511712,2021.6512,t,v,i);
    IF97Region3::PH2TV(43.730208,2584.16687,t,v,i);
    double dp=(100-22.064)/500;
    for(p=22.064;p<=100.00001;p+=dp){
        double maxt=IF97Region3::P2T_B23(p);
        double dt=0.1;
        for(t=350;t<=maxt;t+=dt){
            v=IF97Region3::PT2V(p,t,i);
            h=IF97Region3::TV2H(t,v);
            IF97Region3::PH2TV(p,h,tt,vv,i);
            if(i>max) max=i;
            cout<<setprecision(10)<<i<<'\t'<<h<<"\t"<<v<<"\t"<<vv<<"\t"<<p<<"\t"<<t<<"\t";
            cout<<setprecision(10)<<tt<<"\t"<<abs(100*(tt-t)/t)<<"\t"<<endl;
        }
    }
    cout<<max<<endl;
}
void verifyPS(){
    double p,t,v,s,vv,tt;
    int i=0,max=0;
    double dp=(100-22.064)/500;
    for(p=22.064;p<=100;p+=dp){
        double maxt=IF97Region3::P2T_B23(p);
        double dt=0.1;
        for(t=350;t<=maxt;t+=dt){
            v=IF97Region3::PT2V(p,t,i);
            s=IF97Region3::TV2S(t,v);
            IF97Region3::PS2TV(p,s,tt,vv,i);
            if(i>max) max=i;
            cout<<setprecision(10)<<i<<'\t'<<s<<"\t"<<v<<"\t"<<vv<<"\t"<<p<<"\t"<<t<<"\t";
            cout<<setprecision(10)<<tt<<"\t"<<abs(100*(tt-t)/t)<<"\t"<<endl;
        }
    }
    cout<<max<<endl;
}
void verifyHS(double p0,double p1,int divs){
    double p,t,v,h,s,vv,tt,pp;
    int i=0,max=0;
    double dp=(p1-p0)/divs;
    for(p=p0;p<=p1+0.000001;p+=dp){
        double maxt=IF97Region3::P2T_B23(p);
        double dt=0.01;
        for(t=350;t<=maxt;t+=dt){
            v=IF97Region3::PT2V(p,t,i);
            h=IF97Region3::TV2H(t,v);
            s=IF97Region3::TV2S(t,v);
            IF97Region3::HS2TVP(h,s,tt,vv,pp,i);
            if(i>max) max=i;
            cout<<setprecision(10)<<i<<'\t'<<v<<"\t"<<vv<<"\t"<<abs((vv-v))<<"\t";
            cout<<setprecision(10)<<p<<"\t"<<pp<<"\t"<<abs(pp-p)<<"\t";
            cout<<setprecision(10)<<t<<"\t"<<tt<<"\t"<<abs((tt-t))<<endl;
        }
    }
    cout<<max<<endl;
}
int main(int argc, char *argv[]){ 
    double p,t,v,h,s,tt,vv,pp;
    int i;
    // v=1.0/200;
    // t=650-273.15;
    // pp=IF97Region3::TV2P(t,v);  
    // //verifyPT();
    // //verifyPH();
    // //verifyPS();
    p=97.506048;    
    t=374.1;
    v=IF97Region3::PT2V(p,t,i);
    h=IF97Region3::TV2H(t,v);
    s=IF97Region3::TV2S(t,v);
    IF97Region3::HS2TVP(h,s,tt,vv,pp,i);
    if(argc==4){        
        double p0=atof(argv[1]);
        double p1=atof(argv[2]);
        int divs=atoi(argv[3]);
        time_t start_t, end_t;
        double diff_t;
        time(&start_t);
        verifyHS(p0,p1,divs);
        time(&end_t);
        diff_t = difftime(end_t, start_t);
        cout<<"duration:"<<diff_t<<endl;
    }
    return 0;
} 
