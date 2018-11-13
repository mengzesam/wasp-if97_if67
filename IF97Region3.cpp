#include "IF97Region3.h"
#include "IF97Region4.h"
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
    return 0.0;
}
void verifyHS(){
    // double p,t,h,s,pp,tt;
    // int i=0,max=0;
    // double dp=(100-0.000611657)/200;
    // for(p=0.000611657;p<=100;p+=dp){
    //     double maxt=350;
    //     if(p<=22.064){
    //         maxt=IF97Region4::P2T(p);
    //         maxt=maxt>350?350:maxt;
    //     }
    //     double dt=maxt/200;
    //     for(t=0;t<=maxt;t+=dt){
    //         h=IF97Region3::PT2H(p,t);
    //         s=IF97Region3::PT2S(p,t);
    //         IF97Region3::HS2PT(h,s,pp,tt,i);
    //         if(i>max) max=i;
    //         cout<<setprecision(10)<<i<<'\t'<<h<<"\t"<<s<<"\t"<<p<<"\t"<<t<<"\t";
    //         cout<<setprecision(10)<<pp<<'\t'<<tt<<"\t"<<abs(100*(pp-p)/p)<<"\t";
    //         if(t<1E-9)
    //             cout<<setprecision(10)<<abs(tt-t)<<endl;
    //         else
    //             cout<<setprecision(10)<<abs(100.0*(tt-t)/t)<<endl;
    //     }
    // }
    // cout<<max<<endl;
}
int main(){ 
    double p,t,h,u,s,v,cp,cv,pp,tt,p2;
    int i;
    v=1.0/200;
    t=650-273.15;
    p=IF97Region3::TV2P(t,v);
    h=IF97Region3::TV2H(t,v);
    u=IF97Region3::TV2U(t,v);
    s=IF97Region3::TV2S(t,v);
    cp=IF97Region3::TV2Cp(t,v);
    cv=IF97Region3::TV2Cv(t,v);    
    return 0;
} 
