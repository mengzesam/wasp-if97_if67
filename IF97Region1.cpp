#include "IF97Region1.h"
#include "IF97Region4.h"
#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

double IF97Region1::PT2H(double p,double t){
    double pi=p/p_base;
    double tau=T_base/(t+T0);
    double gamma_tau=0;
    for(int i = 0;i<34; i++)
        gamma_tau+=ni[i]*pow(7.1-pi,Ii[i])*Ji[i]*pow(tau-1.222,Ji[i]-1);
    return tau*gamma_tau*(t+T0)*R;
}
double IF97Region1::PT2S(double p,double t){
    double pi=p/p_base;
    double tau=T_base/(t+T0);
    double gamma_tau=0;
    for(int i = 0;i<34; i++)
        gamma_tau+=ni[i]*pow(7.1-pi,Ii[i])*Ji[i]*pow(tau-1.222,Ji[i]-1);
    double gamma=0;
    for(int i=0;i<34;i++)
		gamma+=ni[i]*pow(7.1-pi,Ii[i])*pow(tau-1.222,Ji[i]);
    return (tau*gamma_tau-gamma)*R;
}
double IF97Region1::PT2V(double p,double t){
    double pi=p/p_base;
    double tau=T_base/(t+T0);
    double gamma_pi=0;
    for(int i=0;i<34;i++)
        gamma_pi+=-ni[i]*Ii[i]*pow(7.1-pi,Ii[i]-1)*pow(tau-1.222,Ji[i]);
    return pi*gamma_pi*(t+T0)*R/(p*1000);
}
double IF97Region1::PT2U(double p,double t){
    double pi=p/p_base;
    double tau=T_base/(t+T0);
    double gamma_tau=0;
    for(int i = 0;i<34; i++)
        gamma_tau+=ni[i]*pow(7.1-pi,Ii[i])*Ji[i]*pow(tau-1.222,Ji[i]-1);
    double gamma_pi=0;
    for(int i=0;i<34;i++)
        gamma_pi+=-ni[i]*Ii[i]*pow(7.1-pi,Ii[i]-1)*pow(tau-1.222,Ji[i]);
	return (tau*gamma_tau-pi*gamma_pi)*(t+T0)*R;
}
double IF97Region1::PT2Cp(double p,double t){
    double pi=p/p_base;
    double tau=T_base/(t+T0);
    double gamma_tau2=0;
    for(int i = 0;i<34; i++)
		gamma_tau2+=ni[i]*pow(7.1-pi,Ii[i])*Ji[i]*(Ji[i]-1)*pow(tau-1.222,Ji[i]-2);    
    return -tau*tau*gamma_tau2*R;
}
double IF97Region1::PT2Cv(double p,double t){//TO check
    double pi=p/p_base;
    double tau=T_base/(t+T0);
    double gamma_tau2=0;
    for(int i = 0;i<34; i++)
		gamma_tau2+=ni[i]*pow(7.1-pi,Ii[i])*Ji[i]*(Ji[i]-1)*pow(tau-1.222,Ji[i]-2);
    double gamma_pi=0;
    for(int i=0;i<34;i++)
        gamma_pi+=-ni[i]*Ii[i]*pow(7.1-pi,Ii[i]-1)*pow(tau-1.222,Ji[i]);
    double gamma_pi2=0;
    for(int i = 0;i<34; i++)
        gamma_pi2+=ni[i]*Ii[i]*(Ii[i]-1)*pow(7.1-pi,Ii[i]-2)*pow(tau-1.222,Ji[i]);
    double gamma_pitau=0;
    for(int i = 0;i<34; i++)
		gamma_pitau+=-ni[i]*Ii[i]*pow(7.1-pi,Ii[i]-1)*Ji[i]*pow(tau-1.222,Ji[i]-1);      
    return (-tau*tau*gamma_tau2+(gamma_pi-tau*gamma_pitau)*(gamma_pi-tau*gamma_pitau)/gamma_pi2)*R;
}
double IF97Region1::PT2W(double p,double t){//TODO
    return 100.0;
}
double IF97Region1::PH2T(double p,double h){
    double P_base=1; //p* 1Mpa
    double T_base=1; //T* 1K
    double H_base=2500; //h* 2500kJ/kg
    double pi=p/P_base;
    double eta=h/H_base;
    double theta=0;
    for(int i =0;i<20; i++) {
        theta+=ni_tab6[i]*pow(pi,Ii_tab6[i])*pow(eta+1,Ji_tab6[i]);
    }
    double t=theta*T_base-T0;//先用backward方程计算t初值，再用弦截法求精确的t
    double hh=PT2H(p,t);
    if(abs(hh-h)>ERR){
        double t0=t;
        double h0=hh;
        t=t0+0.2;//需要判断是否超出region1区域
        hh=PT2H(p,t);
        double t1=t;
        double h1=hh;
        t=t1+(h-h1)/(h0-h1)*(t0-t1);
        hh=PT2H(p,t);
        int i=0;
        while(abs(hh-h)>ERR){
            i++;
            t0=t1;
            h0=h1;
            t1=t;
            h1=hh;
            t=t1+(h-h1)/(h0-h1)*(t0-t1);
            hh=PT2H(p,t);
        }
        cout<<i<<endl;
    }
    return t;    
}
double IF97Region1::PH2S(double p,double h){
    double t=PH2T(p,h);
    return PT2S(p,t);
}
double IF97Region1::PH2V(double p,double h){
    double t=PH2T(p,h);
    return PT2V(p,t);
}
double IF97Region1::PH2U(double p,double h){
    double t=PH2T(p,h);
    return PT2U(p,t);
}
double IF97Region1::PH2Cp(double p,double h){
    double t=PH2T(p,h);
    return PT2Cp(p,t);
}
double IF97Region1::PH2Cv(double p,double h){
    double t=PH2T(p,h);
    return PT2Cv(p,t);
}
double IF97Region1::PH2W(double p,double h){
    double t=PH2T(p,h);
    return PT2W(p,t);
}
double IF97Region1::PS2T(double p,double s){
    double P_base=1; //p* 1Mpa
    double T_base=1; //T* 1K
    double S_base=1; //S* 1kJ/kg.K	
    double pi=p/P_base;
    double sigma=s/S_base;
    double theta=0;
    for(int i =0;i<20; i++) {
        theta+=ni_tab8[i]*pow(pi,Ii_tab8[i])*pow(sigma+2,Ji_tab8[i]);
    }
    double t=theta*T_base-T0;//先用backward方程计算t初值，再用弦截法求精确的t
    double ss=PT2S(p,t);
    if(abs(ss-s)>ERR){
        double t0=t;
        double s0=ss;
        t=t0+0.2;//需要判断是否超出region1区域
        ss=PT2S(p,t);
        double t1=t;
        double s1=ss;
        t=t1+(s-s1)/(s0-s1)*(t0-t1);
        ss=PT2S(p,t);
        //int i=0;
        while(abs(ss-s)>ERR){
        //    i++;
            t0=t1;
            s0=s1;
            t1=t;
            s1=ss;
            t=t1+(s-s1)/(s0-s1)*(t0-t1);
            ss=PT2S(p,t);
        }
        //cout<<i<<endl;
    }
    return t;    
}
double IF97Region1::PS2H(double p,double s){
    double t=PS2T(p,s);
    return PT2H(p,t);
}
double IF97Region1::PS2V(double p,double s){
    double t=PS2T(p,s);
    return PT2V(p,t);
}
double IF97Region1::PS2U(double p,double s){
    double t=PS2T(p,s);
    return PT2U(p,t);
}
double IF97Region1::PS2Cp(double p,double s){
    double t=PS2T(p,s);
    return PT2Cp(p,t);
}
double IF97Region1::PS2Cv(double p,double s){
    double t=PS2T(p,s);
    return PT2Cv(p,t);
}
double IF97Region1::PS2W(double p,double s){
    double t=PS2T(p,s);
    return PT2W(p,t);
} 
double IF97Region1::PV2T(double p,double v){
    double pi=p/p_base;
    double lefts[34];
    for(int i = 0;i<34; i++)
        lefts[i]=-ni[i]*Ii[i]*pow(7.1-pi,Ii[i]-1);    
    double t=175;//TODO:初值选择再优化
    double tl=IF97Region4::P2T(p)+1;
    double tr=350;
    if(p>16.5){
        tl=348;
        tr=350;
    }
    double tau=T_base/(t+T0);
    double gamma_pi=0;
    for(int i=0;i<34;i++)
        gamma_pi+=lefts[i]*pow(tau-1.222,Ji[i]);
    double vv=pi*gamma_pi*(t+T0)*R/(p*1000);
    if(abs(vv-v)>ERR2){
        double t0=t;
        double v0=vv;
        t=t0+80;//TODO:初值选择再优化
        tau=T_base/(t+T0);
        gamma_pi=0;
        for(int i=0;i<34;i++)
            gamma_pi+=lefts[i]*pow(tau-1.222,Ji[i]);
        vv=pi*gamma_pi*(t+T0)*R/(p*1000);
        double t1=t;
        double v1=vv;
        t=t1+(v-v1)/(v0-v1)*(t0-t1);
        tau=T_base/(t+T0);
        gamma_pi=0;
        for(int i=0;i<34;i++)
            gamma_pi+=lefts[i]*pow(tau-1.222,Ji[i]);
        vv=pi*gamma_pi*(t+T0)*R/(p*1000);
        int i=0;
        while(abs(vv-v)>ERR2){
            i++;
            t0=t1;
            v0=v1;
            t1=t;
            v1=vv;
            t=t1+(v-v1)/(v0-v1)*(t0-t1);
            tau=T_base/(t+T0);
            gamma_pi=0;
            for(int i=0;i<34;i++)
                gamma_pi+=lefts[i]*pow(tau-1.222,Ji[i]);
            vv=pi*gamma_pi*(t+T0)*R/(p*1000);
        }
        cout<<i<<endl;  
    }
    return t;
}



int main(){
    double p,t,h,s,v;
    cout<<"p\tt\th\ts\tv"<<endl;
    for(t=5;t<350.1;t+=15){
        p=IF97Region4::T2P(t)+0.00001;
        double dp=(100-p)/20;
        while(p<100.1){
            v=1000*IF97Region1::PT2V(p,t);
            h=IF97Region1::PT2H(p,t);
            s=IF97Region1::PT2S(p,t);
            cout<<setprecision(10)<<p<<"\t"<<t<<"\t"<<h<<"\t"<<s<<"\t"<<v<<endl;
            p+=dp;
        }
    }
}
