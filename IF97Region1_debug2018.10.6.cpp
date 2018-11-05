#include "IF97Region1.h"
#include "IF97Region4.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <time.h>
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
double IF97Region1::PH2T(double p,double h,int& itera){
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
    itera=0;
    if(abs(hh-h)>ERR0){
        double t0=t;
        double h0=hh;
        t=t0-0.15;//需要判断是否超出region1区域
        t=t<0?t0+0.15:t;
        hh=PT2H(p,t);
        double t1=t;
        double h1=hh;
        t=t1+(h-h1)/(h0-h1)*(t0-t1);
        hh=PT2H(p,t);
        while(abs(hh-h)>ERR0){
            itera++;
            t0=t1;
            h0=h1;
            t1=t;
            h1=hh;
            t=t1+(h-h1)/(h0-h1)*(t0-t1);
            hh=PT2H(p,t);
        }
    }
    return t;    
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
    if(abs(hh-h)>ERR0){
        double t0=t;
        double h0=hh;
        t=t0-0.15;//需要判断是否超出region1区域
        t=t<0?t0+0.15:t;
        hh=PT2H(p,t);
        double t1=t;
        double h1=hh;
        t=t1+(h-h1)/(h0-h1)*(t0-t1);
        hh=PT2H(p,t);
        while(abs(hh-h)>ERR0){
            t0=t1;
            h0=h1;
            t1=t;
            h1=hh;
            t=t1+(h-h1)/(h0-h1)*(t0-t1);
            hh=PT2H(p,t);
        }
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
double IF97Region1::PS2T(double p,double s,int& itera){
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
    itera=0;
    if(abs(ss-s)>ERR){
        double t0=t;
        double s0=ss;
        t=t0-0.15;//需要判断是否超出region1区域
        t=t<0?t0+0.15:t;
        ss=PT2S(p,t);
        double t1=t;
        double s1=ss;
        t=t1+(s-s1)/(s0-s1)*(t0-t1);
        ss=PT2S(p,t);
        while(abs(ss-s)>ERR){
            itera++;
            t0=t1;
            s0=s1;
            t1=t;
            s1=ss;
            t=t1+(s-s1)/(s0-s1)*(t0-t1);
            ss=PT2S(p,t);
        }
    }
    return t;    
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
        t=t0-0.15;//需要判断是否超出region1区域
        t=t<0?t0+0.15:t;
        ss=PT2S(p,t);
        double t1=t;
        double s1=ss;
        t=t1+(s-s1)/(s0-s1)*(t0-t1);
        ss=PT2S(p,t);
        while(abs(ss-s)>ERR){
            t0=t1;
            s0=s1;
            t1=t;
            s1=ss;
            t=t1+(s-s1)/(s0-s1)*(t0-t1);
            ss=PT2S(p,t);
        }
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
double IF97Region1::PV2T(double p,double v,int& itera){
    double pi=p/p_base;
 /*.利用拟合的粗略公式设置初值t*/
    //double t=0.00038689898736603837*p*p-1277125254.350867*v*v+209.65767720269834*p*v
    //       +0.334694858543311*p+3762351.523728678*v-2447.7223094011324;
    double t=0.00027955518134947083*p*p-1341269530.928515*v*v+212.1875236957701*p*v
             +0.3412237701411934*p+3905702.229360094*v-2526.3683009658025;
    double lefts[34];
    for(int i = 0;i<34; i++)
        lefts[i]=-ni[i]*Ii[i]*pow(7.1-pi,Ii[i]-1);    
    double tau=T_base/(t+T0);
    double gamma_pi=0;
    for(int i=0;i<34;i++)
        gamma_pi+=lefts[i]*pow(tau-1.222,Ji[i]);
    double vv=pi*gamma_pi*(t+T0)*R/(p*1000);
    itera=0;
    if(abs(vv-v)>ERR2){
        double t0=t;
        double v0=vv;
        t=IF97Region4::P2T(p)-15;//TODO:初值选择再优化
        t=(abs(t-t0)<15)?t0-15:t;    
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
        while(abs(vv-v)>ERR2){
            itera++;
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
    }
    return t;
}
double IF97Region1::PV2T(double p,double v){
    double pi=p/p_base;
 /*.利用拟合的粗略公式设置初值t*/
    //double t=0.00038689898736603837*p*p-1277125254.350867*v*v+209.65767720269834*p*v
    //       +0.334694858543311*p+3762351.523728678*v-2447.7223094011324;
    double t=0.00027955518134947083*p*p-1341269530.928515*v*v+212.1875236957701*p*v
             +0.3412237701411934*p+3905702.229360094*v-2526.3683009658025;
    double lefts[34];
    for(int i = 0;i<34; i++)
        lefts[i]=-ni[i]*Ii[i]*pow(7.1-pi,Ii[i]-1);    
    double tau=T_base/(t+T0);
    double gamma_pi=0;
    for(int i=0;i<34;i++)
        gamma_pi+=lefts[i]*pow(tau-1.222,Ji[i]);
    double vv=pi*gamma_pi*(t+T0)*R/(p*1000);
    if(abs(vv-v)>ERR2){
        double t0=t;
        double v0=vv;
        t=IF97Region4::P2T(p)-15;//TODO:初值选择再优化
        t=(abs(t-t0)<15)?t0-15:t;    
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
        while(abs(vv-v)>ERR2){
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
    }
    return t;
} 
double IF97Region1::PV2H(double p,double v){
    double t=PV2T(p,v);
    return PT2H(p,t);
}
double IF97Region1::PV2S(double p,double v){
    double t=PV2T(p,v);
    return PT2S(p,t);
}
double IF97Region1::PV2U(double p,double v){
    double t=PV2T(p,v);
    return PT2U(p,t);
}
double IF97Region1::PV2Cp(double p,double v){
    double t=PV2T(p,v);
    return PT2Cp(p,t);
}
double IF97Region1::PV2Cv(double p,double v){
    double t=PV2T(p,v);
    return PT2Cv(p,t);
}
double IF97Region1::PV2W(double p,double v){
    double t=PV2T(p,v);
    return PT2W(p,t);
}
double IF97Region1::PU2T(double p,double u,int& itera){
    double pi=p/p_base;
 /*.利用拟合的粗略公式设置初值t*/
    double t=-0.0003143081885165228*p*p+(-1.7372412480879492E-05)*u*u
            +0.0002583034211507768*p*u-0.013044871270756663*p
            +0.24745528414601214*u-0.27562996918626365;
    double lefts_tau[34];
    for(int i = 0;i<34; i++)
        lefts_tau[i]=ni[i]*pow(7.1-pi,Ii[i])*Ji[i];
    double lefts_pi[34];
    for(int i = 0;i<34; i++)
        lefts_pi[i]=-ni[i]*Ii[i]*pow(7.1-pi,Ii[i]-1);    
    double tau=T_base/(t+T0);
    double gamma_tau=0;
    for(int i = 0;i<34; i++)
        gamma_tau+=lefts_tau[i]*pow(tau-1.222,Ji[i]-1);
    double gamma_pi=0;
    for(int i=0;i<34;i++)
        gamma_pi+=lefts_pi[i]*pow(tau-1.222,Ji[i]);
    double uu=(tau*gamma_tau-pi*gamma_pi)*(t+T0)*R;
    itera=0;
    if(abs(uu-u)>ERR){
        double t0=t;
        double u0=uu;
        t=IF97Region4::P2T(p)-15;//TODO:初值选择再优化
        t=(abs(t-t0)<15)?t0-15:t; 
        tau=T_base/(t+T0);
        gamma_tau=0;
        for(int i = 0;i<34; i++)
            gamma_tau+=lefts_tau[i]*pow(tau-1.222,Ji[i]-1);
        gamma_pi=0;
        for(int i=0;i<34;i++)
            gamma_pi+=lefts_pi[i]*pow(tau-1.222,Ji[i]);
        uu=(tau*gamma_tau-pi*gamma_pi)*(t+T0)*R;
        double t1=t;
        double u1=uu;
        t=t1+(u-u1)/(u0-u1)*(t0-t1);
        tau=T_base/(t+T0);
        gamma_tau=0;
        for(int i = 0;i<34; i++)
            gamma_tau+=lefts_tau[i]*pow(tau-1.222,Ji[i]-1);
        gamma_pi=0;
        for(int i=0;i<34;i++)
            gamma_pi+=lefts_pi[i]*pow(tau-1.222,Ji[i]);
        uu=(tau*gamma_tau-pi*gamma_pi)*(t+T0)*R;
        while(abs(uu-u)>ERR){
            itera++;
            t0=t1;
            u0=u1;
            t1=t;
            u1=uu;
            t=t1+(u-u1)/(u0-u1)*(t0-t1);
            tau=T_base/(t+T0);
            gamma_tau=0;
            for(int i = 0;i<34; i++)
                gamma_tau+=lefts_tau[i]*pow(tau-1.222,Ji[i]-1);
            gamma_pi=0;
            for(int i=0;i<34;i++)
                gamma_pi+=lefts_pi[i]*pow(tau-1.222,Ji[i]);
            uu=(tau*gamma_tau-pi*gamma_pi)*(t+T0)*R;
        }
    }
    return t;
}
double IF97Region1::PU2T(double p,double u){
    double pi=p/p_base;
 /*.利用拟合的粗略公式设置初值t*/
    double t=-0.0003143081885165228*p*p+(-1.7372412480879492E-05)*u*u
            +0.0002583034211507768*p*u-0.013044871270756663*p
            +0.24745528414601214*u-0.27562996918626365;
    double lefts_tau[34];
    for(int i = 0;i<34; i++)
        lefts_tau[i]=ni[i]*pow(7.1-pi,Ii[i])*Ji[i];
    double lefts_pi[34];
    for(int i = 0;i<34; i++)
        lefts_pi[i]=-ni[i]*Ii[i]*pow(7.1-pi,Ii[i]-1);    
    double tau=T_base/(t+T0);
    double gamma_tau=0;
    for(int i = 0;i<34; i++)
        gamma_tau+=lefts_tau[i]*pow(tau-1.222,Ji[i]-1);
    double gamma_pi=0;
    for(int i=0;i<34;i++)
        gamma_pi+=lefts_pi[i]*pow(tau-1.222,Ji[i]);
    double uu=(tau*gamma_tau-pi*gamma_pi)*(t+T0)*R;
    if(abs(uu-u)>ERR){
        double t0=t;
        double u0=uu;
        t=IF97Region4::P2T(p)-15;//TODO:初值选择再优化
        t=(abs(t-t0)<15)?t0-15:t; 
        tau=T_base/(t+T0);
        gamma_tau=0;
        for(int i = 0;i<34; i++)
            gamma_tau+=lefts_tau[i]*pow(tau-1.222,Ji[i]-1);
        gamma_pi=0;
        for(int i=0;i<34;i++)
            gamma_pi+=lefts_pi[i]*pow(tau-1.222,Ji[i]);
        uu=(tau*gamma_tau-pi*gamma_pi)*(t+T0)*R;
        double t1=t;
        double u1=uu;
        t=t1+(u-u1)/(u0-u1)*(t0-t1);
        tau=T_base/(t+T0);
        gamma_tau=0;
        for(int i = 0;i<34; i++)
            gamma_tau+=lefts_tau[i]*pow(tau-1.222,Ji[i]-1);
        gamma_pi=0;
        for(int i=0;i<34;i++)
            gamma_pi+=lefts_pi[i]*pow(tau-1.222,Ji[i]);
        uu=(tau*gamma_tau-pi*gamma_pi)*(t+T0)*R;
        while(abs(uu-u)>ERR){
            t0=t1;
            u0=u1;
            t1=t;
            u1=uu;
            t=t1+(u-u1)/(u0-u1)*(t0-t1);
            tau=T_base/(t+T0);
            gamma_tau=0;
            for(int i = 0;i<34; i++)
                gamma_tau+=lefts_tau[i]*pow(tau-1.222,Ji[i]-1);
            gamma_pi=0;
            for(int i=0;i<34;i++)
                gamma_pi+=lefts_pi[i]*pow(tau-1.222,Ji[i]);
            uu=(tau*gamma_tau-pi*gamma_pi)*(t+T0)*R;
        }  
    }
    return t;
}
double IF97Region1::PU2H(double p,double u){
    double t=PU2T(p,u);
    return PT2H(p,t);
}
double IF97Region1::PU2S(double p,double u){
    double t=PU2T(p,u);
    return PT2S(p,t);
}
double IF97Region1::PU2V(double p,double u){
    double t=PU2T(p,u);
    return PT2V(p,t);
}
double IF97Region1::PU2Cp(double p,double u){
    double t=PU2T(p,u);
    return PT2Cp(p,t);
}
double IF97Region1::PU2Cv(double p,double u){
    double t=PU2T(p,u);
    return PT2Cv(p,t);
}
double IF97Region1::PU2W(double p,double u){
    double t=PU2T(p,u);
    return PT2W(p,t);
}
double IF97Region1::PCp2T(double p,double cp,int& itera){
 /*  由于过冷水相同压力下不同温度的Cp相差很小，PCp2T迭代出的t值很难还原回原来的t
 不建议采用P，Cp to求其他参数 */
    double pi=p/p_base;
 /*.利用拟合的粗略公式设置初值t*/
    double t=0.007783423198295928*p*p-39.310983833264864*cp*cp
            +2.9021945139958483*p*cp-12.068151636475518*p
            +492.3397038517826*cp-1269.9566380291237; 
    double tau=T_base/(t+T0);
    double lefts[34];
    for(int i = 0;i<34; i++)
		lefts[i]=ni[i]*pow(7.1-pi,Ii[i])*Ji[i]*(Ji[i]-1);    
    double gamma_tau2=0;
    for(int i = 0;i<34; i++)
        gamma_tau2+=lefts[i]*pow(tau-1.222,Ji[i]-2);    
    double cpcp=-tau*tau*gamma_tau2*R;
    itera=0;
    if(abs(cpcp-cp)>ERR2){
        double t0=t;
        double cp0=cpcp;
        t=IF97Region4::P2T(p)-15;//TODO:初值选择再优化
        t=(abs(t-t0)<15)?t0-15:t; 
        tau=T_base/(t+T0);
        gamma_tau2=0;
        for(int i = 0;i<34; i++)
            gamma_tau2+=lefts[i]*pow(tau-1.222,Ji[i]-2);
        cpcp=-tau*tau*gamma_tau2*R;   
        double t1=t;
        double cp1=cpcp;
        t=t1+(cp-cp1)/(cp0-cp1)*(t0-t1);
        tau=T_base/(t+T0);
        gamma_tau2=0;
        for(int i = 0;i<34; i++)
            gamma_tau2+=lefts[i]*pow(tau-1.222,Ji[i]-2);
        cpcp=-tau*tau*gamma_tau2*R;
        while(abs(cpcp-cp)>ERR2){
            itera++;
            t0=t1;
            cp0=cp1;
            t1=t;
            cp1=cpcp;
            t=t1+(cp-cp1)/(cp0-cp1)*(t0-t1);
            tau=T_base/(t+T0);
            gamma_tau2=0;
            for(int i = 0;i<34; i++)
                gamma_tau2+=lefts[i]*pow(tau-1.222,Ji[i]-2);
            cpcp=-tau*tau*gamma_tau2*R;
        }
    }
    return t;   
}
double IF97Region1::TH2P(double t,double h,int& itera){
    double err=ERR;
    double tau=T_base/(t+T0);
    double lefts[34];
    for(int i = 0;i<34; i++)
        lefts[i]=ni[i]*Ji[i]*pow(tau-1.222,Ji[i]-1);
    double p=10;
    double ps=IF97Region4::T2P(t);
    /*.利用拟合的粗略公式设置初值p*/
    if(t>=0 && t<=50){
        p=-0.00350725380*t*t+0.00045569756*h*h-0.00111040654*t*h
          -4.20897316924*t+1.01047891593*h-0.53520899497;
    }else if(t>50 && t<=100){
        p=-0.00934337710*t*t+0.00010799878*h*h+0.00172808463*t*h
          -4.31944050685*t+1.03981042874*h-1.27932307958;
    }else if(t>100 && t<=150){
        p=-0.02520474643*t*t-0.00049669320*h*h+0.00792833912*t*h
          -3.96834023880*t+0.96961365926*h-2.49893115932;
    }else if(t>150 && t<=200){
        p=-0.07911965905*t*t-0.00252768411*h*h+0.02897189024*t*h
          -1.74587371197*t+0.50532322819*h-14.31227705245;
    }else if(t>200 && t<=247.5711){
        p=-0.38598851231*t*t-0.01455590938*h*h+0.15088800230*t*h
          +14.29549174131*t-2.82596587386*h-156.75413188310;
    }else if(t>247.5711 && t<285.3041857){
        p=-2.1138006643*t*t-0.0929427611*h*h+0.8879602002*t*h
          +67.8651289239*t-14.7530585610*h-281.5752183792;
        double t_300=t/300;
        double p_fle=-7405.88491821*pow(t_300,6)+42720.33203809*pow(t_300,5)
              -102340.68781725*pow(t_300,4)+130341.47388109*pow(t_300,3)
              -93159.58960906*pow(t_300,2)+35800.45288382*t_300-5896.15285931;
        if(p<p_fle+0.001) p=p_fle+0.001;
    }else if(t>=285.3041857 && t<312.5){
        p=2.01660953383*t*t+0.08927982598*h*h-0.85041712566*t*h
          -65.07277755083*t+14.88446578182*h-351.12428675255;
        double t_300=t/300;
        double p_fle=-7405.88491821*pow(t_300,6)+42720.33203809*pow(t_300,5)
              -102340.68781725*pow(t_300,4)+130341.47388109*pow(t_300,3)
              -93159.58960906*pow(t_300,2)+35800.45288382*t_300-5896.15285931;
        if(p>p_fle-0.001) p=p_fle-0.001;
    }else if(t>=312.5 && t<339.8638){
        p=0.5949935025*t*t+0.0256339561*h*h-0.2484843299*t*h
          -16.6045383455*t+4.4809006333*h-695.9067672212;
        double t_300=t/300;
        double p_fle=-7405.88491821*pow(t_300,6)+42720.33203809*pow(t_300,5)
              -102340.68781725*pow(t_300,4)+130341.47388109*pow(t_300,3)
              -93159.58960906*pow(t_300,2)+35800.45288382*t_300-5896.15285931;
        if(p>p_fle-0.001) p=p_fle-0.001;
    }else if(t>=339.8638 && t<=350){
        p=0.29425198222*t*t+0.01287859709*h*h-0.12402004215*t*h
         -3.77056248947*t+1.45474007160*h-544.02491393411;
    }
    if(p<ps) p=ps;
    else if(p>100) p=100;
    double pi=p/p_base;   
    double gamma_tau=0;
    for(int i = 0;i<34; i++)
        gamma_tau+=lefts[i]*pow(7.1-pi,Ii[i]);
    double hh=tau*gamma_tau*(t+T0)*R;
    itera=0;           
    if(abs(hh-h)>err){
        double p0=p;
        double h0=hh;
        if(t>247.5711 && t<285.3041857){
            p=100;
            if(abs(p-p0)<0.2) p=p0-1;
        }else{
            p=ps;
            if(abs(p-p0)<0.2) p=p0+1;        
        }
        pi=p/p_base; 
        gamma_tau=0;
        for(int i = 0;i<34; i++)
            gamma_tau+=lefts[i]*pow(7.1-pi,Ii[i]);
        hh=tau*gamma_tau*(t+T0)*R;
        double p1=p;
        double h1=hh;
        p=p1+(h-h1)/(h0-h1)*(p0-p1);
        pi=p/p_base; 
        gamma_tau=0;
        for(int i = 0;i<34; i++)
            gamma_tau+=lefts[i]*pow(7.1-pi,Ii[i]);
        hh=tau*gamma_tau*(t+T0)*R;
        while(abs(hh-h)>err){
            itera++;
            p0=p1;
            h0=h1;
            p1=p;
            h1=hh;
            p=p1+(h-h1)/(h0-h1)*(p0-p1);
            pi=p/p_base; 
            gamma_tau=0;
            for(int i = 0;i<34; i++)
                gamma_tau+=lefts[i]*pow(7.1-pi,Ii[i]);
            hh=tau*gamma_tau*(t+T0)*R;
        } 
    }
    return p;     
}
/*获取从247.5712-339.8637度过冷水等温线焓值最小时的p
int count=0;
double PatMinH(double t,double left,double right,double left_h,double right_h){
    count++;
    if(right-left<=0.0000005){//求出的p误差范围是+-2*dp
        double dp=(right-left);
        double h0=IF97Region1::PT2H(left-dp,t);
        double h1=left_h;
        double h2=right_h;
        double h3=IF97Region1::PT2H(right+dp,t);
        if((h0>h1 || h0>h2) && (h3>h2 || h3>h1)){
            return (right+left)/2.0;
        }
        return -1;
    }
    double h=IF97Region1::PT2H((left+right)/2.0,t);
    if(left_h<h && h<right_h)
        return PatMinH(t,left,(left+right)/2.0,left_h,h);
    else if(left_h>h && h>right_h)
        return PatMinH(t,(left+right)/2.0,right,h,right_h);
    else{//h<left_h && h<right_h
        srand((int)time(0));
        int select=rand()%2;
        if(select==0){
            double p=PatMinH(t,left,(left+right)/2.0,left_h,h);
            if(p>0)
                return p;
            return PatMinH(t,(left+right)/2.0,right,h,right_h);
        }else{            
            double p=PatMinH(t,(left+right)/2.0,right,h,right_h);
            if(p>0)
                return p;
            return PatMinH(t,left,(left+right)/2.0,left_h,h);
        }
    }
}
void inflectionPs(){
    for(double t=339.8;t<339.9;t+=0.0001){
        count=0;
        double ps=IF97Region4::T2P(t);        
        double p=PatMinH(t,ps,100,IF97Region1::PT2H(ps,t),IF97Region1::PT2H(100,t));
        cout<<setprecision(10)<<t<<'\t'<<p<<'\t'<<count<<'\t';
        cout<<setprecision(13)<<IF97Region1::PT2H(p-0.000001,t)<<'\t';
        cout<<setprecision(13)<<IF97Region1::PT2H(p,t)<<'\t';
        cout<<setprecision(13)<<IF97Region1::PT2H(p+0.000001,t)<<endl;
    }
}
int main(){
    inflectionPs();
    return 0;
}
end 获取从247.6-339.9度过冷水等温线焓值最小时的p*/

void verify(){
    double p,t,h,s,tt,pp;
    int i=0,max=0;
    for(t=1;t<=247.5;t+=0.5){
        p=IF97Region4::T2P(t)+0.000000;
        for(;p<=100;p+=0.5){
            h=IF97Region1::PT2H(p,t);
            pp=IF97Region1::TH2P(t,h,i);
            if(i>max) max=i;
            cout<<setprecision(10)<<i<<'\t'<<h<<"\t"<<p<<"\t"<<pp<<"\t"<<t<<"\t"<<abs(100.0*(pp-p)/p)<<endl;
        }
    }
    double dt=(285.3041857-247.58)/200;  
    for(t=247.58;t<285.3041857;t+=dt){
        p=0.1-7405.88491821289*pow(t/300,6)+42720.3320380859*pow(t/300,5)-102340.687817246*pow(t/300,4)
            +130341.473881092*pow(t/300,3)-93159.5896090619*pow(t/300,2)
            +35800.4528838225*(t/300)-5896.15285930626;
        double p1=100;
        for(;p<=p1;p+=0.5){
            h=IF97Region1::PT2H(p,t);
            pp=IF97Region1::TH2P(t,h,i);
            if(i>max) max=i;
            cout<<setprecision(10)<<i<<'\t'<<h<<"\t"<<p<<"\t"<<pp<<"\t"<<t<<"\t"<<abs(100.0*(pp-p)/p)<<endl;
        }
    }
    dt=(339.8637-285.3041857)/200;    
    for(t=285.3041857;t<339.8638;t+=dt){
        p=IF97Region4::T2P(t)+0.000000;
        double p1=-0.1-7405.88491821289*pow(t/300,6)+42720.3320380859*pow(t/300,5)-102340.687817246*pow(t/300,4)
            +130341.473881092*pow(t/300,3)-93159.5896090619*pow(t/300,2)
            +35800.4528838225*(t/300)-5896.15285930626;
        for(;p<=p1;p+=0.5){
            h=IF97Region1::PT2H(p,t);
            pp=IF97Region1::TH2P(t,h,i);
            if(i>max) max=i;
            cout<<setprecision(10)<<i<<'\t'<<h<<"\t"<<p<<"\t"<<pp<<"\t"<<t<<"\t"<<abs(100.0*(pp-p)/p)<<endl;
        }
    }
    for(t=340;t<=350;t+=0.5){    
        p=IF97Region4::T2P(t)+0.000000;
        for(;p<=100;p+=0.5){
            h=IF97Region1::PT2H(p,t);
            pp=IF97Region1::TH2P(t,h,i);
            if(i>max) max=i;
            cout<<setprecision(10)<<i<<'\t'<<h<<"\t"<<p<<"\t"<<pp<<"\t"<<t<<"\t"<<abs(100.0*(pp-p)/p)<<endl;
        }
    } 
    cout<<max<<endl;
}
void createFitdata(){
    double p,t,h,s,tt,pp;
    //for(t=200;t<=247.5;t+=2.5){
    for(t=248;t<=285;t+=0.5){
        p=0.001-7405.88491821289*pow(t/300,6)+42720.3320380859*pow(t/300,5)-102340.687817246*pow(t/300,4)
            +130341.473881092*pow(t/300,3)-93159.5896090619*pow(t/300,2)
            +35800.4528838225*(t/300)-5896.15285930626;
        double p1=100;
        double dp=(p1-p)/60;
        while(p<=p1){
            h=IF97Region1::PT2H(p,t);
            cout<<setprecision(10)<<t*t<<"\t"<<h*h<<"\t"<<t*h<<"\t"
                <<t<<"\t"<<h<<"\t"<<1<<"\t"<<p<<endl;
            p+=dp;            
        }
    } 
}
void findMinEndPoint(){
    double t,p0,p1,h0,h1,h;
    double dt=(339.8637-285.3041857)/10000;
    for(t=285.3041857;t<=339.8637;t+=dt){
        t=339.8637;
        p0=IF97Region4::T2P(t);
        p1=-7405.88491821289*pow(t/300,6)+42720.3320380859*pow(t/300,5)-102340.687817246*pow(t/300,4)
            +130341.473881092*pow(t/300,3)-93159.5896090619*pow(t/300,2)
            +35800.4528838225*(t/300)-5896.15285930626;
        h=IF97Region1::PT2H(100,t);
        double p=(p0+p1)/2;
        double hh=IF97Region1::PT2H(p,t);
        while(p0<p1 && abs(hh-h)>1E-7){
            if(hh>h)
                p0=p;
            else
                p1=p;
            p=(p0+p1)/2.0;
            hh=IF97Region1::PT2H(p,t);
        }
        cout<<setprecision(15)<<t<<"\t"<<p<<"\t"<<hh<<"\t"<<abs(100*(hh-h)/h)<<endl;
    }
}
int main(){ 
    double p,t,h;
    int i;
    // t=290.214542;
    // h=1290.954507;
    // p=IF97Region1::TH2P(t,h,i);
    //findMinEndPoint();
    verify();
    //createFitdata();
    return 0;
} 

/*int main(){
    double p,t,h,s,cp;
    int i=0;
    p=3;
    t=300-273.15;
    s=IF97Region1::PT2S(p,t);
    h=IF97Region1::PS2H(p,s);
    cout<<setprecision(10)<<t<<'\t'<<s<<"\t"<<p<<"\t"<<h<<"\t"<<endl;
     for(t=5;t<=350;t+=1){
        p=IF97Region4::T2P(t);
        for(;p<=100;p+=1){
            s=IF97Region1::PT2S(p,t);
            double tt=IF97Region1::PS2T(p,s,i);
            //cout<<setprecision(10)<<i<<"\t"<<tt<<"\t"<<t<<"\t"<<100.0*(t-tt)/t<<endl;
        }
    } 
    for(t=5;t<350.1;t+=5){
        p=IF97Region4::T2P(t)+0.00001;
        double dp=(100-p)/25;
        while(p<100.1){
            cp=IF97Region1::PT2Cp(p,t);
            cout<<setprecision(10)<<p*p<<"\t"<<cp*cp<<"\t"<<p*cp<<"\t"
                <<p<<"\t"<<cp<<"\t"<<1<<"\t"<<t<<endl;
            
        }
    }
    return 0; 
}
*/
/* 
int main(){
    double p,t,h,s,cp;
    int i=0;
    p=3;
    t=300-273.15;
    cp=IF97Region1::PT2Cp(p,t);
    double tt=IF97Region1::PCp2T(p,cp,i);
    cout<<setprecision(10)<<i<<'\t'<<cp<<"\t"<<p<<"\t"<<tt<<"\t"<<100.0*(t-tt)/t<<endl;
    for(t=5;t<=350;t+=1){
        p=IF97Region4::T2P(t);
        for(;p<=100;p+=1){
            cp=IF97Region1::PT2Cp(p,t);
            double tt=IF97Region1::PCp2T(p,cp,i);
            cout<<setprecision(10)<<i<<"\t"<<p<<"\t"<<t<<"\t"<<100.0*(t-tt)/t<<endl;
        }
    }
    for(t=5;t<350.1;t+=5){
        p=IF97Region4::T2P(t)+0.00001;
        double dp=(100-p)/25;
        while(p<100.1){
            cp=IF97Region1::PT2Cp(p,t);
            cout<<setprecision(10)<<p*p<<"\t"<<cp*cp<<"\t"<<p*cp<<"\t"
                <<p<<"\t"<<cp<<"\t"<<1<<"\t"<<t<<endl;
            
        }
    }
    return 0;
} */
    
/* int main(){
    double p,t,h,s,u;
    p=3;
    t=500-273.15;
    u=IF97Region1::PT2U(p,t);
    double tt=IF97Region1::PU2T(p,u);
    cout<<setprecision(10)<<u<<"\t"<<p<<"\t"<<tt<<"\t"<<100.0*(t-tt)/t<<endl;

    int i=0;
    for(t=5;t<=350;t+=1){
        p=IF97Region4::T2P(t);
        for(;p<=100;p+=1){
            u=IF97Region1::PT2U(p,t);
            //double tt=IF97Region1::PV2T(p,v);
            double tt=IF97Region1::PU2T(p,u,i);
            cout<<setprecision(10)<<i<<"\t"<<p<<"\t"<<t<<"\t"<<100.0*(t-tt)/t<<endl;
        }
    } 
}

int main(){
    double p,t,h,s,v;
    int i=0;
    for(t=5;t<=350;t+=1){
        p=IF97Region4::T2P(t);
        for(;p<=100;p+=1){
            v=IF97Region1::PT2V(p,t);
            //double tt=IF97Region1::PV2T(p,v);
            double tt=IF97Region1::PV2T(p,v,i);
            cout<<setprecision(10)<<i<<"\t"<<p<<"\t"<<t<<"\t"<<100.0*(t-tt)/t<<endl;
        }
    }
    p=20.54014805;
    t=345;
    v=IF97Region1::PT2V(p,t);
    double tt=IF97Region1::PV2T(p,v,i);
    cout<<setprecision(10)<<i<<"\t"<<p<<"\t"<<tt<<"\t"<<100.0*(t-tt)/t<<endl;
}
int main(){
    double p,t,h,s,u,v,cp,cv;
    cout<<"p\tt\th\ts\tu\tv\tcp\tcv"<<endl;
    for(t=5;t<350.1;t+=5){
        p=IF97Region4::T2P(t)+0.00001;
        double dp=(100-p)/25;
        while(p<100.1){
            v=1000*IF97Region1::PT2V(p,t);
            h=IF97Region1::PT2H(p,t);
            s=IF97Region1::PT2S(p,t);
            u=IF97Region1::PT2U(p,t);
            cp=IF97Region1::PT2Cp(p,t);
            cv=IF97Region1::PT2Cv(p,t);
            cout<<setprecision(10)<<p<<"\t"<<t<<"\t"<<h<<"\t"<<s<<"\t"<<u<<"\t"<<v
                <<"\t"<<cp<<"\t"<<cv<<endl;
            //cout<<setprecision(10)<<p<<"\t"<<v<<"\t"<<t<<endl;
            p+=dp;
        }
    }
}
  */