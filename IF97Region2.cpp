#include "IF97Region2.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <time.h>
using namespace std;

double IF97Region2::PT2H(double p,double t){
    double pi=p/1.0;
    double tau=540.0/(t+T0);
    double gamma0_tau=0;
    for(int i=0;i<9;i++){
        gamma0_tau+=ni_tab10[i]*Ji_tab10[i]*pow(tau,Ji_tab10[i]-1);
    }
    double gammar_tau=0;
    for(int i=0;i<43;i++){
        gammar_tau+=ni_tab11[i]*pow(pi,Ii_tab11[i])*Ji_tab11[i]*pow(tau-0.5,Ji_tab11[i]-1);
    }
    double h=tau*(gamma0_tau+gammar_tau)*(t+T0)*R;
    return h;
}
double IF97Region2::PT2S(double p,double t){
    double pi=p/1.0;
    double tau=540.0/(t+T0);
    double gamma0=log(pi);
    double gamma0_tau=0;
    for(int i=0;i<9;i++){
        gamma0+=ni_tab10[i]*pow(tau,Ji_tab10[i]);
        gamma0_tau+=ni_tab10[i]*Ji_tab10[i]*pow(tau,Ji_tab10[i]-1);
    }
    double gammar=0;
    double gammar_tau=0;
    for(int i=0;i<43;i++){
        gammar+=ni_tab11[i]*pow(pi,Ii_tab11[i])*pow(tau-0.5,Ji_tab11[i]);
        gammar_tau+=ni_tab11[i]*pow(pi,Ii_tab11[i])*Ji_tab11[i]*pow(tau-0.5,Ji_tab11[i]-1);
    }
    double s=(tau*(gamma0_tau+gammar_tau)-(gamma0+gammar))*R;
    return s;
}
double IF97Region2::PT2V(double p,double t){    
    double pi=p/1.0;
    double tau=540.0/(t+T0);
    double gamma0_pi=1/pi;
    double gammar_pi=0;
    for(int i=0;i<43;i++) {
        gammar_pi+=ni_tab11[i]*Ii_tab11[i]*pow(pi,Ii_tab11[i]-1)*pow(tau-0.5,Ji_tab11[i]);
    }
	double v=(pi*(gamma0_pi+gammar_pi))*(t+T0)*R/(p*1000);
    return v;
}
double IF97Region2::PT2U(double p,double t){   
    double pi=p/1.0;
    double tau=540.0/(t+T0);
    double gamma0_pi=1/pi;
    double gamma0_tau=0;
    for(int i=0;i<9;i++){
        gamma0_tau+=ni_tab10[i]*Ji_tab10[i]*pow(tau,Ji_tab10[i]-1);
    }
    double gammar_pi=0;    
    double gammar_tau=0;
    for(int i=0;i<43;i++){
        gammar_pi+=ni_tab11[i]*Ii_tab11[i]*pow(pi,Ii_tab11[i]-1)*pow(tau-0.5,Ji_tab11[i]);
        gammar_tau+=ni_tab11[i]*pow(pi,Ii_tab11[i])*Ji_tab11[i]*pow(tau-0.5,Ji_tab11[i]-1);
    }
	double u=(tau*(gamma0_tau+gammar_tau)-pi*(gamma0_pi+gammar_pi))*(t+T0)*R;
    return u;
}
double IF97Region2::PT2Cp(double p,double t){ 
    double pi=p/1.0;
    double tau=540.0/(t+T0);
    double gamma0_tau2=0;
    for(int i=0;i<9;i++){
        gamma0_tau2+=ni_tab10[i]*Ji_tab10[i]*(Ji_tab10[i]-1)*pow(tau,Ji_tab10[i]-2);
    }
    double gammar_tau2=0;
    for(int i=0;i<43;i++){
        gammar_tau2+=ni_tab11[i]*pow(pi,Ii_tab11[i])*Ji_tab11[i]*(Ji_tab11[i]-1)
                    *pow(tau-0.5,Ji_tab11[i]-2);
    }
	double cp=(-tau*tau*(gamma0_tau2+gammar_tau2))*R;  
    return cp;
} 
double IF97Region2::PT2Cv(double p,double t){ 
    double pi=p/1.0;
    double tau=540.0/(t+T0);
    double gamma0_tau2=0;
    for(int i=0;i<9;i++){
        gamma0_tau2+=ni_tab10[i]*Ji_tab10[i]*(Ji_tab10[i]-1)*pow(tau,Ji_tab10[i]-2);
    }
    double gammar_tau2=0;
	double gammar_pi=0;
	double gammar_pi2=0;
	double gammar_pitau=0;    
    for(int i=0;i<43;i++){
        gammar_tau2+=ni_tab11[i]*pow(pi,Ii_tab11[i])*Ji_tab11[i]*(Ji_tab11[i]-1)
                    *pow(tau-0.5,Ji_tab11[i]-2);
        gammar_pi+=ni_tab11[i]*Ii_tab11[i]*pow(pi,Ii_tab11[i]-1)*pow(tau-0.5,Ji_tab11[i]);
		gammar_pi2+=ni_tab11[i]*Ii_tab11[i]*(Ii_tab11[i]-1)*pow(pi,Ii_tab11[i]-2)
                    *pow(tau-0.5,Ji_tab11[i]);
		gammar_pitau+=ni_tab11[i]*Ii_tab11[i]*pow(pi,Ii_tab11[i]-1)
                    *Ji_tab11[i]*pow(tau-0.5,Ji_tab11[i]-1);
    }
	double cv=(-tau*tau*(gamma0_tau2+gammar_tau2)
            -pow(1.0+pi*gammar_pi-tau*pi*gammar_pitau,2)/(1.0-pi*pi*gammar_pi2))*R;
	return cv;
}
double IF97Region2::PT2W(double p,double t){
    return 0.0;
}
double IF97Region2::PH2T(double p,double h,int& itera){
    double err=ERR;
    double pi=p/1.0;
    double left0[9];
    for(int i=0;i<9;i++){
        left0[i]=ni_tab10[i]*Ji_tab10[i];
    }
    double leftr[43];    
    for(int i=0;i<43;i++){
        leftr[i]=ni_tab11[i]*pow(pi,Ii_tab11[i])*Ji_tab11[i];//*pow(tau-0.5,Ji_tab11[i]-1);
    }
    double t;
    if(p<=4){
        t=PH2T2a(p,h);
    }else if(h>=H2bc(p)){
        t=PH2T2b(p,h);
    }else{
        t=PH2T2c(p,h);
    }
    double tau=540.0/(t+T0);
    double gamma0_tau=0;
    for(int i=0;i<9;i++){
        gamma0_tau+=left0[i]*pow(tau,Ji_tab10[i]-1);
    }
    double gammar_tau=0;    
    for(int i=0;i<43;i++){
        gammar_tau+=leftr[i]*pow(tau-0.5,Ji_tab11[i]-1);
    }
    double hh=tau*(gamma0_tau+gammar_tau)*(t+T0)*R;
    itera=0;
    if(abs(hh-h)>err){
        double t0=t;
        double h0=hh;
        double t1=t+0.1>800.0?t-0.1:t+0.1;
        tau=540.0/(t1+T0);
        gamma0_tau=0;
        for(int i=0;i<9;i++){
            gamma0_tau+=left0[i]*pow(tau,Ji_tab10[i]-1);
        }
        gammar_tau=0;    
        for(int i=0;i<43;i++){
            gammar_tau+=leftr[i]*pow(tau-0.5,Ji_tab11[i]-1);
        }
        double h1=tau*(gamma0_tau+gammar_tau)*(t1+T0)*R;
        t=t1+(h-h1)/(h0-h1)*(t0-t1);
        tau=540.0/(t+T0);
        gamma0_tau=0;
        for(int i=0;i<9;i++){
            gamma0_tau+=left0[i]*pow(tau,Ji_tab10[i]-1);
        }
        gammar_tau=0;    
        for(int i=0;i<43;i++){
            gammar_tau+=leftr[i]*pow(tau-0.5,Ji_tab11[i]-1);
        }
        hh=tau*(gamma0_tau+gammar_tau)*(t+T0)*R;
        while(abs(hh-h)>err){
            itera++;
            t0=t1;
            h0=h1;
            t1=t;
            h1=hh;
            t=t1+(h-h1)/(h0-h1)*(t0-t1);
            tau=540.0/(t+T0);
            gamma0_tau=0;
            for(int i=0;i<9;i++){
                gamma0_tau+=left0[i]*pow(tau,Ji_tab10[i]-1);
            }
            gammar_tau=0;    
            for(int i=0;i<43;i++){
                gammar_tau+=leftr[i]*pow(tau-0.5,Ji_tab11[i]-1);
            }
            hh=tau*(gamma0_tau+gammar_tau)*(t+T0)*R;
        }
    }
    return t;
}
double IF97Region2::PS2T(double p,double s,int& itera){
    double err=ERR2;
    double pi=p/1.0;
    double lnpi=log(pi);  
    double left[43];
    double left_tau[43];    
    for(int i=0;i<43;i++){
        left[i]=ni_tab11[i]*pow(pi,Ii_tab11[i]);//*pow(tau-0.5,Ji_tab11[i]);
        left_tau[i]=ni_tab11[i]*pow(pi,Ii_tab11[i])*Ji_tab11[i];//*pow(tau-0.5,Ji_tab11[i]-1);
    }
    double t;
    if(p<=4){
        t=PS2T2a(p,s);
    }else if(s>=S2bc){
        t=PS2T2b(p,s);
    }else{
        t=PS2T2c(p,s);
    }
    double tau=540.0/(t+T0);
    double gamma0=lnpi;
    double gamma0_tau=0;    
    for(int i=0;i<9;i++){
        gamma0+=ni_tab10[i]*pow(tau,Ji_tab10[i]);
        gamma0_tau+=ni_tab10[i]*Ji_tab10[i]*pow(tau,Ji_tab10[i]-1);
    }
    double gammar=0;
    double gammar_tau=0;       
    for(int i=0;i<43;i++){
        gammar+=left[i]*pow(tau-0.5,Ji_tab11[i]);
        gammar_tau+=left_tau[i]*pow(tau-0.5,Ji_tab11[i]-1);
    }
    double ss=(tau*(gamma0_tau+gammar_tau)-(gamma0+gammar))*R;
    itera=0;
    if(abs(ss-s)>err){
        double t0=t;
        double s0=ss;
        double t1=t+0.1>800.0?t-0.1:t+0.1;
        tau=540.0/(t1+T0);
        gamma0=lnpi;
        gamma0_tau=0;    
        for(int i=0;i<9;i++){
            gamma0+=ni_tab10[i]*pow(tau,Ji_tab10[i]);
            gamma0_tau+=ni_tab10[i]*Ji_tab10[i]*pow(tau,Ji_tab10[i]-1);
        }
        gammar=0;
        gammar_tau=0;       
        for(int i=0;i<43;i++){
            gammar+=left[i]*pow(tau-0.5,Ji_tab11[i]);
            gammar_tau+=left_tau[i]*pow(tau-0.5,Ji_tab11[i]-1);
        }
        double s1=(tau*(gamma0_tau+gammar_tau)-(gamma0+gammar))*R;
        t=t1+(s-s1)/(s0-s1)*(t0-t1);
        tau=540.0/(t+T0);
        gamma0=lnpi;
        gamma0_tau=0;    
        for(int i=0;i<9;i++){
            gamma0+=ni_tab10[i]*pow(tau,Ji_tab10[i]);
            gamma0_tau+=ni_tab10[i]*Ji_tab10[i]*pow(tau,Ji_tab10[i]-1);
        }
        gammar=0;
        gammar_tau=0;       
        for(int i=0;i<43;i++){
            gammar+=left[i]*pow(tau-0.5,Ji_tab11[i]);
            gammar_tau+=left_tau[i]*pow(tau-0.5,Ji_tab11[i]-1);
        }
        ss=(tau*(gamma0_tau+gammar_tau)-(gamma0+gammar))*R;
        while(abs(ss-s)>err){
            itera++;
            t0=t1;
            s0=s1;
            t1=t;
            s1=ss;
            t=t1+(s-s1)/(s0-s1)*(t0-t1);
            tau=540.0/(t+T0);
            gamma0=lnpi;
            gamma0_tau=0;    
            for(int i=0;i<9;i++){
                gamma0+=ni_tab10[i]*pow(tau,Ji_tab10[i]);
                gamma0_tau+=ni_tab10[i]*Ji_tab10[i]*pow(tau,Ji_tab10[i]-1);
            }
            gammar=0;
            gammar_tau=0;       
            for(int i=0;i<43;i++){
                gammar+=left[i]*pow(tau-0.5,Ji_tab11[i]);
                gammar_tau+=left_tau[i]*pow(tau-0.5,Ji_tab11[i]-1);
            }
            ss=(tau*(gamma0_tau+gammar_tau)-(gamma0+gammar))*R;
        }
    }
    return t;
}
/*P,H or P,S 辅助函数 */
double IF97Region2::P2bc(double h){
    double eta=h/1.0;
    double pi=0.90584278514723E3-0.67955786399241*eta+0.12809002730136E-3*eta*eta;
    return 1.0*pi;
}
double IF97Region2::H2bc(double p){
    double pi=p/1.0;
    double eta=0.26526571908428E4+pow((pi-0.45257578905948E1)/0.12809002730136E-3,0.5);
    return 1.0*eta;
}
double IF97Region2::PH2T2a(double p,double h){
    double ni[34]={//if97 table20 page22
			 0.10898952318288E4,
			 0.84951654495535E3,
			-0.10781748091826E3,
			 0.33153654801263E2,
			-0.74232016790248E1,
			 0.11765048724356E2,
			 0.18445749355790E1,
			-0.41792700549624E1,
			 0.62478196935812E1,
			-0.17344563108114E2,
			-0.20058176862096E3,
			 0.27196065473796E3,
			-0.45511318285818E3,
			 0.30919688604755E4,
			 0.25226640357872E6,
			-0.61707422868339E-2,
			-0.31078046629583E0,
			 0.11670873077107E2,
			 0.12812798404046E9,
			-0.98554909623276E9,
			 0.28224546973002E10,
			-0.35948971410703E10,
			 0.17227349913197E10,
			-0.13551334240775E5,
			 0.12848734664650E8,
			 0.13865724283226E1,
			 0.23598832556514E6,
			-0.13105236545054E8,
			 0.73999835474766E4,
			-0.55196697030060E6,
			 0.37154085996233E7,
			 0.19127729239660E5,
			-0.41535164835634E6,
			-0.62459855192507E2
    };
	int Ii[34]={
			0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,4,4,4,5,5,5,6,6,7
    };
	int	Ji[34]={
		0,1,2,3,7,20,0,1,2,3,7,9,11,18,44,0,2,7,36,38,
        40,42,44,24,44,12,32,44,32,36,42,34,44,28
	};
    double pi=p/1.0;
    double eta=h/2000.0;
    double theta=0;
    for(int i=0;i<34;i++){
        theta+=ni[i]*pow(pi,Ii[i])*pow(eta-2.1,Ji[i]);
    }
    return 1.0*theta-T0;   
}
double IF97Region2::PH2T2b(double p,double h){
    double ni[38]={//if97 table21 page23
			 0.14895041079516E4,
			 0.74307798314034E3,
			-0.97708318797837E2,
			 0.24742464705674E1,
			-0.63281320016026E0,
			 0.11385952129658E1,
			-0.47811863648625E0,
			 0.85208123431544E-2,
			 0.93747147377932E0,
			 0.33593118604916E1,
			 0.33809355601454E1,
			 0.16844539671904E0,
			 0.73875745236695E0,
			-0.47128737436186E0,
			 0.15020273139707E0,
			-0.21764114219750E-2,
			-0.21810755324761E-1,
			-0.10829784403677E0,
			-0.46333324635812E-1,
			 0.71280351959551E-4,
			 0.11032831789999E-3,
			 0.18955248387902E-3,
			 0.30891541160537E-2,
			 0.13555504554949E-2,
			 0.28640237477456E-6,
			-0.10779857357512E-4,
			-0.76462712454814E-4,
			 0.14052392818316E-4,
			-0.31083814331434E-4,
			-0.10302738212103E-5,
			 0.28217281635040E-6,
			 0.12704902271945E-5,
			 0.73803353468292E-7,
			-0.11030139238909E-7,
			-0.81456365207833E-13,
			-0.25180545682962E-10,
			-0.17565233969407E-17,
			 0.86934156344163E-14
    };
	int Ii[38]={
		0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,
        2,3,3,3,3,4,4,4,4,4,4,5,5,5,6,7,7,9,9
    };
	int	Ji[38]={
        0,1,2,12,18,24,28,40,0,2,6,12,18,24,28,40,2,8,
        18,40,1,2,12,24,2,12,18,24,28,40,18,24,40,28,2,28,1,40
	};
    double pi=p/1.0;
    double eta=h/2000.0;
    double theta=0;
    for(int i=0;i<38;i++){
        theta+=ni[i]*pow(pi-2.0,Ii[i])*pow(eta-2.6,Ji[i]);
    }
    return 1.0*theta-T0;   
}
double IF97Region2::PH2T2c(double p,double h){
    double ni[23]={//if97 table22 page24
			-0.32368398555242E13,
			 0.73263350902181E13,
			 0.35825089945447E12,
			-0.58340131851590E12,
			-0.10783068217470E11,
			 0.20825544563171E11,
			 0.61074783564516E6,
			 0.85977722535580E6,
			-0.25745723604170E5,
			 0.31081088422714E5,
			 0.12082315865936E4,
			 0.48219755109255E3,
			 0.37966001272486E1,
			-0.10842984880077E2,
			-0.45364172676660E-1,
			 0.14559115658698E-12,
			 0.11261597407230E-11,
			-0.17804982240686E-10,
			 0.12324579690832E-6,
			-0.11606921130984E-5,
			 0.27846367088554E-4,
			-0.59270038474176E-3,
			 0.12918582991878E-2
    };
	int Ii[23]={
		-7,-7,-6,-6,-5,-5,-2,-2,-1,-1,0,0,1,1,2,6,6,6,6,6,6,6,6
    };
	int	Ji[23]={
        0,4,0,2,0,2,0,1,0,2,0,1,4,8,4,0,1,4,10,12,16,20,22
	};
    double pi=p/1.0;
    double eta=h/2000.0;
    double theta=0;
    for(int i=0;i<23;i++){
        theta+=ni[i]*pow(pi+25.0,Ii[i])*pow(eta-1.8,Ji[i]);
    }
    return 1.0*theta-T0;   
}
double IF97Region2::PS2T2a(double p,double s){
	double ni[46]={//if97 table25 page26
			-0.39235983861984E6,
			 0.51526573827270E6,
			 0.40482443161048E5,
			-0.32193790923902E3,
			 0.96961424218694E2,
			-0.22867846371773E2,
			-0.44942914124357E6,
			-0.50118336020166E4,
			 0.35684463560015E0,
			 0.44235335848190E5,
			-0.13673388811708E5,
			 0.42163260207864E6,
			 0.22516925837475E5,
			 0.47442144865646E3,
			-0.14931130797647E3,
			-0.19781126320452E6,
			-0.23554399470760E5,
			-0.19070616302076E5,
			 0.55375669883164E5,
			 0.38293691437363E4,
			-0.60391860580567E3,
			 0.19363102620331E4,
			 0.42660643698610E4,
			-0.59780638872718E4,
			-0.70401463926862E3,
			 0.33836784107553E3,
			 0.20862786635187E2,
			 0.33834172656196E-1,
			-0.43124428414893E-4,
			 0.16653791356412E3,
			-0.13986292055898E3,
			-0.78849547999872E0,
			 0.72132411753872E-1,
			-0.59754839398283E-2,
			-0.12141358953904E-4,
			 0.23227096733871E-6,
			-0.10538463566194E2,
			 0.20718925496502E1,
			-0.72193155260427E-1,
			 0.20749887081120E-6,
			-0.18340657911379E-1,
			 0.29036272348696E-6,
			 0.21037527893619E0,
			 0.25681239729999E-3,
			-0.12799002933781E-1,
			-0.82198102652018E-5
    };
	double Ii[46]={
        -1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.25,-1.25,-1.25,-1,-1,-1,-1,-1,-1,
        -0.75,-0.75,-0.5,-0.5,-0.5,-0.5,-0.25,-0.25,-0.25,-0.25,0.25,0.25,
        0.25,0.25,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.75,0.75,0.75,0.75,1,1,1.25,1.25,1.5,1.5
    };
    int	Ji[46]={
        -24,-23,-19,-13,-11,-10,-19,-15,-6,-26,-21,-17,-16,-9,-8,-15,-14,-26,-13,
        -9,-7,-27,-25,-11,-6,1,4,8,11,0,1,5,6,10,14,16,0,4,9,17,7,18,3,15,5,18
	};
    double pi=p/1.0;
    double sigma=s/2.0;
    double theta=0;
    for(int i=0;i<46;i++){
        theta+=ni[i]*pow(pi,Ii[i])*pow(sigma-2.0,Ji[i]);
    } 
    return 1*theta-T0;
}
double IF97Region2::PS2T2b(double p,double s){
	double ni[44]={ //if97 table26 page27
			 0.31687665083497E6,
			 0.20864175881858E2,
			-0.39859399803599E6,
			-0.21816058518877E2,
			 0.22369785194242E6,
			-0.27841703445817E4,
			 0.99207436071480E1,
			-0.75197512299157E5,
			 0.29708605951158E4,
			-0.34406878548526E1,
			 0.38815564249115E0,
			 0.17511295085750E5,
			-0.14237112854449E4,
			 0.10943803364167E1,
			 0.89971619308495E0,
			-0.33759740098958E4,
			 0.47162885818355E3,
			-0.19188241993679E1,
			 0.41078580492196E0,
			-0.33465378172097E0,
			 0.13870034777505E4,
			-0.40663326195838E3,
			 0.41727347159610E2,
			 0.21932549434532E1,
			-0.10320050009077E1,
			 0.35882943516703E0,
			 0.52511453726066E-2,
			 0.12838916450705E2,
			-0.28642437219381E1,
			 0.56912683664855E0,
			-0.99962954584931E-1,
			-0.32632037778459E-2,
			 0.23320922576723E-3,
			-0.15334809857450E0,
			 0.29072288239902E-1,
			 0.37534702741167E-3,
			 0.17296691702411E-2,
			-0.38556050844504E-3,
			-0.35017712292608E-4,
			-0.14566393631492E-4,
			 0.56420857267269E-5,
			 0.41286150074605E-7,
			-0.20684671118824E-7,
			 0.16409393674725E-8
    };
	int Ii[44]={
		-6,-6,-5,-5,-4,-4,-4,-3,-3,-3,-3,-2,-2,-2,-2,-1,-1,-1,
        -1,-1,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,3,3,3,4,4,5,5,5
    };
	int	Ji[44]={
		0,11,0,11,0,1,11,0,1,11,12,0,1,6,10,0,1,5,8,9,
        0,1,2,4,5,6,9,0,1,2,3,7,8,0,1,5,0,1,3,0,1,0,1,2
    };
    double pi=p/1.0;
    double sigma=s/0.7853;
    double theta=0;
    for(int i=0;i<44;i++){
        theta+=ni[i]*pow(pi,Ii[i])*pow(10.0-sigma,Ji[i]);
    } 
    return 1*theta-T0;
}
double IF97Region2::PS2T2c(double p,double s){
	double ni[30]={//if97 table27 page28
			 0.90968501005365E3,
			 0.24045667088420E4,
			-0.59162326387130E3,
			 0.54145404128074E3,
			-0.27098308411192E3,
			 0.97976525097926E3,
			-0.46966772959435E3,
			 0.14399274604723E2,
			-0.19104204230429E2,
			 0.53299167111971E1,
			-0.21252975375934E2,
			-0.31147334413760E0,
			 0.60334840894623E0,
			-0.42764839702509E-1,
			 0.58185597255259E-2,
			-0.14597008284753E-1,
			 0.56631175631027E-2,
			-0.76155864584577E-4,
			 0.22440342919332E-3,
			-0.12561095013413E-4,
			 0.63323132660934E-6,
			-0.20541989675375E-5,
			 0.36405370390082E-7,
			-0.29759897789215E-8,
			 0.10136618529763E-7,
			 0.59925719692351E-11,
			-0.20677870105164E-10,
			-0.20874278181886E-10,
			 0.10162166825089E-9,
			-0.16429828281347E-9
    };
	int Ii[30]={
        -2,-2,-1,0,0,0,0,1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,7,7,7,7,7
    };
    int	Ji[30]={
        0,1,0,0,1,2,3,0,1,3,4,0,1,2,0,1,5,0,1,4,0,1,2,0,1,0,1,3,4,5
	};
    double pi=p/1.0;
    double sigma=s/2.9251;
    double theta=0;
    for(int i=0;i<30;i++){
        theta+=ni[i]*pow(pi,Ii[i])*pow(2.0-sigma,Ji[i]);
    } 
    return 1*theta-T0;
}
int main(){ 
    double p,t,h,u,s,v,cp,cv,pp,tt,p2;
    int i;
    p=30;
    t=700-273.15;
    // h=IF97Region2::PT2H(p,t);
    // s=IF97Region2::PT2S(p,t);
    // v=IF97Region2::PT2V(p,t);
    // u=IF97Region2::PT2U(p,t);
    // cp=IF97Region2::PT2Cp(p,t);
    // cv=IF97Region2::PT2Cv(p,t);
    p=0.0035;
    s=8.52238967;  
    t=IF97Region2::PS2T(p,s,i);
    p=0.0035;
    s=10.1749996;  
    t=IF97Region2::PS2T(p,s,i);
    p=30;
    s=5.17540298;  
    t=IF97Region2::PS2T(p,s,i);     
    return 0;
} 