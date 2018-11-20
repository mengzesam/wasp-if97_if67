/*
    reference book:IAPWS R7-97(2012)
*/
#include "IF97Base.h"
class IF97Region3:public IF97Base{
private://构造函数为private，禁止实例化IF97Region3
    IF97Region3(){};
public: //static function: T,V to 
    static double TV2P(double t,double v);
    static double TV2H(double t,double v);
    static double TV2S(double t,double v);
    static double TV2U(double t,double v);
    static double TV2Cp(double t,double v);
    static double TV2Cv(double t,double v);
    static double TV2W(double t,double v);
public://static function:P,T to
    static void testPT();
    static double PT2V(double p,double t, int& itera);
public://static function:P,H to
    static double PH2T(double p,double h, int& itera);
public://static function:P,S to
    static double PS2T(double p,double s, int& itera);
private:////static function:P,H or P,S to 辅助函数，根据Supp-Tv(ph,ps)3-2014.pdf
    const static double S3ab; //==Sc
    static double H3ab(double p);
    static double PH2T3a(double p,double h);
    static double PH2T3b(double p,double h);
    static double PH2V3a(double p,double h);
    static double PH2V3b(double p,double h);
    static double PS2T3a(double p,double s);
    static double PS2T3b(double p,double s);
    static double PS2V3a(double p,double s);
    static double PS2V3b(double p,double s);
private://static function:P,T to 辅助函数，根据Supp-VPT3-2016.pdf
    static double T3ab(double p);
    static double T3cd(double p);
    static double T3ef(double p);
    static double T3gh(double p);
    static double T3ij(double p);
    static double T3jk(double p);
    static double T3mn(double p);
    static double T3op(double p);
    static double T3qu(double p);
    static double T3rx(double p);
    static double T3uv(double p);
    static double T3wx(double p);
    static double PT2V3a(double p,double t);
    static double PT2V3b(double p,double t);
    static double PT2V3c(double p,double t);
    static double PT2V3d(double p,double t);
    static double PT2V3e(double p,double t);
    static double PT2V3f(double p,double t);
    static double PT2V3g(double p,double t);
    static double PT2V3h(double p,double t);
    static double PT2V3i(double p,double t);
    static double PT2V3j(double p,double t);
    static double PT2V3k(double p,double t);
    static double PT2V3l(double p,double t);
    static double PT2V3m(double p,double t);
    static double PT2V3n(double p,double t);
    static double PT2V3o(double p,double t);
    static double PT2V3p(double p,double t);
    static double PT2V3q(double p,double t);
    static double PT2V3r(double p,double t);
    static double PT2V3s(double p,double t);
    static double PT2V3t(double p,double t);
    static double PT2V3u(double p,double t);
    static double PT2V3v(double p,double t);
    static double PT2V3w(double p,double t);
    static double PT2V3x(double p,double t);
    static double PT2V3y(double p,double t);
    static double PT2V3z(double p,double t);



    static void HS2PT(double h,double s,double& p,double& t,int& itera);
    static void HS2PT(double h,double s,double& p,double& t);
private:
    const static double ni[40];
    const static int Ii[39];
    const static int Ji[39];
};
const double IF97Region3::ni[40]={ //page30 table30
			-0.15732845290239E2,
			 0.20944396974307E2,
			-0.76867707878716E1,
			 0.26185947787954E1,
			-0.28080781148620E1,
			 0.12053369696517E1,
			-0.84566812812502E-2,
			-0.12654315477714E1,
			-0.11524407806681E1,
			 0.88521043984318E0,
			-0.64207765181607E0,
			 0.38493460186671E0,
			-0.85214708824206E0,
			 0.48972281541877E1,
			-0.30502617256965E1,
			 0.39420536879154E-1,
			 0.12558408424308E0,
			-0.27999329698710E0,
			 0.13899799569460E1,
			-0.20189915023570E1,
			-0.82147637173963E-2,
			-0.47596035734923E0,
			 0.43984074473500E-1,
			-0.44476435428739E0,
			 0.90572070719733E0,
			 0.70522450087967E0,
			 0.10770512626332E0,
			-0.32913623258954E0,
			-0.50871062041158E0,
			-0.22175400873096E-1,
			 0.94260751665092E-1,
			 0.16436278447961E0,
			-0.13503372241348E-1,
			-0.14834345352472E-1,
			 0.57922953628084E-3,
			 0.32308904703711E-2,
			 0.80964802996215E-4,
			-0.16557679795037E-3,
			-0.44923899061815E-4,
			 0.10658070028513E1
 };
const int IF97Region3::Ii[39]={//page30 table30
    0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,6,6,6,7,8,9,9,10,10,11
 };
const int IF97Region3::Ji[39]={//page30 table30
    0,1,2,7,10,12,23,2,6,15,17,0,2,6,7,22,26,0,2,4,16,26,0,2,4,26,1,3,26,0,2,26,2,26,2,26,0,1,26
 };
const double S3ab=4.41202148223476;//==Sc
