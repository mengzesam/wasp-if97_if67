/*
    reference book:IAPWS R7-97(2012)
*/
class IF97Region1{
private://构造函数为private，禁止实例化IF97Region1
    IF97Region1(){};
public: //static function: P,T to 
    static double PT2H(double p,double t);
    static double PT2S(double p,double t);
    static double PT2V(double p,double t);
    static double PT2U(double p,double t);
    static double PT2Cp(double p,double t);
    static double PT2Cv(double p,double t);
    static double PT2W(double p,double t);
public: //static function: P,H to 
    static double PH2T(double p,double h,int& itera);
    static double PH2T(double p,double h);
    static double PH2S(double p,double h);
    static double PH2V(double p,double h);
    static double PH2U(double p,double h);
    static double PH2Cp(double p,double h);
    static double PH2Cv(double p,double h);
    static double PH2W(double p,double h);
public: //static function: P,S to 
    static double PS2T(double p,double s,int& itera);
    static double PS2T(double p,double s);
    static double PS2H(double p,double s);
    static double PS2V(double p,double s);
    static double PS2U(double p,double s);
    static double PS2Cp(double p,double s);
    static double PS2Cv(double p,double s);
    static double PS2W(double p,double s);
public: //static function: P,V to 
    static double PV2T(double p,double v,int&);
    static double PV2T(double p,double v);
    static double PV2H(double p,double v);
    static double PV2S(double p,double v);
    static double PV2U(double p,double v);
    static double PV2Cp(double p,double v);
    static double PV2Cv(double p,double v);
    static double PV2W(double p,double v);
public: //static function: P,U to 
    static double PU2T(double p,double u,int&);
    static double PU2T(double p,double u);
    static double PU2H(double p,double u);
    static double PU2S(double p,double u);
    static double PU2V(double p,double u);
    static double PU2Cp(double p,double u);
    static double PU2Cv(double p,double u);
    static double PU2W(double p,double u);
public: //static function: P,Cp to 
        //由于过冷水相同压力下不同温度的Cp相差很小，PCp2T迭代出的t值很难还原回原来的t
        //不建议采用P，Cp to求其他参数
    static double PCp2T(double p,double cp,int&);
    static double PCp2T(double p,double cp);
    static double PCp2H(double p,double cp);
    static double PCp2S(double p,double cp);
    static double PCp2V(double p,double cp);
    static double PCp2U(double p,double cp);
    static double PCp2Cv(double p,double cp);
    static double PCp2W(double p,double cp);
public: //static function: T,H to 
    static double TH2P(double t,double h,double& p2, int& itera);
    static double TH2P(double t,double h,double& p2);
    static double TH2S(double t,double h);
    static double TH2V(double t,double h);
    static double TH2U(double t,double h);
    static double TH2Cp(double t,double h);
    static double TH2Cv(double t,double h);
    static double TH2W(double t,double h);
public: //static function: T,S to 
    static double TS2P(double t,double s,double& p2,int& itera);
    static double TS2P(double t,double s,double& p2);
    static double TS2H(double t,double s);
    static double TS2V(double t,double s);
    static double TS2U(double t,double s);
    static double TS2Cp(double t,double s);
    static double TS2Cv(double t,double s);
    static double TS2W(double t,double s);
public: //static function: T,V to 
    static double TV2P(double t,double v,int& itera);
    static double TV2P(double t,double v);
    static double TV2H(double t,double v);
    static double TV2S(double t,double v);
    static double TV2U(double t,double v);
    static double TV2Cp(double t,double v);
    static double TV2Cv(double t,double v);
    static double TV2W(double t,double v);
public: //static function: T,U to 
    static double TU2P(double t,double u,int& itera);
    static double TU2P(double t,double u);
    static double TU2H(double t,double u);
    static double TU2S(double t,double u);
    static double TU2V(double t,double u);
    static double TU2Cp(double t,double u);
    static double TU2Cv(double t,double u);
    static double TU2W(double t,double u);
public: //static function: T,Cp to 
    static double TCp2P(double t,double Cp,int& itera);
    static double TCp2P(double t,double Cp);
    static double TCp2H(double t,double Cp);
    static double TCp2S(double t,double Cp);
    static double TCp2V(double t,double Cp);
    static double TCp2U(double t,double Cp);
    static double TCp2Cv(double t,double Cp);
    static double TCp2W(double t,double Cp);
private:
    const static double ERR0;
    const static double ERR;
    const static double ERR2;
    const static double T0;
    const static double R;
    const static double Tc;
    const static double pc;
    const static double rhoc;
    const static double p_base;
    const static double T_base;
    const static double ni[34];
    const static double Ii[34];
    const static double Ji[34];
    const static double ni_tab6[20];
    const static double Ii_tab6[20];
    const static double Ji_tab6[20];
    const static double ni_tab8[20];
    const static double Ii_tab8[20];
    const static double Ji_tab8[20];
};
const double IF97Region1::ERR0=1E-5;
const double IF97Region1::ERR=1E-7;
const double IF97Region1::ERR2=1E-10;
const double IF97Region1::T0=273.15;
const double IF97Region1::R=0.461526;//page5 (1)
const double IF97Region1::Tc=647.096;//page5 (2)
const double IF97Region1::pc=22.064;//page5 (3)
const double IF97Region1::rhoc=322;//page5 (4)
const double IF97Region1::p_base=16.53;//page6
const double IF97Region1::T_base=1386;//page6
const double IF97Region1::ni[34]={ //page7 table2
             0.14632971213167E0,
            -0.84548187169114E0,
            -0.37563603672040E1,
             0.33855169168385E1,
            -0.95791963387872E0,
             0.15772038513228E0,
            -0.16616417199501E-1,
             0.81214629983568E-3,
             0.28319080123804E-3,
            -0.60706301565874E-3,
            -0.18990068218419E-1,
            -0.32529748770505E-1,
            -0.21841717175414E-1,
            -0.52838357969930E-4,
            -0.47184321073267E-3,
            -0.30001780793026E-3,
             0.47661393906987E-4,
            -0.44141845330846E-5,
            -0.72694996297594E-15,
            -0.31679644845054E-4,
            -0.28270797985312E-5,
            -0.85205128120103E-9,
            -0.22425281908000E-5,
            -0.65171222895601E-6,
            -0.14341729937924E-12,
            -0.40516996860117E-6,
            -0.12734301741641E-8,
            -0.17424871230634E-9,
            -0.68762131295531E-18,
             0.14478307828521E-19,
             0.26335781662795E-22,
            -0.11947622640071E-22,
             0.18228094581404E-23,
            -0.93537087292458E-25
 };
const double IF97Region1::Ii[34]={//page7 table2
            0,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,3,4,4,4,5,8,8,21,23,29,30,31,32
 };
const double IF97Region1::Ji[34]={//page7 table2
            -2,-1,0,1,2,3,4,5,-9,-7,-1,0,1,3,-3,0,1,3
            ,17,-4,0,6,-5,-2,10,-8,-11,-6,-29,-31,-38,-39,-40,-41
 };
const double IF97Region1::ni_tab6[20]={//page10 table6
			-0.23872489924521E3,
			 0.40421188637945E3,
			 0.11349746881718E3,
			-0.58457616048039E1,
			-0.15285482413140E-3,
			-0.10866707695377E-5,
			-0.13391744872602E2,
			 0.43211039183559E2,
			-0.54010067170506E2,
			 0.30535892203916E2,
			-0.65964749423638E1,
			 0.93965400878363E-2,
			 0.11573647505340E-6,
			-0.25858641282073E-4,
			-0.40644363084799E-8,
			 0.66456186191635E-7,
			 0.80670734103027E-10,
			-0.93477771213947E-12,
			 0.58265442020601E-14,
			-0.15020185953503E-16
 };
const double IF97Region1::Ii_tab6[20]={//page10 table6
			0,0,0,0,0,0,1,1,1,1,1,1,1,2,2,3,3,4,5,6
 };
const double IF97Region1::Ji_tab6[20]={//page10 table6
			0,1,2,6,22,32,0,1,2,3,4,10,32,10,32,10,32,32,32,32
 };
const double IF97Region1::ni_tab8[20]={//page11 table8
			 0.17478268058307E3,
			 0.34806930892873E2,
			 0.65292584978455E1,
			 0.33039981775489E0,
			-0.19281382923196E-6,
			-0.24909197244573E-22,
			-0.26107636489332E0,
			 0.22592965981586E0,
			-0.64256463395226E-1,
			 0.78876289270526E-2,
			 0.35672110607366E-9,
			 0.17332496994895E-23,
			 0.56608900654837E-3,
			-0.32635483139717E-3,
			 0.44778286690632E-4,
			-0.51322156908507E-9,
			-0.42522657042207E-25,
			 0.26400441360689E-12,
			 0.78124600459723E-28,
			-0.30732199903668E-30
 };
const double IF97Region1::Ii_tab8[20]={//page11 table8
			0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,4
 };
const double IF97Region1::Ji_tab8[20]={//page11 table8
			0,1,2,3,11,31,0,1,2,3,12,31,0,1,2,9,31,10,32,32
 };
