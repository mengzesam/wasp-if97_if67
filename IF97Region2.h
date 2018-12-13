/*
    reference book:IAPWS R7-97(2012)
*/
#include "IF97Base.h"
class IF97Region2:public IF97Base{
private://构造函数为private,禁止实例化IF97Region3
    IF97Region2(){};
public: //static function: P,T to 
    static double PT2H(double p,double t);
    static double PT2S(double p,double t);
    static double PT2V(double p,double t);
    static double PT2U(double p,double t);
    static double PT2Cp(double p,double t);
    static double PT2Cv(double p,double t);
    static double PT2W(double p,double t);
public: //TODO:static function: P,T to  for metastable-vapor region
public://static function:P,H to
	static double PH2T(double p,double h,int& itera);
	static double PH2T(double p,double h);
public://static function:P,S to
	static double PS2T(double p,double s,int& itera);
	static double PS2T(double p,double s);
public://static function:P,V to
	static double PV2T(double p,double v,int & itera);
public://static function:P,U to
	static double PU2T(double p,double u,int & itera);
public://static function:P,Cp to
	static double PCp2T(double p,double cp,int & itera);
public://static function:P,Cv to
	static double PCv2T(double p,double cv,int & itera);
public://static function:T,H to
	static double TH2P(double t,double h,int & itera);
	static double TH2Pbeta(double t,double h,int & itera);
public://static function:T,S to
	static double TS2P(double t,double s,int & itera);
public://static function:T,V to
	static double TV2P(double t,double v,int & itera);
public://static function:T,Cp to
	static double TCp2P(double t,double ch,int & itera);
public://static function:T,Cv to
	static double TCv2P(double t,double cv,int & itera);
public://static function:H,S to
	static void HS2PT(double h,double s,double& p,double& t,int& itera);
private://static function:P,H or P,S to 辅助函数
	const static double S2bc;
	static double P2bc(double h);
	static double H2bc(double p);
	static double PH2T2a(double p,double h);
	static double PH2T2b(double p,double h);
	static double PH2T2c(double p,double h);
	static double PS2T2a(double p,double s);
	static double PS2T2b(double p,double s);
	static double PS2T2c(double p,double s);
private://static function:H,S to 辅助函数
	static double H2ab(double s);
	static double HS2P2a(double h,double s);
	static double HS2P2b(double h,double s);
	static double HS2P2c(double h,double s);
private:
    const static double ni_tab10[9];
    const static int Ji_tab10[9];
    const static double ni_tab11[43];
    const static int Ii_tab11[43];
    const static int Ji_tab11[43];
};
const double IF97Region2::ni_tab10[9]={//if97 table10 page13
			-0.96927686500217E1,
			 0.10086655968018E2,
			-0.56087911283020E-2,
			 0.71452738081455E-1,
			-0.40710498223928E0,
			 0.14240819171444E1,
			-0.43839511319450E1,
			-0.28408632460772E0,
			 0.21268463753307E-1
};
const int IF97Region2::Ji_tab10[9]={
			0,1,-5,-4,-3,-2,-1,2,3
};
const double IF97Region2::ni_tab11[43]={//if97 table11 page14
			-0.17731742473213E-2,
			-0.17834862292358E-1,
			-0.45996013696365E-1,
			-0.57581259083432E-1,
			-0.50325278727930E-1,
			-0.33032641670203E-4,
			-0.18948987516315E-3,
			-0.39392777243355E-2,
			-0.43797295650573E-1,
			-0.26674547914087E-4,
			 0.20481737692309E-7,
			 0.43870667284435E-6,
			-0.32277677238570E-4,
			-0.15033924542148E-2,
			-0.40668253562649E-1,
			-0.78847309559367E-9,
			 0.12790717852285E-7,
			 0.48225372718507E-6,
			 0.22922076337661E-5,
			-0.16714766451061E-10,
			-0.21171472321355E-2,
			-0.23895741934104E2,
			-0.59059564324270E-17,
			-0.12621808899101E-5,
			-0.38946842435739E-1,
			 0.11256211360459E-10,
			-0.82311340897998E1,
			 0.19809712802088E-7,
			 0.10406965210174E-18,
			-0.10234747095929E-12,
			-0.10018179379511E-8,
			-0.80882908646985E-10,
			 0.10693031879409E0,
			-0.33662250574171E0,
			 0.89185845355421E-24,
			 0.30629316876232E-12,
			-0.42002467698208E-5,
			-0.59056029685639E-25,
			 0.37826947613457E-5,
			-0.12768608934681E-14,
			 0.73087610595061E-28,
			 0.55414715350778E-16,
			-0.94369707241210E-6
};
const int IF97Region2::Ii_tab11[43]={
			1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,5,6,6,6,7,7,7,
            8,8,9,10,10,10,16,16,18,20,20,20,21,22,23,24,24,24
};
const int IF97Region2::Ji_tab11[43]={
			0,1,2,3,6,1,2,4,7,36,0,1,3,6,35,1,2,3,7,3,16,35,0,11,25,
            8,36,13,4,10,14,29,50,57,20,35,48,21,53,39,26,40,58
};
const double IF97Region2::S2bc=5.85;//if97 page21 fig2