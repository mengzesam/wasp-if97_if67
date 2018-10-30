/*
    reference book:IAPWS R7-97(2012)
*/
#include <math.h>
class IF97Region4{
private://构造函数为private，禁止实例化IF97Region1
    IF97Region4(){};
public: //static function
    static double P2T(double ps);
    static double T2P(double ts);
   
private:
    const static double ni[10];
};
const double IF97Region4::ni[10]={ //page34 table34
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
double IF97Region4::P2T(double ps){
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
double IF97Region4::T2P(double ts){
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


