#include "AGA892DC.h"

void AGA892DC::sortInput(){
    if(status<1) return;
    ID[0]=ID_raw[0];
    Index2raw[0]=0;
    for(int i=1;i<NUM_x;i++){
        int key=ID_raw[i];
        int j=i-1;
        while(j>=0 && key<ID[j]){
            ID[j+1]=ID[j];
            Index2raw[j+1]=Index2raw[j];
            j--;
        }
        ID[j+1]=key;
        Index2raw[j+1]=i;
    }
    status=2;
}

double AGA892DC::calcB(){
    if(status<1) return -999999.0;
    if(status<2)
        sortInput();
    double T=t+273.15;
    double B=0.0;
    for(int n=0;n<18;n++){
        double sum_r=0.0;
        for(int i=0;i<NUM_x;i++){
            int ii=ID[i];
            int Fi_sqrt=0;
            int W_i=0;
            double Q_i=0.0;
            double S_i=0.0;
            double K_i=Ki[ii-1];
            if(ii==3)
                Q_i=Qi3;
            else if(ii==6){
                Q_i=Qi6;
                S_i=Si6;
                W_i=1;                
            }else if(ii==7){
                Q_i=Qi7;
                S_i=Si7;  
            }else if(ii==8)
                Fi_sqrt=1;
            double G_i=Gi[ii];
            double E_i=Ei[ii];
            double xi=xi_raw[Index2raw[i]];
            for(int j=0;i<NUM_x;j++){
                double Bnij;
                double Gij;
                double Eij;
                int jj=ID[j];
                int Fj_sqrt=0;
                int W_j=0;
                double Q_j=0.0;
                double S_j=0.0;
                double K_j=Ki[jj-1];
                if(jj==3)
                    Q_j=Qi3;
                else if(jj==6){
                    Q_j=Qi6;
                    S_j=Si6;
                    W_j=1;                
                }else if(jj==7){
                    Q_j=Qi7;
                    S_j=Si7;  
                }else if(jj==8)
                    Fj_sqrt=1;
                double xj=xi_raw[Index2raw[j]];
                double G_j=Gi[ii];
                double E_j=Ei[ii];
                double E_ij=1;     
                double G_ij=1;    
                int ki=ii;
                int kj=jj;
                if(ii>jj){
                    ki=jj;
                    kj=ii;
                }    
                if(ki==1 && kj>=2 && kj<=19){
                    E_ij=E1j[kj-2];
                    if(kj==3)
                        G_ij=G13;
                    else if(kj==8)
                        G_ij=G18;
                }else if(ki==2 && kj>=3 && kj<=14){
                    E_ij=E2j[kj-3];
                    if(kj==3)
                        G_ij=G23;
                }else if(ki==3 && kj>=4 && kj<=19){
                    E_ij=E3j[kj-4];
                    if(kj<=6)
                        G_ij=G3j[kj-4];
                }else if(ki==4 && kj>=5 && kj<=14){
                    E_ij=E4j[kj-5];
                }else if(ki==5){
                    if(kj==8)
                        E_ij=E58;
                    else if(kj==12)
                        E_ij=E512;
                }else if(ki==7 && kj>=15 && kj<=19){
                    E_ij=E7j[kj-15];
                }else if(ki==5){
                    if(kj==9)
                        E_ij=E89;
                    else if(kj==11)
                        E_ij=E811;
                    else if(kj==12)
                        E_ij=E812; 
                }
                Gij=G_ij*(G_i+G_j)/2.0;
                Eij=E_ij*sqrt(E_i*E_j);
                Bnij=pow(1.0+Gij-gn[n],gn[n])*pow(Q_i+Q_j+1.0-qn[n],qn[n])*pow(1.0+Fi_sqrt*Fj_sqrt-fn[n],fn[n])
                     *(S_i*S_j+1.0-sn[n],sn[n])*pow(1.0+W_i*W_j-wn[n],wn[n]);
                sum_r+=xi*xj*Bnij*pow(Eij,un[n])*pow(K_i*K_j,1.5);       
            }
        }
        B+=an[n]*pow(T,-un[n])*sum_r;
    }
    status=3;
    return B;
}

double AGA892DC::calcK(){
    if(status<1) return -999999.0;
    if(status<2)
        sortInput();
    double K=0.0;
    double left=0.0;
    for(int i=0;i<NUM_x;i++){
        int ii=ID[i];
        double xi=xi_raw[Index2raw[i]];
        double K_i=Ki[ii-1];
        left+=xi*pow(K_i,2.5);
    }
    double right=0.0;
    for(int i=0;i<NUM_x-1;i++){
        int ii=ID[i];
        if((ii>=1 && ii<=4) || ii==7){
            double xi=xi_raw[Index2raw[i]];
            double K_i=Ki[ii-1];
            for(int j=i+1;j<NUM_x;j++){
                int jj=ID[j];
                double K_ij=1.0;
                if(ii==1 && jj<=19){
                    K_ij=K1j[jj-2];
                }else if(ii==2 && jj<=8){
                    K_ij=K2j[jj-3];
                }else if(ii==3 && jj<=19){
                    K_ij=K3j[jj-4];
                }else if(ii==4 && jj<=8){
                    K_ij=K3j[jj-5];
                }else if(ii==7 && jj>=15 && jj<=19){
                    K_ij=K7j[jj-15];
                }
                if(abs(K_ij-1.0)<1E-7){//K_ij!=1.0
                    double xj=xi_raw[Index2raw[j]];
                    double K_j=Ki[jj-1];
                    right+=xi*xj*(K_ij*K_ij*K_ij*K_ij*K_ij-1.0)*pow(K_i*K_j,2.5);
                }//else right+=0; 
            }
        }//else right+=0;
    }
    K=pow(left*left+2*right,0.2);
    status=4;
    return K;
}

double AGA892DC::calcF(){
    if(status<1) return -999999.0;
    if(status<2)
        sortInput();
    double F=0.0;
    status=5;
    for(int i=0;i<NUM_x;i++){
        int ii=ID[i];
        if(ii==8){
            F=xi_raw[Index2raw[i]];
            return F*F;
        }
    }
    return F;
}

double AGA892DC::calcQ(){
    if(status<1) return -999999.0;
    if(status<2)
        sortInput();
    double Q=0.0;
    status=6;
    for(int i=0;i<NUM_x;i++){
        int ii=ID[i];
        if(ii==3){
            Q+=xi_raw[Index2raw[i]]*Qi3;
        }else if(ii==6){
            Q+=xi_raw[Index2raw[i]]*Qi6;
        }else if(ii==7){
            Q+=xi_raw[Index2raw[i]]*Qi7;
        }
    }
    return Q;
}

double AGA892DC::calcG(){
    if(status<1) return -999999.0;
    if(status<2)
        sortInput();
    double G=0.0;    
    for(int i=0;i<NUM_x;i++){
        int ii=ID[i];
        double xi=xi_raw[Index2raw[i]];
        double G_i=Gi[ii-1];
        G+=xi*G_i;
    }
    for(int i=0;i<NUM_x-1;i++){
        int ii=ID[i];
        if(ii>=1 && ii<=3){
            double xi=xi_raw[Index2raw[i]];
            double G_i=Gi[ii-1];            
            for(int j=i+1;j<NUM_x;j++){
                int jj=ID[j];
                double G_ij=1.0;
                if(ii==1){
                    if(jj==3)
                        G_ij=G13;
                    else if(jj==8){
                        G_ij=G18;
                    }
                }else if(ii==2 && jj==3){
                    G_ij=G23;
                }else if(ii==3){
                    if(jj==4){
                        G_ij=G3j[0];
                    }else if(jj==6){
                        G_ij=G3j[2];
                    }
                }            
                if(abs(G_ij-1.0)<1E-7){//G_ij!=1.0
                    double xj=xi_raw[Index2raw[j]];
                    double G_j=Gi[jj-1];
                    G+=xi*xj*(G_ij-1.0)*(G_i+G_j);
                }//else G+=0;    
            }
        }//else G+=0;
    }
    status=7;
    return G;
}    

double AGA892DC::calcU(){
    if(status<1) return -999999.0;
    if(status<2)
        sortInput();
    double U=0.0;
    double left=0.0;        
    for(int i=0;i<NUM_x;i++){
        int ii=ID[i];
        double xi=xi_raw[Index2raw[i]];
        double E_i=Ei[ii-1];
        left+=xi*E_i;
    }    
    double right=0.0;
    for(int i=0;i<NUM_x-1;i++){
        int ii=ID[i];
        if((ii>=1 && ii<=4) || ii==7){
            double xi=xi_raw[Index2raw[i]];
            double E_i=Ki[ii-1];
            for(int j=i+1;j<NUM_x;j++){
                int jj=ID[j];
                double U_ij=1.0;
                if(ii==1 && jj<=19){
                    U_ij=U1j[jj-2];
                }else if(ii==2 && jj<=12){
                    U_ij=U2j[jj-3];
                }else if(ii==3 && jj<=19){
                    U_ij=U3j[jj-4];
                }else if(ii==4 && jj<=14){
                    U_ij=U3j[jj-5];
                }else if(ii==7 && jj>=15 && jj<=19){
                    U_ij=U7j[jj-15];
                }
                if(abs(U_ij-1.0)<1E-7){//U_ij!=1.0
                    double xj=xi_raw[Index2raw[j]];
                    double E_j=Ei[jj-1];
                    right+=xi*xj*(U_ij*U_ij*U_ij*U_ij*U_ij-1.0)*pow(E_i*E_j,2.5);
                }//else right+=0; 
            }
        }//else right+=0;
    }
    U=pow(left*left+right,0.2);
    status=8;
    return U;
}

double AGA892DC::rhom2P(double rhom,double tt,double B,double K,double F,double Q,double G,double U){    
    double rhor=K*K*K*rhom;
    double T=t+273.15;
    double R=8.314510;
    double sum1=0.0;
    double sum2=0.0;
    for(int n=12;n<18;n++){
        double Cn;
        Cn=an[n]*pow(G+1.0-gn[n],gn[n])*pow(Q*Q+1.0-qn[n],qn[n])*pow(F+1.0-fn[n],fn[n])
            *pow(U,un[n])*pow(T,-un[n]);
        sum1+=Cn;
        sum2+=Cn*(bn[n]-cn[n]*kn[n]*pow(rhor,kn[n]))*pow(rhor,bn[n])*exp(-cn[n]*pow(rhor,kn[n]));
    }
    for(int n=18;n<58;n++){
        double Cn;
        Cn=an[n]*pow(G+1.0-gn[n],gn[n])*pow(Q*Q+1.0-qn[n],qn[n])*pow(F+1.0-fn[n],fn[n])
            *pow(U,un[n])*pow(T,-un[n]);
        sum2+=Cn*(bn[n]-cn[n]*kn[n]*pow(rhor,kn[n]))*pow(rhor,bn[n])*exp(-cn[n]*pow(rhor,kn[n]));
    }
    p=rhom*R*T*(1.0+B*rhom-rhor*sum1+sum2);
    return p;
}

double AGA892DC::calcRhom(){
    if(status<1) return -999999.0;
    if(status<2)
        sortInput();
    double err=1E-6;
    double B=calcB();
    double K=calcK();
    double F=calcF();
    double Q=calcQ(); 
    double G=calcG();  
    double U=calcU();
    double rhom=0.6;
    double pp=rhom2P(rhom,t,B,K,F,Q,G,U);
    if(abs(pp-p)>err){
        double rhom0=rhom;
        double p0=pp;
        double rhom1=1.0;
        double p1=rhom2P(rhom1,t,B,K,F,Q,G,U);
        rhom=rhom1+(p-p1)/(p0-p1)*(rhom0-rhom1);
        pp=rhom2P(rhom,t,B,K,F,Q,G,U);
        while(abs(pp-p)>err){
            rhom0=rhom1;
            p0=p1;
            rhom1=rhom;
            p1=pp;
            rhom=rhom1+(p-p1)/(p0-p1)*(rhom0-rhom1);
            pp=rhom2P(rhom,t,B,K,F,Q,G,U);
        }
    }
    return rhom;
}



int main(int argc, char const *argv[])
{
    //AGA892DC::sortInput();
    return 0;
}
