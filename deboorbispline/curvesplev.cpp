#include <math.h>

void sortX(const double x[],double xs[],int index[],const int m){
    xs[0]=x[0];
    index[0]=0;
    for(int i=1;i<m;i++){
        double key=x[i];
        int j=i-1;
        while(j>=0 && key<xs[j]){
            xs[j+1]=xs[j];
            index[j+1]=index[j];
            j--;
        }
        xs[j+1]=key;
        index[j+1]=i;
    }
}
int findXPos_(const double x,const double knots[],int left,int right){    
    double err=1E-8;
    while(left<=right){
        int mid=left+(right-left)/2;
        if(abs(knots[mid]-x)<err)//==
            return mid;
        else if(x<knots[mid])
            return findXPos_(x,knots,left,mid-1);
        else
            return findXPos_(x,knots,mid+1,right);
    }
    if(left==0) return left;
    return left-1;
}
int findXPos(const double x,const double knots[],const int n,int start=0){
    //find the first position if knots[i] knots[i+1]... == x
    double err=1E-8;
    if(start>=n) start=0;
    int pos=findXPos_(x,knots,start,n-1);
    while(pos>start && abs(knots[pos]-knots[pos-1])<err){
        pos--;
    }
    return pos;
}
void bsplvb(const double knots[],const int n,const int k,const double x,const int j,double bj[]){
    //output:bj[] length k+1
    //pls refer to De Boor:"A Practical Guide to Spline"  Chapter X
    double err=1E-8;
    bj[0]=1.0;
    for(int i=1;i<=k;i++){
        double saved=0.0;
        for(int r=1;r<=i;r++){  
            double deltaR=knots[j+r]-x;
            double deltaL=x-knots[j+r-i];//deltaL(i+1-r): deltaL(s)=x-knots[j+1-s]         
            double term=0.0;
            if(abs(deltaL+deltaR)>err)//non-zero
                term=bj[r-1]/(deltaL+deltaR);
            bj[r-1]=saved+deltaR*term;
            saved=deltaL*term;
        }
        bj[i]=saved;
    }
}
int curvesplev(const double knots[],const double BCoeff[],const int n,const int k,
               const double x[],double y[],const int m){
/*
    evaluates in a number of points x(i),i=1,2,...,m
    a spline s(x) of degree k(order k+1), given in its b-spline representation
    return: 0 success
            -1 NOT 2*k>n
            -2 NOT knots[i]<=knots[i+1]
            less than -10: x[i]<tb or x[i]>te (i=-ret-11)
*/
    int ret=0;
    if(2*k>n) return -1;
    for(int i=1;i<n;i++)
        if(knots[i-1]>knots[i])
            return -2;
    double err=1E-8;
    double* x_sort=new double[m];
    int* index=new int[m];
    double* bj=new double[k+1];
    sortX(x,x_sort,index,m);
    double tb=knots[k]-err;
    double te=knots[n-k-1]+err;
    int j=0;
    for(int i=0;i<m;i++){
        double xx=x_sort[i];
        if(xx<tb){
            y[index[i]]=-999999.0;
            ret=-(10+index[i]+1);
            continue;
        }else if(xx>te){
            y[index[i]]=-999999.0;
            ret=-(10+index[i]+1);
            continue;
        }
        j=findXPos(xx,knots,n,j);
        if(j<k) 
            j=k;
        else if(j>=n-k-1)
            j=n-k-2;
        bsplvb(knots,n,k,xx,j,bj);//before calling, be sure that bj legth must be:k+1       
        y[index[i]]=0.0;
        for(int r=0;r<k+1;r++){
            y[index[i]]+=BCoeff[j-k+r]*bj[r];
        }
    }
    delete bj;
    delete index;
    delete x_sort;
    return ret;
}

int main(){
/*     double x[]={1,1,1,3,3,3,3,5,5,5,7,7,7,7,7,7,7,9,9,9};
    int n=20;
    double xx=11;
    int i=findXPos(xx,x,n,10); */
    double BCoeff[]={
        4.71828182846,5.16616728369,6.15131614325,7.92149863373,10.01334554157,12.46264981465,15.31312905827,18.61818007610,22.44302187144,26.86731311477,31.98834912596,37.92496667703,44.82231332913,52.85767271375,62.24757954739,73.25650993001,86.20749569959,101.49508883399,119.60119620720,141.11442020432,166.75368140194,197.39707137523,234.11709359598,278.22370676283,331.31689804331,395.35089617692,472.71260153612,566.31738081930,679.72607095294,817.28788798036,984.31497637275,1187.29560403894,1434.15455929846,1734.57120046262,2100.36792247009,2545.98463110700,3089.05826712812,3751.13063861680,4558.51296938186,5543.34186075389,6744.86904623604,8211.03670139155,10000.40153168560,12184.48485868220,14850.64302180680,18105.57329497630,22079.59602284130,26931.88483380220,32856.85483661370,40091.96518133660,48927.24912880590,56121.59893418460,60116.14171519780,60116.14171519780,60116.14171519780,60116.14171519780,60116.14171519780
    };
    double knots[]={
        1,1,1,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.2,4.4,4.6,4.8,5,5.2,5.4,5.6,5.8,6,6.2,6.4,6.6,6.8,7,7.2,7.4,7.6,7.8,8,8.2,8.4,8.6,8.8,9,9.2,9.4,9.6,9.8,10,10.2,10.4,10.6,10.8,11,11,11,11
    };
    int k=3;
    int n=57;
    int m=5;
    double x[5]={
        10.5,11,5.5,3.1,1.5
    };
    double y[10]={};
    int ret=curvesplev(knots,BCoeff,n,k,x,y,m);

    return 0;
}