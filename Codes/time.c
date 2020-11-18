#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define sec_year 3.154e7
#define pc 3.086e18
#define Msun 2.e33

double h=0.7;
double omega_m=0.25;
double omega_l=0.75;
double M0=1.e11*Msun;
double zf=1.82228585400769;
// double zf=2.129832443766769;
double z_start=6.0;

double t(double z)
{
    double H0=100.*1.e5/3.086e24*h;
    double omega_l=1.-omega_m;
    // double time=1./H0*2./3./sqrt(omega_l)*log((sqrt(omega_l/pow((1.+z),3.))+sqrt(omega_l/pow(1.+z,3.)+omega_m))/sqrt(omega_m));
    double time=2./3./H0/sqrt(1.-omega_m)*asinh(sqrt(1.-omega_m)/(sqrt(omega_m*pow((1.+z),3.))));
    return time;
}

double t_inv(double z,double t0)
{
    double t_start=t(z_start);
    return t(z)-t(z_start)-t0;
}

double M(double z)
{
    double nu=1.211+1.858*log10(1.+zf)+0.308*pow(omega_l,2.)-0.032*log10(M0/(1.e11/h*Msun));
    double rhs=-0.301*pow((log10(1.+z)/log10(1.+zf)),nu);
    // printf("%g %g\n",z,rhs);
    return M0*pow(10.,(rhs));
}

double solver(double t0)
{
    // double t_start=t(z_start);
    // t0=t_start-t0;
    double left,right,mid,eps=0.0001;
    double midval,rightval,root;
    left=0.;
    right=10.;
    do
    {
        mid=(left+right)/2.;
        rightval=t_inv(right,t0);
        midval=t_inv(mid,t0);
        // if(abs(midval)<eps) break;
        if(rightval*midval>0) 
        {right=mid;}
        else {left=mid;}

    } while (right-left>eps);
    root=(right+left)/2.;
    return root;
}

double ct(double time)
{
    double nu=1.211+1.858*log10(1.+zf)+0.308*pow(omega_l,2.)-0.032*log10(M0/(1.e11/h*Msun));
    double z4=pow(10.,(pow((log10(0.04)/-0.301),(1./nu))*log10(1.+zf)))-1.;
    double t4=t(z4);
    // printf("%g \t %g \n",time/1e9/sec_year,t4/1e9/sec_year);
    double c=4.*pow((1.+pow((time/3.75/t4),(8.4))),(1./8.));
    return c;
}

double z4()
{
    double nu=1.211+1.858*log10(1.+zf)+0.308*pow(omega_l,2.)-0.032*log10(M0/(1.e11/h*Msun));
    double z04=pow(10.,pow((log10(0.04)/(-0.301)),(1./nu))*log10(1.+zf))-1.;
    return z04;
}

void main()
{
    double t_start=t(0)-t(z_start);
    // printf("%g \n",t_start/sec_year/1e9);
    // printf("%g \n",t(0)/sec_year/1e9);
    double z[100];
    int i;
    for(i=0;i<100;i++)
    {
        z[i]=(double)i/10.;
        double time=t(z[i]);
        // printf("%g \t %g \n ",z[i],(t(z[i]))/1e9/sec_year);
        printf(" %g \t %g \t %g \n ",z[i],time/1e9/sec_year,ct(time));
        // printf("%g \t %g \n",time[i]/sec_year/1e9,solver(time[i]));
    }
    // printf("%g \n",solver(t(z_start)-t(z_start)));
}