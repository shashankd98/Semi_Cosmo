#include <stdio.h>
#include <math.h>
#define pi 3.14159

double Om=0.3089;
double H0=70.e5/3.086e24;
double G=1.0;
// double dcr0 = 3.0*H0*H0/(8.*pi*G);
double Msun=2.e33; 
double M200=1;
double R200m=1.48;
double R200c=1; 
double c=4.0;
double nu=4.0;
double be = 1.5; 
double se=1.5; 
// double rt=(1.9-0.18*nu)*R200m; 


double phi_nfw(double r)
{
    double x=r/R200c; double rs=R200c/c;
    return -G*M200*log(1.+c*x)/(r*(log(1.+c)-c/(1.+c)));
}

double phi_outer(double r)
{
    double const_outer;
    const_outer=-0.1*phi_nfw(R200c)/((be*pow((5.*R200m),(se))*(pow(R200c,(2.-se))/((3.-se)*(2.-se)))+pow(R200c,2.)/6.));
    return const_outer*((be*pow((5.*R200m),(se))*(pow(r,(2.-se))/((3.-se)*(2.-se)))+pow(r,2.)/6.));
}

double sigmoid(double x,double ce)
{
    return (1.+erf(x))/ce;
}

double phi_full(double r)
{
    double rt=(1.9-0.18*nu)*R200m;
    double w1 =1.-sigmoid(log(r/rt),7.); 
    double w2 =pow((1.-pow(w1,4.)),(1./4.));
    double phi_fit = phi_nfw(r)*w1+phi_outer(r)*w2;
    return phi_fit;
}

void main()
{
    // double a=1.12344e-15;
    // double b=1e-15;
    // double c=1e20; 
    // printf("%g \n",a*b*c);
    int len=1000;
    double r[len];
    double phi[len];
    double g[len];
    double rho[len];
    double slope[len];
    int i;
    for(i=0;i<len;i++)
    {
        double ind=i/(double)len*5.-2.;
        r[i]=R200m*pow(10.,ind);
        phi[i]=phi_full(r[i]);
        // printf("%f \t %f \t %f \t %f \n",r[i],-phi_nfw(r[i]),phi_outer(r[i]),phi_full(r[i]));
    }
    for(i=0;i<len;i++)
    {
        if(i==0) g[i]=(phi[i+1]-phi[i])/(r[i+1]-r[i]);
        else if(i==len-1) g[i]=(phi[i]-phi[i-1])/(r[i]-r[i-1]);
        else g[i]=((phi[i+1]-phi[i])/(r[i+1]-r[i])+(phi[i]-phi[i-1])/(r[i]-r[i-1]))/2.;
    }
    for(i=0;i<len;i++)
    {
        if(i==0) rho[i]=(g[i+1]*pow(r[i+1],2)-g[i]*pow(r[i],2))/(pow(r[i+1],3)-pow(r[i],3));
        else if(i==len-1) rho[i]=(g[i]*pow(r[i],2)-g[i-1]*pow(r[i-1],2))/(pow(r[i],3)-pow(r[i-1],3));
        else rho[i]=((g[i+1]*pow(r[i+1],2)-g[i]*pow(r[i],2))/(pow(r[i+1],3)-pow(r[i],3))+(g[i]*pow(r[i],2)-g[i-1]*pow(r[i-1],2))/(pow(r[i],3)-pow(r[i-1],3)))/2.;
    }
    for(i=0;i<len;i++)
    {
        if(i==0) slope[i]=(log(rho[i+1])-log(rho[i]))/(log(r[i+1])-log(r[i]));
        else if(i==len-1) slope[i]=(log(rho[i])-log(rho[i-1]))/(log(r[i])-log(r[i-1]));
        else slope[i]=((log(rho[i+1])-log(rho[i]))/(log(r[i+1])-log(r[i]))+(log(rho[i])-log(rho[i-1]))/(log(r[i])-log(r[i-1])))/2.;
    }
    for(i=0;i<len;i++)
    {
        // printf("%f \t %f \t %f \t %f \t %f \n",r[i],phi[i],g[i],rho[i],slope[i]);
    }
}