#include <stdio.h>
#include <math.h>
#define pi 3.14159
#define pc 3.086e18
#define Msun 2.e33

double Om=0.3089;
double H0=70.e5/3.086e24;
double G=6.67e-8;
double M200=1e14*Msun;
double c=4.0;
double nu=4.0;
double be = 1.5; 
double se=1.5; 
double redshift=0.0;
double gam=4.0;
double bet=6.0;

double Ez(double z)
{
    return sqrt( Om*pow((1.+z),3.)+(1.-Om) );
}

double Ezsqr(double z)
{
    return Om*pow((1.+z),3.)+(1.-Om);
}

double delc(double x) 
{
    return 18.*pow(pi,2) + 82.*x - 39.*x*x;
}

double d_NFW(double r)
{
    double dcr0=3.*H0*H0/(8.*pi*G);
    double dcr=dcr0*Ezsqr(redshift);
    double dm=dcr0*Om*pow((1.+redshift),3.);
    double R200c = pow((3.*M200/(800.*pi*dcr)),(1./3.));
    double R200m = pow((3.*M200/(800.*pi*dm)),(1./3.));
    double rs=R200c/c;
    return 200.*c*c*c*dcr/(3.*(log(1.+c)-c/(1.+c)))/(r/rs*pow((1+r/rs),2.));
}

double f_tr(double r)
{
    double dcr0=3.*H0*H0/(8.*pi*G);
    double dcr=dcr0*Ezsqr(redshift);
    double dm=dcr0*Om*pow((1.+redshift),3.);
    double R200c = pow((3.*M200/(800.*pi*dcr)),(1./3.));
    double R200m = pow((3.*M200/(800.*pi*dm)),(1./3.));
    double rt=(1.9-0.18*nu)*R200m; 
    return pow((1. + pow((r/rt),bet)),(-gam/bet));
}

double d_out(double r)
{
    double dcr0=3.*H0*H0/(8.*pi*G);
    double dcr=dcr0*Ezsqr(redshift);
    double dm=dcr0*Om*pow((1.+redshift),3.);
    double R200c = pow((3.*M200/(800.*pi*dcr)),(1./3.));
    double R200m = pow((3.*M200/(800.*pi*dm)),(1./3.));
    return dm*(be*pow((r/(5.*R200m)),(-se)) + 1.);
}

double d_rho(double r)
{
    return d_NFW(r)*f_tr(r)+d_out(r);
}

double phi_nfw(double r)
{
    double dcr0=3.*H0*H0/(8.*pi*G);
    double dm=dcr0*Om*pow((1.+redshift),3.);
    double dcr=dcr0*Ezsqr(redshift);
    double R200c = pow((3.*M200/(800.*pi*dcr)),(1./3.));
    double x=r/R200c;
    return -G*M200*log(1.+c*x)/(r*(log(1.+c)-c/(1.+c)));
}

double phi_outer(double r)
{
    double dcr0=3.*H0*H0/(8.*pi*G);
    double dcr=dcr0*Ezsqr(redshift);
    double dm=dcr0*Om*pow((1.+redshift),3.);
    double R200c = pow((3.*M200/(800.*pi*dcr)),(1./3.));
    double R200m = pow((3.*M200/(800.*pi*dm)),(1./3.));
    return (4.*pi*G)*dm*((be*pow((5.*R200m),(se))*(pow(r,(2.-se))/((3.-se)*(2.-se)))+pow(r,2.)/6.));
}

double sigmoid(double x,double ce)
{
    return (1.+erf(x))/ce;
}

double phi_full(double r)
{
    double dcr0=3.*H0*H0/(8.*pi*G);
    double dcr=dcr0*Ezsqr(redshift);
    double dm=dcr0*Om*pow((1.+redshift),3.);
    double R200c = pow((3.*M200/(800.*pi*dcr)),(1./3.));
    double R200m = pow((3.*M200/(800.*pi*dm)),(1./3.));

    double rt=(1.9-0.18*nu)*R200m;
    double w1 =1.-sigmoid(log(r/rt),7.); 
    double w2 =pow((1.-pow(w1,4.)),(1./8.));
    double phi_fit = phi_nfw(r)*w1+(phi_outer(r)-2./3.*pi*G*dm*pow(r,2))*w2;
    return phi_fit;
}

void main()
{
    double dcr0=3.*H0*H0/(8.*pi*G);
    double dm=dcr0*Om*pow((1.+redshift),3.);
    double dcr=dcr0*Ezsqr(redshift);
    double R200c = pow((3.*M200/(800.*pi*dcr)),(1./3.));
    double R200m = pow((3.*M200/(800.*pi*dm)),(1./3.));
    // printf("%g \n ",R200m/1.e6/pc);

    int len=1000;
    double r[len];
    double phi[len];
    double rho_true[len];
    double g[len];
    double rho[len];
    double slope[len];
    int i;
    for(i=0;i<len;i++)
    {
        double ind=i/(double)len*5.-2.;
        r[i]=R200m*pow(10.,ind);
        phi[i]=phi_full(r[i]);
        rho_true[i]=d_rho(r[i]);
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
        if(i==0) rho[i]=(g[i+1]*pow(r[i+1],2)-g[i]*pow(r[i],2))/(pow(r[i+1],3)/3.-pow(r[i],3)/3.)/(4.*pi*G);
        else if(i==len-1) rho[i]=(g[i]*pow(r[i],2)-g[i-1]*pow(r[i-1],2))/(pow(r[i],3)/3.-pow(r[i-1],3)/3.)/(4.*pi*G);
        else rho[i]=((g[i+1]*pow(r[i+1],2)-g[i]*pow(r[i],2))/(pow(r[i+1],3)/3.-pow(r[i],3)/3.)+(g[i]*pow(r[i],2)-g[i-1]*pow(r[i-1],2))/(pow(r[i],3)/3.-pow(r[i-1],3)/3.))/2./(4.*pi*G);
    }
    for(i=0;i<len;i++)
    {
        if(i==0) slope[i]=(log(rho[i+1])-log(rho[i]))/(log(r[i+1])-log(r[i]));
        else if(i==len-1) slope[i]=(log(rho[i])-log(rho[i-1]))/(log(r[i])-log(r[i-1]));
        else slope[i]=((log(rho[i+1])-log(rho[i]))/(log(r[i+1])-log(r[i]))+(log(rho[i])-log(rho[i-1]))/(log(r[i])-log(r[i-1])))/2.;
    }
    for(i=0;i<len;i++)
    {
        printf("%g \t %g \t %g \t %g \t %g \t %g \n",r[i]/R200m,phi[i]-phi[0],g[i],rho[i]+dm,slope[i],rho_true[i]);
    }
}