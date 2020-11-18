#include <stdio.h>
#include <math.h>

void main()
{
    double x[10];
    double y[10];
    double z[10];
    int i;
    for(i=0;i<10;i++)
    {
        x[i]=i/2.;
        y[i]=i*i;
    }
    for(i=0;i<10;i++)
    {
        if(i==0) z[i]=(y[i+1]-y[i])/(x[i+1]-x[i]);
        else z[i]=(y[i]-y[i-1])/(x[i]-x[i-1]);
    }
    for(i=0;i<10;i++)
    {
        printf("%f %f %f \n",x[i],y[i],z[i]);
    }
}