#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double func(double x)
{
    return cos(x)-x;
}

void main()
{
    double left,right,mid,eps=0.0001;
    double midval,rightval,root;
    left=0.;
    right=3.;
    do
    {
        mid=(left+right)/2.;
        rightval=func(right);
        midval=func(mid);
        if(abs(midval)<eps) break;
        if(rightval*midval>0) right=mid;
        else left=mid;
        
        printf("%g %g \n",left,right);
    } while (abs(right-left)>eps);
    root=(right+left)/2.;
    printf("%g %g \n",root,func(root));
}