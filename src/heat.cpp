/* 
 * File:   heat.cpp
 * Author: Nikhil
 * 
 * Created on May 18, 2015, 1:16 PM
 */

#include "../include/heat.hpp"
#include <cmath>

heat::heat(int nodes, double dX):ODESolver(nodes)
{
    dx = dX;
}


void heat::fcn(double *fvalue, const double *y, const double time){
    /*fvalue[0] = -1000*y[0]+3000-2000*exp(-time);
    for(int i=1; i<n-1; i++){
        fvalue[i] = -1000*y[i]+3000-2000*exp(-time);
    }
    fvalue[n-1] = -1000*y[n-1]+3000-2000*exp(-time);*/
    
    for(int i=0; i<n; i++){
        fvalue[i] = -50*(y[i] - cos(time));
    }
    
    /*fvalue[0] = (y[1] - 2.0*y[0] + 500)/(dx*dx);
    for(int i=1; i<n-1; i++){
        fvalue[i] = (y[i+1] - 2.0*y[i] + y[i-1])/(dx*dx);;
    }
    fvalue[n-1] = (200 - 2.0*y[n-1] + y[n-2])/(dx*dx);*/
    
}

