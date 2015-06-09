/* 
 * File:   ODESolver.cpp
 * Author: Nikhil
 * 
 * Created on May 16, 2015, 3:58 PM
 */

#include "../include/ODESolver.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "../include/lapackpp/lapackpp.h"


ODESolver::ODESolver(int size)
{
	n = size;
	jacob = new double [n*n];
	yp = new double [n];
	ym = new double [n];
	a = new double [n*n];
	b = new double [n];
        x = new double [n];
        dummyp = new double [n];
}

void ODESolver::initialize(double *y_var, double *y_dot, double *resifun)
{
	y = y_var; yd = y_dot, resifcn = resifun;
	for(int k=0; k<n; k++)
	{
		x[k] = y[k];
	}
}

void ODESolver::integrate(const int totalSteps, const double dt)
{
		tStep = 0;
                std::cout<< "value of n = "<<n<< "\ttStep :\t" << tStep << std::endl;
                LaGenMatDouble A(n,n); LaVectorDouble B(n), X(n);
	do
	{
                
		time = dt*tStep;
		jacobian(jacob, y, time);
		for(int i=0; i<n; i++)
		{
			for(int j =0; j<n; j++)
			{
				if(i==j){dummy = 1;}else{dummy = 0;}
				
				A(i,j) = dummy - dt*(*(jacob + i*n + j));
			}
		}
		for(int k=0; k<n; k++){resifcn[k] = fcn(y, tStep*dt, k);}
		for(int i=0; i<n; i++)
		{
			B(i) = dt*(*(resifcn + i)); 
		}
		
		LaLinearSolve(A,X,B);
		
		for(int k=0; k<n; k++){	y[k] = X(k)+y[k];}
		
		tStep++;
                console();
                write();
	}while(tStep < totalSteps);
}

void ODESolver::jacobian(double* jacob, double* yVar, const double time)
{
    double fluxPlus, fluxMinus;
	
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<n; j++)
		{
			yVar[j] = yVar[j] + 0.0001;
                        fluxPlus = fcn(y, time, i);
                        yVar[j] = yVar[j] - 0.0002;
                        fluxMinus = fcn(y, time, i);
                        yVar[j] = yVar[j] + 0.0001;
                        jacob[i*n + j] = (fluxPlus - fluxMinus)/0.0002; 
		}
	}
	
}


void ODESolver::gseid(const double *a, double *x, const double *b, const double err, const int n )
{
    double residue, error, rms, lhs, dummy, sum, rmsold;
    int i, j, check;

     //========= Check for diagonal dominance ========================//
    
    for(int i=0; i<n; i++){
        sum = 0;
        for(int k=0; k<n; k++){
            sum = sum + fabs(a[i*n+k]);
        }
        if(a[i*n+i] < sum-a[i*n+i]){check = 0;}else{check=1;}
    }
    
    if(check ==0){std::cout<<"System is not diagonally dominant \n exiting code"<<std::endl; exit(0);}
    do
    {

        rms = 0;
        for(i =0; i<n; i++)
        {
            lhs = 0;
            for(j=0; j<n; j++)
            {
                lhs = lhs + *(a +i*n +j)*x[j];
            }
            dummy = x[i]; //cout<<"lhs = "<<lhs<<endl;
            x[i] = (b[i] - lhs)/(*(a + i*n +i)) + x[i]; //cout<<"x = "<<x[i]<<endl;
            rms = rms + (x[i] - dummy)*(x[i] - dummy)/abs(dummy)*1000;
        }
        rms = rms/n;
        if(rmsold > rms) {std::cout<<"Divergence detected \t rms value = \t" << rms <<std::endl;}
        rmsold = rms;
    }while(rms > 0.0001);
}

void ODESolver::write()
{
    std::ofstream results("bin/results.dat", std::ios::out);
    for(int k=0; k<n/3; k++)
    {
        results << 0.01*k<<"\t" << y[k*3+0] <<"\t" << y[k*3+1] <<"\t" << y[k*3+2] << "\n";
        //results << k*0.01 << "\t" << y[k] << "\n";
    }
}
void ODESolver::console()
{
    std::cout<<"TimeStep = \t " << tStep << std::endl;
}
