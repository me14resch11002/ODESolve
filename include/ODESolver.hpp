/* 
 * File:   ODESolver.hpp
 * Author: Nikhil
 *
 * Created on May 16, 2015, 3:58 PM
 */

#ifndef ODESOLVER_HPP
#define	ODESOLVER_HPP


#include <iostream>
#include <fstream>

class ODESolver {
public:
    
    ODESolver(int );
    void initialize(double *, double *, double *);
    void integrate(const int , const double);
    //void integrate(const double , const double);
    void write();

    
    
protected:
    void jacobian(double*, double *, const double);
    virtual double fcn(const double *, const double, const int) =0;
    void gseid(const double *a, double *x, const double *b, const double err, const int n );
    void console();
		
    int n;
    double *y, *yd, *resifcn, *dummyp;
    int tStep;
    double dummy, *jacob, *yp, *ym, time, fp, fm, *a, *x, *b;
    

private:

};

#endif	/* ODESOLVER_HPP */

