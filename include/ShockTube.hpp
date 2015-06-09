/* 
 * File:   ShockTube.hpp
 * Author: Nikhil
 *
 * Created on May 16, 2015, 4:30 PM
 */

#ifndef SHOCKTUBE_HPP
#define	SHOCKTUBE_HPP

#include "ODESolver.hpp"

class ShockTube:public ODESolver {
public:
    //ShockTube();
    ShockTube(int, double);
    ShockTube(const ShockTube& orig);
    virtual ~ShockTube();

private:
                virtual double fcn(const double*, const double, const int );
		double calcFlux(double w1, double w2, double w3, int varTag, int flag);
		double calcTemp( double w1, double w2, double w3);
		
		double dx;
};

#endif	/* SHOCKTUBE_HPP */

