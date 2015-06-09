/* 
 * File:   heat.hpp
 * Author: Nikhil
 *
 * Created on May 18, 2015, 1:16 PM
 */

#ifndef HEAT_HPP
#define	HEAT_HPP
#include "ODESolver.hpp"
class heat:public ODESolver
{
public:
    heat(int, double);

private:
    double dx;
    virtual void fcn(double*, const double*, const double );

};

#endif	/* HEAT_HPP */

