/* 
 * File:   ShockTube.cpp
 * Author: Nikhil
 * 
 * Created on May 16, 2015, 4:30 PM
 */

#include "../include/ShockTube.hpp"
#include <cmath>

/*ShockTube::ShockTube()
{
}*/

/*ShockTube::ShockTube(const ShockTube& orig)
{
}*/

ShockTube::~ShockTube()
{
}

ShockTube::ShockTube(int nodes, double dX):ODESolver(nodes*3)
{
		dx = dX;
}

/* void ShockTube::fcn(double *fValue, const double *y, const double time)
{
	double w1, w2, w3, ww1, ww2, ww3, we1, we2, we3;
	for(int i=1; i<(n/3-1); i++)
	{
		w1 = *(y + i*3 + 0); ww1 = *(y + (i-1)*3 + 0); we1 = *(y + (i+1)*3 + 0);
		w2 = *(y + i*3 + 1); ww2 = *(y + (i-1)*3 + 1); we2 = *(y + (i+1)*3 + 1);
		w3 = *(y + i*3 + 2); ww3 = *(y + (i-1)*3 + 2); we3 = *(y + (i+1)*3 + 2);
		
		fValue[i*3 + 0] = -1*(calcFlux(w1, w2, w3, 1, 1) - calcFlux(ww1, ww2, ww3, 1, 1) + calcFlux(we1, we2, we3, 1, -1) - calcFlux(w1, w2, w3, 1, -1))/dx;
		fValue[i*3 + 1] = -1*(calcFlux(w1, w2, w3, 2, 1) - calcFlux(ww1, ww2, ww3, 2, 1) + calcFlux(we1, we2, we3, 2, -1) - calcFlux(w1, w2, w3, 2, -1))/dx;
		fValue[i*3 + 2] = -1*(calcFlux(w1, w2, w3, 3, 1) - calcFlux(ww1, ww2, ww3, 3, 1) + calcFlux(we1, we2, we3, 3, -1) - calcFlux(w1, w2, w3, 3, -1))/dx;
		
	}
	
}*/

double ShockTube::fcn(const double *y, const double time, const int varIndex)
{
	double w1, w2, w3, ww1, ww2, ww3, we1, we2, we3;
	int i;
        
        i = varIndex/3;
        if(i==0)
        {
            w1 = y[0]; ww1 = 11.6144; we1 = y[3];
            w2 = y[1]; ww2 = 0;       we2 = y[4];
            w3 = y[2]; ww3 = 2500000; we3 = y[5];
        }
        else if(i==n/3-1)
        {
            w1 = y[n-3]; ww1 = y[n-6]; we1 = 1.16144;
            w2 = y[n-2]; ww2 = y[n-5]; we2 = 0;
            w3 = y[n-1]; ww3 = y[n-4]; we3 = 250000;
        }
        else
        {
            w1 = *(y + i*3 + 0); ww1 = *(y + (i-1)*3 + 0); we1 = *(y + (i+1)*3 + 0);
            w2 = *(y + i*3 + 1); ww2 = *(y + (i-1)*3 + 1); we2 = *(y + (i+1)*3 + 1);
            w3 = *(y + i*3 + 2); ww3 = *(y + (i-1)*3 + 2); we3 = *(y + (i+1)*3 + 2);
        }
		
        switch(varIndex%3)
        {
        case 0:
		return -1*(calcFlux(w1, w2, w3, 1, 1) - calcFlux(ww1, ww2, ww3, 1, 1) + calcFlux(we1, we2, we3, 1, -1) - calcFlux(w1, w2, w3, 1, -1))/dx;
                break;
        case 1:
		return -1*(calcFlux(w1, w2, w3, 2, 1) - calcFlux(ww1, ww2, ww3, 2, 1) + calcFlux(we1, we2, we3, 2, -1) - calcFlux(w1, w2, w3, 2, -1))/dx;
                break;
        case 2:
		return -1*(calcFlux(w1, w2, w3, 3, 1) - calcFlux(ww1, ww2, ww3, 3, 1) + calcFlux(we1, we2, we3, 3, -1) - calcFlux(w1, w2, w3, 3, -1))/dx;
                break;		
	}
	
}

double ShockTube::calcFlux(double w1, double w2, double w3, int varTag, int flag)
{
		double gamma, gasConstant, temp, sonic;
	
		gamma = 1.4; gasConstant = 287.0;
		temp = calcTemp(w1, w2, w3);
		sonic = sqrt(gamma*gasConstant*temp);


		switch (varTag)
		{
			case 1:
					if(flag > 0) {return w1/(2*gamma)*( 2*(gamma-1)*w2/w1+(w2/w1 + sonic) );}
					else if(flag< 0) {return w1/(2*gamma)*(w2/w1 - sonic);}
					else{std::cout << "Error in flag value flag is neither greater than zero or less than zero" <<std::endl;}
					break;
			case 2: 
					if(flag > 0) {return w1/(2*gamma)*(2*(gamma -1)*w2*w2/(w1*w1) + (w2/w1 + sonic)*(w2/w1 + sonic));}
					else if(flag < 0){return w1/(2*gamma)*( (w2/w1 - sonic)*(w2/w1 - sonic) );}
					else{std::cout << "Error in flag value flag is neither greater than zero or less than zero" <<std::endl;}
					break;
			case 3:
					if(flag > 0) {return w1/(2*gamma)*( (gamma-1)*w2*w2*w2/(w1*w1*w1) + (w2/w1 + sonic)*(w2/w1 + sonic)*(w2/w1 + sonic)/2 + (3-gamma)/(2*(gamma-1))*(w2/w1 + sonic)*sonic*sonic );}
					else if( flag < 0) {return w1/(2*gamma)*( (w2/w1 - sonic)*(w2/w1 - sonic)*(w2/w1 - sonic)/2.0 + (3-gamma)/(2*(gamma-1))*(w2/w1 - sonic)*sonic*sonic  );}
					else{std::cout << "Error in flag value flag is neither greater than zero or less than zero" <<std::endl;}
					break;
			default :
					std::cout << "Error in variable tag in calcFlux" << std::endl;
		}
	
	
}

double ShockTube::calcTemp( double w1, double w2, double w3)
{
		double gasConstant, gamma;
	
		gamma = 1.4; gasConstant = 287.0;
		return (gamma -1)/gasConstant*( w3/w1 - w2*w2/(2.0*w1*w1));
}

