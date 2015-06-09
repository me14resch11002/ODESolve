/* 
 * File:   main.cpp
 * Author: Nikhil
 *
 * Created on May 27, 2015, 4:53 AM
 */

#include <cstdlib>
#include "include/ShockTube.hpp"


using namespace std;

/*
 * 
 */
int main(int argc, char** argv)
{
        double gamma, gConstant, rho, uVel, temp, dx;
	double *w, *wd, *resifcn;
	int nodes;
	
	nodes = 1001; dx = 0.01;
	w = new double [3*nodes]; wd = new double [3*nodes]; resifcn = new double [3*nodes];
	
	for(int i=0; i<nodes; i++)
	{
		if(i<nodes/2){rho = 11.6144;}else{rho = 1.16144;}
		w[i*3 + 0] = rho;
		w[i*3 + 1] = 0;
		w[i*3 + 2] = rho*287.0/(1.4 - 1)*300;
		/*wd[i*3 + 0] = 0; wd[i*3 + 1] = 0; wd[i*3 + 2] = 0;
		resifcn[i*3 + 0] = 0;resifcn[i*3 + 1] = 0;resifcn[i*3 + 2] = 0;*/
	}
        //wd[0] = 300/(dx*dx);
	
		ShockTube stube(nodes, dx);
		stube.initialize(w, wd, resifcn);
                stube.integrate(60, 0.0001);
	
                stube.write();
    return 0;
}

