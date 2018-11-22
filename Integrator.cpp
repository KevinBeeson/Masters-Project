// Integrator.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;
// weifjiwef
double potential(vector<double> r)
{
	double gConstant, mSun, V;
	gConstant = 1;
	mSun = 1;
	V = -gConstant*mSun / r[0];
	return V;
}
vector<double> changeVarCP(vector<double> r)
{
	vector<double> rNew(3);
	rNew[0] = pow(r[0] * r[0] + r[1] * r[1] + r[2] * r[2], 0.5);
	rNew[1] = atan(r[1] / r[0]);
	rNew[2] = atan(pow(r[0] * r[0] + r[1] * r[1], 0.5) / r[2]);
	return rNew;
}
vector<double> changeVarPC(vector<double> r)
{
	vector<double> rNew(3);
	rNew[0] = r[0]*sin(r[2])*cos(r[1]);
	rNew[1] = r[0] * sin(r[2])*sin(r[1]);
	rNew[2] = r[0] * cos(r[2]);
	return rNew;
}
vector<double> Dq(vector<double> r, vector<double> p)
{
	double gConstant, mSun,R;
	vector<double> Dv(3);
	gConstant = 1;
	mSun = 1;
	R = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
	Dv[0] = -gConstant*mSun*r[0] / R;
	Dv[1] = -gConstant*mSun*r[1] / R;
	Dv[2] = -gConstant*mSun*r[2] / R;
	return Dv;
}
int main()
{
	vector<double> r(3), p(3), Dv(3);
	double  historyR[700][3], historyP[700][3], step;
	ofstream myfile;
	myfile.open("R.txt");
	myfile << "x \t px \t y \t py \t z \t pz \t t \n";
	step = 100;
	r[0] = 1;
	r[1] = 0;
	r[2] = 0;
	p[0] = 0;
	p[1] = 1;
	p[2] = 0;
	for (int t = 0; t < 700 ;  ++t) {
		//
		
		cout << "here \n ";
		Dv = Dq(r,p);
		for (int i = 0; i < 3; ++i) {
			historyR[t][i] = r[i];
			historyP[t][i] = p[i];
			p[i] = p[i] + 1 / step* Dv[i];
			r[i] = r[i] + 1 / step*p[i];
			// cout << r[i]<< "\t"<< p[i]<<"\n";
			myfile << historyR[t][i] << "\t"<< historyP[t][i] << "\t";
		}

		myfile << t/step << "\n";
	}
	myfile.close();
	system("PAUSE");
    return 0;
}


