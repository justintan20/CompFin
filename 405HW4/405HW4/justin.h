//
//  justin.hpp
//  405HW4
//
//  Created by Justin Tan on 2/5/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#pragma once

#include <stdio.h>
#include <vector>
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;

//binomial methods
double binomEuroCallA(double r, double sigma, double s, double k, double t, int n);
double binomEuroCallB(double r, double sigma, double s, double k, double t, int n);
double binomEuroCallJR(double r, double sigma, double s, double k, double t, int n);
double binomEuroCallCRR(double r, double sigma, double s, double k, double t, int n);
double binomEuroPutCRR(double r, double sigma, double s, double k, double t, int n);
double binomAmPutCRR(double r, double sigma, double s, double k, double t, int n);
//trinomial methods
double trinomEuroCall(double r, double sigma, double s, double k, double t, int n);
double trinomEuroCallLog(double r, double sigma, double s, double k, double t, int n);
//halton
double euroCall(double r, double t, double sigma, double s_0, double k, int n, int base1, int base2);



//greeks
double delta(double r, double t, double sigma, double s_0, double k);
double vega(double r, double t, double sigma, double s_0, double k);
double theta(double r, double t, double sigma, double s_0, double k);
double rho(double r, double t, double sigma, double s_0, double k);
double gamma(double r, double t, double sigma, double s_0, double k);

//halton sequence
vector<double> halton(int k, int m);
//halton generated normals
vector<double> normGenBMHalton(int num, int base1, int base2);
//sim brownian motion
vector<double> brownianSimEnd(int base1, int base2, double t, int n);

//math helpers
double power(double base, double p);
double calcStdDev(const vector<double>& v);
double calcMean(const vector<double>& v);
double dotProd(const vector<double>& a, const vector<double>& b);
//file output
void toCSV(const vector<vector<double>>& data, const string& name);
