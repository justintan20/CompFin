//
//  justin.hpp
//  405HW3
//
//  Created by Justin Tan on 1/27/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#pragma once

#include <stdio.h>
#include <vector>
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;

//simulation
vector<double> randUniformGen(int num, int seed);
vector<double> normGenBM(int num, int seed);
vector<double> brownianSimEnd(int seed, double t, int n);
void bivariateNormal(int seed, int n, double a, vector<double>& x, vector<double>& y);
double euroCallVR(double r, double t, double sigma, double s_0, double k, int seed, double n);
double euroCallBS(double r, double t, double sigma, double s_0, double k);
//halton sequence
vector<double> halton(int k, int m);
//greeks
double delta(double r, double t, double sigma, double s_0, double k);
double vega(double r, double t, double sigma, double s_0, double k);
double theta(double r, double t, double sigma, double s_0, double k);
double rho(double r, double t, double sigma, double s_0, double k);
double gamma(double r, double t, double sigma, double s_0, double k);
//heston
double hestonReflect(double rho, double t, double k, double r, double s_0, double v_0, double sigma, double alpha, double beta, int n, int m);
double hestonPartTrun(double rho, double t, double k, double r, double s_0, double v_0, double sigma, double alpha, double beta, int n, int m);
double hestonFullTrun(double rho, double t, double k, double r, double s_0, double v_0, double sigma, double alpha, double beta, int n, int m);
//stat helpers
double calcMean(const vector<double>& v);
double calcStdDev(const vector<double>& v);
double calcVar(const vector<double>& v);
double calcCov(const vector<double>& a, const vector<double>& b);
double cdfNorm(double z);
//math helpers
double dotProd(const vector<double>& a, const vector<double>& b);
double power(double base, double p);
//file output
void toCSV(const vector<vector<double>>& data, const string& name);
