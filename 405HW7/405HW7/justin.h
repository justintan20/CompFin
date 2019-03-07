//
//  justin.hpp
//  405HW7
//
//  Created by Justin Tan on 2/26/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#pragma once

#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

using namespace std;

double euroPutEFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_x, double k);
double euroPutIFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_x, double k);
double euroPutCNFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_x, double k);

double amPutEFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_s, double k);
double amPutIFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_s, double k);
double amPutCNFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_s, double k);

double amCallEFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_s, double k);
double amCallIFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_s, double k);
double amCallCNFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_s, double k);

double euroCallBS(double r, double t, double sigma, double s_0, double k);
double euroPutBS(double r, double t, double sigma, double s_0, double k);

//math helpers
double cdfNorm(double z);
double power(double base, double p);
vector<double> vectorMultiplication(const vector<double>& v, double times_num);
vector<double> vectorSubtract(const vector<double>& v1, const vector<double>& v2);
vector<double> solveLinearSystemLU(vector<vector<double>>& A, vector<double>& b);

//file output
void toCSV(const vector<vector<double>>& data, const string& name);
