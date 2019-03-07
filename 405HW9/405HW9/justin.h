//
//  justin.hpp
//  405HW9
//
//  Created by Justin Tan on 3/6/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#pragma once

#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>

using namespace std;

double priceMBSNumerix(unsigned int T_years, double loan_amount, double wac, double r0, double kappa, double r_bar, double sigma, int num_paths);
double cprNumerixCIR(double r, double mort_rate, double pv0, double pv_tMinus1, unsigned int curr_period, unsigned int month);

vector<vector<double>> ratesCIR(double r0, double kappa, double r_bar, double sigma, int num_paths, int num_steps, int T_years);
double ZCBClosedFormCIR(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T);

vector<double> normGenBM(int num, int seed);
double power(double base, double p);
//file output
void toCSV(const vector<vector<double>>& data, const string& name);
