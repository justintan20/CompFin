//
//  justin.hpp
//  405HW6
//
//  Created by Justin Tan on 2/20/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#pragma once

#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

using namespace std;

vector<vector<double>> stockPaths(int m, int n, double t, double s0, double sigma, double r);
vector<vector<double>> stockPathsJumps(int m, int n, double t, double s0, double sigma, double r, double lambda, double gamma);

vector<double> normGenBM(int num);
vector<double> expGen(double lambda, int num, int seed);

double fixedStrikeLookbackCall(double r, double s0, double k, double sigma, double t, int m, int n);
double fixedStrikeLookbackPut(double r, double s0, double k, double sigma, double t, int m, int n);

void Proj6_2function(double lambda1, double lambda2, double T, double& value, double& prob, double& Et);

double vec_max(const vector<double>& v);
double vec_min(const vector<double>& v);
double power(double base, double p);
double calcMean(const vector<double>& v);
double calcVar(const vector<double>& v);

vector<double> poisson(double lambda, int num, int seed);

//file output
void toCSV(const vector<vector<double>>& data, const string& name);
