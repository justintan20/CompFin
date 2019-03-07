//
//  justin.h
//  405HW5
//
//  Created by Justin Tan on 2/14/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#pragma once

#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

//least-square monte carlo for put option
double LSMC(double s0, double k, double r, double sigma, double t, int m, int n, int num_func, const string& type);
//forward start options
double forwardStartEuroPut(double s0, double sigma, double r, double t, double T, int n, int m);
double forwardStartAmericanPut(double s0, double sigma, double r, double t, double t_ex, int m, int n, int num_func, const string& type);

//simulate stock paths
vector<vector<double>> stockPaths(int m, int n, double t, double s0, double sigma, double r);

//math helpers
vector<double> normGenBM(int num);
vector<double> solveLinearSystem(vector<vector<double>>& A, vector<double>& b);
vector<double> vectorMultiplication(const vector<double>& v, double times_num);
vector<double> vectorSubtract(const vector<double>& v1, const vector<double>& v2);
double dotProd(const vector<double>& a, const vector<double>& b);

//update
vector<double> solveLinearSystemLU(vector<vector<double>>& A, vector<double>& b);


//file output
void toCSV(const vector<vector<double>>& data, const string& name);
