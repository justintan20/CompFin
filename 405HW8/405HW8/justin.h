//
//  justin.hpp
//  405HW8
//
//  Created by Justin Tan on 3/3/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#pragma one

#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

double ZCBVasicek(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T, int steps, int paths);
double ZCBClosedFormVasicek(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T);

double couponBondVasicek(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T, const vector<double>& c_i, const vector<double>& T_i, int steps, int paths);

double callZCBVasicek(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T, double S, double k, int steps, int paths);
double callZCBClosedFormVasicek(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T, double S, double k);

double callCouponBondVasicek(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T, double S, double k, const vector<double>& c_i, const vector<double>& T_i, int steps, int paths);

double ZCBCIR(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T, int steps, int paths);
double ZCBClosedFormCIR(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T);

double callZCBCIR(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T, double S, double k, int steps, int paths);
double callZCBCIR2(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T, double S, double k, int steps, int paths);

double ZCBG2(double r0, double sigma, double faceValue, double t, double T, double a, double b, double x0, double y0, double phi0, double eta, double rho, int steps, int paths);
double ZCBClosedFormG2(double r0, double sigma, double faceValue, double t, double T, double a, double b, double x0, double y0, double phi0, double eta, double rho);

double putZCBG2(double r0, double sigma, double faceValue, double t, double T, double S, double k, double a, double b, double x0, double y0, double phi0, double eta, double rho, int steps, int paths);
double putZCBClosedFormG2(double r0, double sigma, double faceValue, double t, double T, double S, double k, double a, double b, double x0, double y0, double phi0, double eta, double rho);

//math helpers
double cdfNorm(double z);
double power(double base, double p);
vector<double> normGenBM(int num, int seed);
