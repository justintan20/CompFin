//
//  main.cpp
//  405HW8
//
//  Created by Justin Tan on 3/3/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#include "justin.h"

int main(int argc, const char * argv[]) {
    //problem 1
    cout << "Problem 1" << endl;
    double r0 = 0.05;
    double sigma = 0.18;
    double kappa = 0.82;
    double r_bar = 0.05;
    //a
    double faceValue = 1000;
    double t = 0;
    double T = 0.5;
    int steps = T*365;
    int paths = 10000;
    double price1a_mc = ZCBVasicek(r0, r_bar, sigma, kappa, faceValue, t, T, steps, paths);
//    double price1a_closedform = ZCBClosedFormVasicek(r0, r_bar, sigma, kappa, faceValue, t, T);
    cout << "(a) Price of pure discount bond: " << price1a_mc << endl;
//    cout << price1a_closedform << endl;
    //b
    T = 4;
    steps = 100;
    paths = 10000;
    vector<double> T_i{0.5,1,1.5,2,2.5,3,3.5,4};
    vector<double> c_i{30,30,30,30,30,30,30,1030};
    double price1b_mc = couponBondVasicek(r0, r_bar, sigma, kappa, faceValue, t, T, c_i,T_i, steps, paths);
    cout << "(b) Price of coupon bond: " << price1b_mc << endl;
    //c
    T = 3.0/12.0;
    double S = 0.5;
    double k = 980;
    double price1c_mc = callZCBVasicek(r0, r_bar, sigma, kappa, faceValue, t, T, S, k, steps, paths);
    cout << "(c) Price of call on pure discount bond: " << price1c_mc << endl;
//    double price1c_closed = callZCBClosedFormVasicek(r0, r_bar, sigma, kappa, faceValue, t, T, S, k);
//    cout << price1c_closed << endl;
    
    //d
    S = 4;
    steps = 128;
    paths = 1000;
    for(int i = 0; i < T_i.size(); i++){
        T_i[i] -= T;
    }

    double price1d = callCouponBondVasicek(r0, r_bar, sigma, kappa, faceValue, t, T, S, k, c_i, T_i, steps, paths);
    cout << "(d) Price of call on coupon bond: " << price1d << endl;
    
    
    //p2
    cout << "\nProblem 2" << endl;
    kappa = 0.92;
    r_bar = 0.055;
    S = 1;
    steps = 100;
    paths = 1000;
    T = 0.5;
    
//    double price2a_mc = ZCBCIR(r0, r_bar, sigma, kappa, faceValue, t, S, steps, paths);
//    cout << "(a) Price of pure discount bond: " << price2a_mc << endl;
//    double price2a_closed = ZCBClosedFormCIR(r0, r_bar, sigma, kappa, faceValue, t, S);
//    cout << price2a_closed << endl;
    
    double price2a_mc = callZCBCIR2(r0, r_bar, sigma, kappa, faceValue, t, T, S, k, steps, paths);
    cout << "Price of call on pure discount bond (Monte Carlo): " << price2a_mc << endl;
    paths = 10000;
    double price2a_ex = callZCBCIR(r0, r_bar, sigma, kappa, faceValue, t, T, S, k, steps, paths);
    cout << "Price of call on pure discount bond (explicit formula): " << price2a_ex << endl;
    
    //p3
    cout << "\nProblem 3" << endl;
    r0 = 0.03;
    double x0 = 0;
    double y0 = 0;
    double phi0 = 0.03;
    double rho = 0.7;
    double a = 0.1;
    double b = 0.3;
    sigma = 0.03;
    double eta = 0.08;
    T = 0.5;
    k = 985;
    S = 1;
    steps = 100;
    paths = 1000;
//    double price3a = ZCBG2(r0, sigma, faceValue, t, S, a, b, x0, y0, phi0, eta, rho, steps, paths);
//    cout << price3a << endl;
//    double price3a_closed = ZCBClosedFormG2(r0, sigma, faceValue, t, S, a, b, x0, y0, phi0, eta, rho);
//    cout << price3a_closed << endl;
    double price3 = putZCBG2(r0, sigma, faceValue, t, T, S, k, a, b, x0, y0, phi0, eta, rho, steps, paths);
    cout << "Price of put on pure discount bond: " << price3 << endl;
    double price3_closed = putZCBClosedFormG2(r0, sigma, faceValue, t, T, S, k, a, b, x0, y0, phi0, eta, rho);
    cout << price3_closed << endl;
    return 0;
}
