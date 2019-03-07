//
//  main.cpp
//  405HW3
//
//  Created by Justin Tan on 1/27/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#include "justin.h"

//helper function for question 5
double integral_q5(const vector<double>& h1, const vector<double>& h2){
    double sum = 0;
    for(int i = 0; i < h1.size(); i++){
        double value = exp(-h1[i]*h2[i]) *(sin(6*M_PI*h1[i])+ power(cos(2*M_PI*h2[i]), 1/3.0));
        sum += value;
    }
    double result = sum / h1.size();
    return result;
}



int main(int argc, const char * argv[]) {
    int seed = 3000;
    //problem 1
    cout << "Problem 1" << endl;
    double n = 100.0;
    double m = n*n;
    double T = 2.0;
    vector<double> z = normGenBM(n*m*3, seed);
    //Prob
    vector<double> y_2;
    int count5 = 0;
    int zIndex = 0;
    for(int i = 0; i < m; i++){
        double y_t = 3.0/4.0;
        double t = T/n;
        for(int j = 0; j < n; j++){
            y_t = y_t + (2.0/(1.0+t)*y_t + (1+t*t*t)/3.0)*(T/n) + (1+t*t*t)/3.0*sqrt(T/n)*z[zIndex];
            t += T/n;
            zIndex++;
        }
        y_2.push_back(y_t);
        if(y_t > 5){
            count5++;
        }
    }
    cout << "Prob: " << (double)count5 / y_2.size() << endl;
    
    
    //E1
    vector<double> x_2;
    for(int j = 0; j < m; j++){
        double x_t = 1;
        for(int i = 0; i < n; i++){
            x_t = x_t + (0.2 - 0.5*x_t)*(T/n) + 2.0/3.0*sqrt(T/n)*z[zIndex];
            zIndex++;
        }
        x_2.push_back(x_t);
        
    }
    vector<double> x_2_third;
    for(double i : x_2){
        x_2_third.push_back(power(i, 1.0/3.0));
    }
    double expected_x_2_third = calcMean(x_2_third);
    cout << "E1: " << expected_x_2_third << endl;
    
    //E2
    T = 3.0;
    vector<double> y_3;
    for(int i = 0; i < m; i++){
        double y_t = 3.0/4.0;
        double t = T/n;
        for(int j = 0; j < n; j++){
            y_t = y_t + (2.0/(1.0+t)*y_t + (1+t*t*t)/3.0)*(T/n) + (1+t*t*t)/3.0*sqrt(T/n)*z[zIndex];
            t += T/n;
            zIndex++;
        }
        y_3.push_back(y_t);
    }
    cout << "E2: " << calcMean(y_3) << endl;
    
    //E3
    vector<double> x2y2;
    for(int i = 0; i < x_2.size(); i++){
        if(x_2[i] > 1){
            x2y2.push_back(x_2[i]*y_2[i]);
        }
        else{
            x2y2.push_back(0);
        }
    }
    cout << "E3: " << calcMean(x2y2) << endl;
    
    //problem 2
    cout << "\nProblem 2" << endl;
    //E1
    int seed2 = 12345;
    vector<double> wnormal = normGenBM(n*m, seed);
    vector<double> znormal = normGenBM(n*m, seed2);
    vector<double> x3;
    zIndex = 0;
    T = 3.0;
    for(int i = 0; i < m; i++){
        long double x_t = 1.0;
        for(int j = 0; j < n; j++){
            x_t = x_t + 0.25*x_t*(T/n) + 1.0/3.0*x_t*sqrt(T/n)*wnormal[zIndex] -0.75*x_t*sqrt(T/n)*znormal[zIndex];
            zIndex++;
        }
        x_t = 1+x_t;
        x_t = power(x_t, 1.0/3.0);
        x3.push_back(x_t);
    }
    cout << "E1: " << calcMean(x3) << endl;
    
    //E2
    vector<double> w_3 = brownianSimEnd(seed, T, m);
    vector<double> z_3 = brownianSimEnd(seed2, T, m);
    vector<double> y3_p2;
    for(int i = 0; i < m; i++){
        double y_t = exp(-0.08 * T + 1.0/3.0*w_3[i] + 0.75*z_3[i]);
        y_t += 1;
        y_t = power(y_t, 1.0/3.0);
        y3_p2.push_back(y_t);
    }
    cout << "E2: " << calcMean(y3_p2) << endl;
    
    //problem 3
    cout << "\nProblem 3" << endl;
    double k = 20.0;
    double r = 0.04;
    double sigma = 0.25;
    double time = 0.5;
    //a
    cout << "Price of call assuming initial stock price is 15:" << endl;
    double c_vr = euroCallVR(r, time, sigma, 15, k, seed, 1000);
    cout << "Monte Carlo: " << c_vr << endl;
    //b
    double c_bs = euroCallBS(r, time, sigma, 15, k);
    cout << "Black-Scholes: " << c_bs << endl;
    //c
    vector<double> s_0;
    for(int i = 0; i < 11; i++){
        s_0.push_back(i + 15);
    }
    vector<double> deltas;
    vector<double> vegas;
    vector<double> thetas;
    vector<double> rhos;
    vector<double> gammas;
    for(double s : s_0){
        deltas.push_back(delta(r, time, sigma, s, k));
        vegas.push_back(vega(r, time, sigma, s, k));
        thetas.push_back(theta(r, time, sigma, s, k));
        rhos.push_back(rho(r, time, sigma, s, k));
        gammas.push_back(gamma(r, time, sigma, s, k));
    }
    vector<vector<double>> data;
    data.push_back(deltas);
    data.push_back(vegas);
    data.push_back(thetas);
    data.push_back(rhos);
    data.push_back(gammas);
    toCSV(data, "greeks.csv");
    cout << "All computed greeks values are outputted in file \"greeks.csv\". All plots of greeks are in PDF." << endl;
    
    //problem 4
    cout << "\nProblem 4" << endl;
    //use t = 0.5, k = 50
    double rho = -0.6;
    r = 0.03;
    double s0 = 48.0;
    double v_0 = 0.05;
    sigma = 0.42;
    double alpha = 5.8;
    double beta = 0.0625;
    T = 5;
    k = 50.0;
    
    double price1 = hestonReflect(rho, T, k, r, s0, v_0, sigma, alpha, beta, n, m);
    double price2 = hestonPartTrun(rho, T, k, r, s0, v_0, sigma, alpha, beta, n, m);
    double price3 = hestonFullTrun(rho, T, k, r, s0, v_0, sigma, alpha, beta, n, m);
    cout << "C1: " << price1 << "\nC2: " << price2 << "\nC3: " << price3 << endl;
    
    //problem 5
    cout << "\nProblem 5" << endl;
    vector<double> u = randUniformGen(200, seed);
    vector<double> u1;
    vector<double> u2;
    for(int i = 0; i < 100; i++){
        u1.push_back(u[i]);
        u2.push_back(u[i+100]);
    }
    
    vector<vector<double>> sequences;
    sequences.push_back(u1);
    sequences.push_back(u2);
    
    vector<double> halton_2 = halton(100, 2);
    vector<double> halton_7 = halton(100, 7);
    vector<double> halton_4 = halton(100, 4);
    sequences.push_back(halton_2);
    sequences.push_back(halton_7);
    sequences.push_back(halton_4);
    
    toCSV(sequences, "halton.csv");
    cout << "Plots for sequences are in PDF." << endl;
    
    vector<double> h2 = halton(10000,2);
    vector<double> h4 = halton(10000, 4);
    vector<double> h7 = halton(10000, 7);
    vector<double> h5 = halton(10000, 5);
    
    double integral_2_4 = integral_q5(h2, h4);
    double integral_2_7 = integral_q5(h2, h7);
    double integral_5_7 = integral_q5(h5, h7);
    
    cout << "Integral values: " << endl;
    cout << "Base (2,4): " << integral_2_4 << "\nBase (2,7): " << integral_2_7 << "\nBase (5,7): " << integral_5_7 << endl;
    return 0;
}
