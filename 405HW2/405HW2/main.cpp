//
//  main.cpp
//  405HW2
//
//  Created by Justin Tan on 1/19/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#include <iostream>
#include "justinComp.h"
#include <fstream>

using namespace std;

//outputs double matrix to csv file
void toCSV(const vector<vector<double>>& data, const string& name){
    ofstream output(name);
    for(auto i : data){
        for(int j = 0; j < i.size() - 1; j++){
            output << i[j] << ",";
        }
        output << i[i.size()-1] << endl;
    }
    output.close();
}

//calculates variance of vector of doubles
double calcVar(const vector<double>& v){
    double size = v.size();
    double mean = calcMean(v);
    double sum = 0.0;
    for(double i : v){
        sum += pow(i-mean,2);
    }
    double var = sum / size;
    return var;
}

//calculates covariance of two vectors of doubles
double calcCov(const vector<double>& a, const vector<double>& b){
    double size = a.size();
    double meanA = calcMean(a);
    double meanB = calcMean(b);
    double sum = 0.0;
    for(int i = 0; i < size; i++){
        double num = (a[i] - meanA)*(b[i] - meanB);
        sum += num;
    }
    double cov = sum / size;
    return cov;
}

//creates bivariate normal distributions with certain given covariance (a) and returns rho
double bivariateNormalRho(int seed, int n, double a){
    double rho = -1;
    vector<double> u = randUniformGen(n*2, seed);
    vector<double> z = normGenBM(u);
    vector<double> x;
    vector<double> y;
    for(int i = 0; i < z.size(); i++){
        x.push_back(z[i]);
        y.push_back(a*z[i] + sqrt(1-a*a)*z[i+1]);
    }
    double varX = calcVar(x);
    double varY = calcVar(y);
    double cov = calcCov(x, y);
    double size = n;
    rho = (cov * size / (size-1))/(sqrt(varX * size / (size - 1))*sqrt(varY * size / (size - 1)));
    return rho;
}

//creates bivariate normal distributions with certain given covariance (a)
void bivariateNormal(int seed, int n, double a, vector<double>& x, vector<double>& y){
    vector<double> u = randUniformGen(n*2, seed);
    vector<double> z = normGenBM(u);
    for(int i = 0; i < z.size(); i++){
        x.push_back(z[i]);
        y.push_back(a*z[i] + sqrt(1-a*a)*z[i+1]);
    }
}

//simulates ending value of brownian motion at given time t
vector<double> brownianSimEnd(int seed, double t, int n){
    vector<double> u = randUniformGen(n, seed);
    vector<double> z = normGenBM(u);
    vector<double> brownian;
    for(double i : z){
        brownian.push_back(sqrt(t)*i);
    }
    return brownian;
}

//simulates european call option payoff at final time (t)
vector<double> euroCall(double r, double t, double sigma, double s_0, double k){
    vector<double> result;
    vector<double> w = brownianSimEnd(400, t, 10000);
    vector<double> s;
    for(double i : w){
        double pay = s_0*exp((r-sigma*sigma/2.0)*t+sigma*i) - k;
        if(pay > 0){
            s.push_back(pay);
        }
        else{
            s.push_back(0);
        }
    }
    for(double j : s){
        result.push_back(j * exp(-1*r*t));
    }
    return result;
}

//simulates stock prices at specified time period (t)
vector<double> stockPrice(double r, double t, double sigma, double s_0, int n){
    vector<double> w = brownianSimEnd(400, t, n);
    vector<double> s;
    for(double i : w){
        double price = s_0*exp((r-sigma*sigma/2.0)*t+sigma*i);
        s.push_back(price);
    }
    return s;
}

//simulates num_paths amount of brownian motion paths till time (t)
vector<vector<double>> brownianPaths(int num_paths, int num_divisions, double t){
    double delta = t / num_divisions;
    vector<vector<double>> result;
    vector<double> u = randUniformGen(num_paths*num_divisions, 2000);
    vector<double> z = normGenBM(u);
    int index = 0;
    for(int i = 0; i < num_paths; i++){
        vector<double> vnow;
        vnow.push_back(sqrt(delta)*z[index]);
        index++;
        for(int j = 1; j < num_divisions; j++){
            vnow.push_back(vnow[j-1]+sqrt(delta) * z[index]);
            index++;
        }
        result.push_back(vnow);
    }
    return result;
}

int main(int argc, const char * argv[]) {
    //question 1
    cout << "Problem 1" << endl;
    int seed = 123456;
    int n = 1000;
    double a = -0.7;
    double rho = bivariateNormalRho(seed, n, a);
    cout << "Rho = " << rho << "\n" << endl;
    
    //question 2
    cout << "Problem 2" << endl;
    vector<double> x2;
    vector<double> y2;
    bivariateNormal(seed, 5000, 0.6, x2, y2);
    vector<double> e;
    for(int i = 0; i < x2.size(); i++){
        double value = pow(x2[i],3) + sin(y2[i]) + x2[i]*x2[i]*y2[i];
        if(value > 0){
            e.push_back(value);
        }
        else{
            e.push_back(0);
        }
    }
    double expected = calcMean(e);
    cout << "Expected value: " << expected << endl;
    
    //question 3
    cout << "\nProblem 3" << endl;
    vector<double> wa1 = brownianSimEnd(seed, 5, 20000);
    vector<double> a1;
    for(double i : wa1){
        a1.push_back(i*i+sin(i));
    }
    double ea1 = calcMean(a1);
    cout << "Ea1: " << ea1 << "\t";
    
    vector<double> wa2 = brownianSimEnd(seed, 0.5, 20000);
    vector<double> a2;
    for(double i : wa2){
        a2.push_back(exp(0.5/2)*cos(i));
    }
    double ea2 = calcMean(a2);
    cout << "Ea2: " << ea2 << "\t";
    
    vector<double> wa3 = brownianSimEnd(seed, 3.2, 20000);
    vector<double> a3;
    for(double i : wa3){
        a3.push_back(exp(3.2/2)*cos(i));
    }
    double ea3 = calcMean(a3);
    
    cout << "Ea3: " << ea3 << "\t";
    vector<double> wa4 = brownianSimEnd(seed, 6.5, 20000);
    vector<double> a4;
    for(double i : wa4){
        a4.push_back(exp(6.5/2)*cos(i));
    }
    double ea4 = calcMean(a4);
    cout << "Ea4: " << ea4 << endl;
    double deva1 = calcVar(a1);
    double deva2 = calcVar(a2);
    double deva3 = calcVar(a3);
    double deva4 = calcVar(a4);
    cout << "Variances: " << deva1 << " " << deva2 << " " << deva3 << " " << deva4 << endl;
    
    
    vector<double> ucontrol = randUniformGen(20000, seed);
    vector<double> ucontrolb;
    for(double u : ucontrol){
        ucontrolb.push_back(u*M_PI);
    }

    vector<double> controlb1;
    for(double z : ucontrolb){
        controlb1.push_back(sin(z));
    }
    vector<double> tb1;
    double gammab1 = calcCov(controlb1, a1) / calcVar(controlb1);
    for(int i = 0; i < a1.size(); i++){
        tb1.push_back(a1[i] - gammab1*(controlb1[i]-2/M_PI));
    }
    cout << "Eb1: " << calcMean(tb1) << "\t";
    
    
    vector<double> controlb2;
    for(double u : ucontrolb){
        controlb2.push_back(cos(u));
    }
    vector<double> tb2;
    double gammab2 = calcCov(controlb2, a2) / calcVar(controlb2);
    for(int i = 0; i < a2.size(); i++){
        tb2.push_back(a2[i] - gammab2*controlb2[i]);
    }
    cout << "Eb2: " << calcMean(tb2) << "\t";
    
    vector<double> tb3;
    double gammab3 = calcCov(controlb2, a3) / calcVar(controlb2);
    for(int i = 0; i < a3.size(); i++){
        tb3.push_back(a3[i] - gammab3*controlb2[i]);
    }
    cout << "Eb3: " << calcMean(tb3) << "\t";
    
    vector<double> tb4;
    double gammab4 = calcCov(controlb2, a4) / calcVar(controlb2);
    for(int i = 0; i < a4.size(); i++){
        tb4.push_back(a4[i] - gammab4*controlb2[i]);
    }
    cout << "Eb4: " << calcMean(tb4) << endl;
    cout <<"Variances: " << calcVar(tb1) << " " << calcVar(tb2) << " " << calcVar(tb3) << " "<< calcVar(tb4) << endl;
    
    //question 4
    cout << "\nProblem 4" << endl;
    vector<double> s = euroCall(0.04, 5, 0.2, 88, 100);
    cout << "Price of call from Monte Carlo Simulation: " << calcMean(s) << endl;
    cout << "Variance of results: " << calcVar(s) << endl;
    vector<double> controlw = brownianSimEnd(400, 5, 10000);
    vector<double> tcall;
    double gammaCall = calcCov(s, controlw) / calcVar(controlw);
    for(int i = 0; i < controlw.size(); i++){
        tcall.push_back(s[i]-gammaCall*controlw[i]);
    }
    cout << "Price of call after variance reduction (control variates): " << calcMean(tcall) << endl;
    cout << "Variance of results: " << calcVar(tcall) << endl;
    
    //question 5
    cout << "\nProblem 5" << endl;
    cout << "All simulated data in this problem are exported as .csv files and all plots are generated from R code, details in PDF.\n" << endl;
    vector<vector<double>> s_n;
    for(int i = 0; i < 10; i++){
        s_n.push_back(stockPrice(0.04, i+1, 0.18, 88, 1000));
    }
    toCSV(s_n, "s_n.csv");
    
    vector<vector<double>> wPaths = brownianPaths(6, 1000, 10);
    vector<vector<double>> sPaths;
    for(auto i : wPaths){
        vector<double> here;
        double t = 10.0/1000.0;
        for(auto j : i){
            here.push_back(88*exp(0.18*j+(0.04-0.18*0.18/2)*t));
            t += (10.0/1000.0);
        }
        sPaths.push_back(here);
    }
    toCSV(sPaths, "s_paths.csv");
    
    vector<vector<double>> s_n2;
    for(int i = 0; i < 10; i++){
        s_n2.push_back(stockPrice(0.04, i+1, 0.35, 88, 1000));
    }
    toCSV(s_n2, "s_n2.csv");
    
    vector<vector<double>> wPaths2 = brownianPaths(6, 1000, 10);
    vector<vector<double>> sPaths2;
    for(auto i : wPaths2){
        vector<double> here;
        double t = 10.0/1000.0;
        for(auto j : i){
            here.push_back(88*exp(0.35*j+(0.04-0.35*0.35/2)*t));
            t += (10.0/1000.0);
        }
        sPaths2.push_back(here);
    }
    toCSV(sPaths2, "s_paths2.csv");
    
    //question 6
    cout << "Problem 6" << endl;
    double division = 1.0/1000.0;
    double sum = 0;
    while (division <= 1){
        sum += sqrt(1-division * division);
        division += (1.0/1000.0);
    }
    double integralEuler = 4*sum / 1000;
    cout << "Estimate (Euler's Method): " << integralEuler << endl;
    
    vector<double> integralU = randUniformGen(1000, seed);
    double sumMC = 0;
    for (double i  : integralU) {
        sumMC += (4*sqrt(1-i * i));
    }
    double integralMC = sumMC / 1000;
    cout << "Estimate (Monte Carlo): " << integralMC << endl;
    
    double alpha = 0.74;
    vector<double> importance;
    for(double i : integralU){
        importance.push_back(4* sqrt(1.0-i*i)/((1.0-alpha*i*i)/(1.0-alpha/3.0)));
    }
    double integralImportance = calcMean(importance);
    cout << "Estimate (Importance Sampling): " << integralImportance << endl;
    return 0;
}
