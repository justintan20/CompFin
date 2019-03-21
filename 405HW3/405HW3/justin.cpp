//
//  justin.cpp
//  405HW3
//
//  Created by Justin Tan on 1/27/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#include "justin.h"

//generates normal values
vector<double> normGenBM(int num, int seed){
    int a = pow(7,5);
    int m = pow(2,31) - 1;
    int b = 0;
    vector<double> uniform;
    long int xNow = seed;
    long int xNew = 0;
    for(int i = 0; i < num; i++){
        xNew = (a*xNow + b)%m;
        uniform.push_back(xNew / (double)m);
        xNow = xNew;
    }
    vector<double> result;
    int i = 0;
    while(i < uniform.size()){
        result.push_back(sqrt(-2*log(uniform[i]))*cos(2*M_PI*uniform[i+1]));
        result.push_back(sqrt(-2*log(uniform[i]))*sin(2*M_PI*uniform[i+1]));
        i += 2;
    }
    return result;
}

//returns mean of vector of doubles
double calcMean(const vector<double>& v){
    double size = v.size();
    double sum = 0.0;
    for(double i : v){
        sum += i;
    }
    return (sum / size);
}

//returns standard deviation of vector of doubles
double calcStdDev(const vector<double>& v){
    double size = v.size();
    double mean = calcMean(v);
    double sum = 0.0;
    for(double i : v){
        sum += pow(i-mean,2);
    }
    double stdDev = sqrt(sum / size);
    return stdDev;
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

//simulates ending brownian motion value at time T
vector<double> brownianSimEnd(int seed, double t, int n){
    vector<double> z = normGenBM(n, seed);
    vector<double> brownian;
    for(double i : z){
        brownian.push_back(sqrt(t)*i);
    }
    return brownian;
}

//calculates price of european call with variance reduction (antithetic variables)
double euroCallVR(double r, double t, double sigma, double s_0, double k, int seed, double n){
    vector<double> result;
    vector<double> z1;
    vector<double> z2;
    bivariateNormal(seed, n, -1, z1, z2);
    vector<double> s;
    for(int i = 0; i < n; i++){
        double pay1 = s_0*exp((r-sigma*sigma/2.0)*t+sigma*sqrt(t)*z1[i]) - k;
        if(pay1 < 0){
            pay1 = 0;
        }
        double pay2 = s_0*exp((r-sigma*sigma/2.0)*t+sigma*sqrt(t)*z2[i]) - k;
        if(pay2 < 0){
            pay2 = 0;
        }
        double pay = (pay1 + pay2)/2.0;
        s.push_back(pay);
    }
    for(double j : s){
        result.push_back(j * exp(-1.0*r*t));
    }
    return calcMean(result);
}

//creates bivariate normal distributions with certain given covariance (a)
void bivariateNormal(int seed, int n, double a, vector<double>& x, vector<double>& y){
    vector<double> z = normGenBM(n*2, seed);
    for(int i = 0; i < z.size(); i += 2){
        x.push_back(z[i]);
        y.push_back(a*z[i] + sqrt(1-a*a)*z[i+1]);
    }
}

//calculates price of european call by Black-Scholes
double euroCallBS(double r, double t, double sigma, double s_0, double k){
    double d1 = 1.0/(sigma * sqrt(t))*(log(s_0/k)+(r+0.5*sigma*sigma)*t);
    double d2 = d1 - sigma*sqrt(t);
    double c = s_0*cdfNorm(d1) - exp(-1.0*r*t)*k*cdfNorm(d2);
    return c;
}

//numerically estimates cdf of normal
double cdfNorm(double z){
    bool isNeg = false;
    if(z < 0){
        z *= -1;
        isNeg = true;
    }
    double d1 = 0.0498673470;
    double d2 = 0.0211410061;
    double d3 = 0.0032776323;
    double d4 = 0.0000380036;
    double d5 = 0.0000488906;
    double d6 = 0.0000053830;
    double n = 1-0.5*pow((1+d1*z+d2*z*z+d3*pow(z, 3.0)+d4*pow(z, 4.0)+d5*pow(z, 5.0)+d6*pow(z, 6)),-1*16);
    if(isNeg){
        n = 1 - n;
    }
    return n;
}

//greeks
double delta(double r, double t, double sigma, double s_0, double k){
    double epsilon = 0.0001;
    double c1 = euroCallVR(r, t, sigma, s_0 - epsilon, k,4000,10000);
    double c2 = euroCallVR(r, t, sigma, s_0 + epsilon, k,4000,10000);
    double result = (c2 - c1)/2.0/epsilon;
    return result;
}

double vega(double r, double t, double sigma, double s_0, double k){
    double epsilon = 0.0001;
    double c1 = euroCallVR(r, t, sigma - epsilon, s_0, k,4000,10000);
    double c2 = euroCallVR(r, t, sigma + epsilon, s_0, k,4000,10000);
    double result = (c2 - c1)/2.0/epsilon;
    return result;
}

double theta(double r, double t, double sigma, double s_0, double k){
    double epsilon = 0.0001;
    double c1 = euroCallVR(r, t - epsilon, sigma, s_0, k,4000,10000);
    double c2 = euroCallVR(r, t + epsilon, sigma, s_0, k,4000,10000);
    double result = (c2 - c1)/(-1.0)/2.0/epsilon;
    return result;
}

double rho(double r, double t, double sigma, double s_0, double k){
    double epsilon = 0.00001;
    double c1 = euroCallVR(r-epsilon, t, sigma, s_0, k,4000,10000);
    double c2 = euroCallVR(r+epsilon, t, sigma, s_0, k,4000,10000);
    double result = (c2 - c1)/2/epsilon;
    return result;
}

double gamma(double r, double t, double sigma, double s_0, double k){
    double epsilon = 0.0001;
    double delta1 = delta(r, t, sigma, s_0-epsilon, k);
    double delta2 = delta(r, t, sigma, s_0 + epsilon, k);
    double result = (delta2 - delta1)/2/epsilon;
    return result;
}

//heston reflection
double hestonReflect(double rho, double t, double k, double r, double s_0, double v_0, double sigma, double alpha, double beta, int n, int m){
    vector<double> z_s;
    vector<double> z_v;
    bivariateNormal(50000, n*m, rho, z_s, z_v);
    int zIndex = 0;
    vector<double> s_T;
    double delta = t / n;
    for(int i = 0; i < m; i++){
        double s = s_0;
        double v = v_0;
        for(int j = 0; j < n; j++){
            if(v < 0){
                s = s + r*s*delta + sqrt((-1.0)*v)*s*sqrt(delta)*z_s[zIndex];
                v = (-1.0)*v + alpha*(beta + v)*delta+sigma*sqrt((-1.0)*v)*sqrt(delta)*z_v[zIndex];
            }
            else{
                s = s + r*s*delta + sqrt(v)*s*sqrt(delta)*z_s[zIndex];
                v = v + alpha*(beta - v)*delta+sigma*sqrt(v)*sqrt(delta)*z_v[zIndex];
            }
            zIndex++;
        }
        s = s-k;
        if(s<0){
            s = 0;
        }
        s_T.push_back(s);
    }
    return (exp(-1.0*r*t)*calcMean(s_T));
}

//heston partial truncate
double hestonPartTrun(double rho, double t, double k, double r, double s_0, double v_0, double sigma, double alpha, double beta, int n, int m){
    vector<double> z_s;
    vector<double> z_v;
    bivariateNormal(4000, n*m, rho, z_s, z_v);
    int zIndex = 0;
    vector<double> s_T;
    for(int i = 0; i < m; i++){
        double s = s_0;
        double v = v_0;
        double delta = t / n;
        for(int j = 0; j < n; j++){
            if(v < 0){
                s = s + r*s*delta;
                v = v + alpha*(beta - v)*delta;
            }
            else{
                s = s + r*s*delta + sqrt(v)*s*sqrt(delta)*z_s[zIndex];
                v = v + alpha*(beta - v)*delta+sigma*sqrt(v)*sqrt(delta)*z_v[zIndex];
            }
            zIndex++;
        }
        s = s-k;
        if(s<0){
            s = 0;
        }
        s_T.push_back(s);
    }
    return (exp(-1.0*r*t)*calcMean(s_T));
}

//heston fully truncate
double hestonFullTrun(double rho, double t, double k, double r, double s_0, double v_0, double sigma, double alpha, double beta, int n, int m){
    vector<double> z_s;
    vector<double> z_v;
    bivariateNormal(4000, n*m, rho, z_s, z_v);
    int zIndex = 0;
    vector<double> s_T;
    for(int i = 0; i < m; i++){
        double s = s_0;
        double v = v_0;
        double delta = t / n;
        for(int j = 0; j < n; j++){
            if(v < 0){
                s = s + r*s*delta;
                v = v + alpha*beta*delta;
            }
            else{
                s = s + r*s*delta + sqrt(v)*s*sqrt(delta)*z_s[zIndex];
                v = v + alpha*(beta - v)*delta+sigma*sqrt(v)*sqrt(delta)*z_v[zIndex];
            }
            zIndex++;
        }
        s = s-k;
        if(s<0){
            s = 0;
        }
        s_T.push_back(s);
    }
    return (exp(-1.0*r*t)*calcMean(s_T));
}

//generate halton sequence
vector<double> halton(int k, int m){
    int num = (log(k) / log(m)) + 1;
    vector<double> seq;
    vector<double> vetBase;
    vector<double> workVet;
    for(int i = 0; i < k; i++){
        seq.push_back(0);
    }
    for(int i = 0; i < num; i++){
        vetBase.push_back(pow(m, i+1));
        workVet.push_back(0);
    }
    for(int i = 0; i < k; i++){
        int j = 0;
        bool ok = false;
        while(!ok){
            workVet[j] = workVet[j] + 1;
            if(workVet[j] < m){
                ok = true;
            }
            else{
                workVet[j] = 0;
                j++;
            }
        }
        seq[i] = dotProd(workVet, vetBase);
    }
    return seq;
}
//calculates dot product of vectors
double dotProd(const vector<double>& a, const vector<double>& b){
    double sum = 0;
    for(int i = 0; i < a.size(); i++){
        sum += a[i]/b[i];
    }
    return sum;
}

/*
 generates numbers from uniform distribution
 inputs: parameters a, b, m, number of numbers to generate, seed
 output: vector of generated random numbers
 */
vector<double> randUniformGen(int num, int seed){
    int a = pow(7,5);
    int m = pow(2,31) - 1;
    int b = 0;
    vector<double> result;
    long int xNow = seed;
    long int xNew = 0;
    for(int i = 0; i < num; i++){
        xNew = (a*xNow + b)%m;
        result.push_back(xNew / (double)m);
        xNow = xNew;
    }
    return result;
}

//power of number, accounts for negative values
double power(double base, double p){
    double x = base;
    double result;
    if(x < 0){
        x *= -1;
        result = pow(x, p);
        result *= -1;
    }
    else{
        result = pow(x, p);
    }
    return result;
}


