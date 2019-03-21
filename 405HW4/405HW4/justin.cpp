//
//  justin.cpp
//  405HW4
//
//  Created by Justin Tan on 2/5/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#include "justin.h"

double binomEuroCallA(double r, double sigma, double s, double k, double t, int n){
    double delta = t / n;
    double c = 0.5*(exp(-1.0*r*delta)+exp((r+sigma*sigma)*delta));
    double d = c - sqrt(c*c - 1.0);
    double u = 1.0/d;
    double p = (exp(r*delta) - d) / (u - d);

    vector<vector<double>> tree;
    vector<double> endVals;
    for(int i = 0; i < n + 1; i++){
        double val = s*power(u, n-i)*power(d, i) - k;
        if(val < 0){
            val = 0;
        }
        endVals.push_back(val);
    }
    tree.push_back(endVals);
    for(int i = 1; i < n + 1; i++){
        vector<double> valHere;
        for(int j = 0; j < n-i+1; j++){
            valHere.push_back((tree[i-1][j] * p + (1-p) * tree[i-1][j+1])*exp(-1*r*delta));
        }
        tree.push_back(valHere);
    }
    double result = tree[n][0];
    return result;
}

double binomEuroCallB(double r, double sigma, double s, double k, double t, int n){
    double delta = t/n;
    double u = exp(delta*r)*(1.0+sqrt(exp(sigma*sigma*delta)-1.0));
    double d = exp(delta*r)*(1.0-sqrt(exp(sigma*sigma*delta)-1.0));
    double p = 0.5;
    vector<vector<double>> tree;
    vector<double> endVals;
    for(int i = 0; i < n + 1; i++){
        double val = s*power(u, n-i)*power(d, i) - k;
        if(val < 0){
            val = 0;
        }
        endVals.push_back(val);
    }
    tree.push_back(endVals);
    for(int i = 1; i < n + 1; i++){
        vector<double> valHere;
        for(int j = 0; j < n-i+1; j++){
            valHere.push_back((tree[i-1][j] * p + (1-p) * tree[i-1][j+1])*exp(-1*r*delta));
        }
        tree.push_back(valHere);
    }
    double result = tree[n][0];
    return result;
}

double binomEuroCallCRR(double r, double sigma, double s, double k, double t, int n){
    double delta = t/n;
    double d = exp(-1.0*sigma*sqrt(delta));
    double u = exp(1.0*sigma*sqrt(delta));
    double p = 0.5 * (1.0 + (r - sigma * sigma / 2.0)*sqrt(delta)/sigma);
    vector<vector<double>> tree;
    vector<double> endVals;
    for(int i = 0; i < n + 1; i++){
        double val = s*power(u, n-i)*power(d, i) - k;
        if(val < 0){
            val = 0;
        }
        endVals.push_back(val);
    }
    tree.push_back(endVals);
    for(int i = 1; i < n + 1; i++){
        vector<double> valHere;
        for(int j = 0; j < n-i+1; j++){
            valHere.push_back((tree[i-1][j] * p + (1-p) * tree[i-1][j+1])*exp(-1*r*delta));
        }
        tree.push_back(valHere);
    }
    double result = tree[n][0];
    return result;
}

double binomEuroCallJR(double r, double sigma, double s, double k, double t, int n){
    double delta = t/n;
    double u = exp((r - sigma * sigma / 2.0)*delta + sigma * sqrt(delta));
    double d = exp((r - sigma * sigma / 2.0)*delta - sigma * sqrt(delta));
    double p = 0.5;
    vector<vector<double>> tree;
    vector<double> endVals;
    for(int i = 0; i < n + 1; i++){
        double val = s*power(u, n-i)*power(d, i) - k;
        if(val < 0){
            val = 0;
        }
        endVals.push_back(val);
    }
    tree.push_back(endVals);
    for(int i = 1; i < n + 1; i++){
        vector<double> valHere;
        for(int j = 0; j < n-i+1; j++){
            valHere.push_back((tree[i-1][j] * p + (1-p) * tree[i-1][j+1])*exp(-1*r*delta));
        }
        tree.push_back(valHere);
    }
    double result = tree[n][0];
    return result;
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

//returns mean of vector of doubles
double calcMean(const vector<double>& v){
    double size = v.size();
    double sum = 0.0;
    for(double i : v){
        sum += i;
    }
    return (sum / size);
}

//greeks
double delta(double r, double t, double sigma, double s_0, double k){
    double epsilon = 0.0001;
    double c1 = binomEuroCallA(r, sigma, s_0 - epsilon, k, t, 500);
    double c2 = binomEuroCallA(r, sigma, s_0 + epsilon, k, t, 500);
    double result = (c2 - c1)/2/epsilon;
    return result;
}

double vega(double r, double t, double sigma, double s_0, double k){
    double epsilon = 0.0001;
    double c1 = binomEuroCallA(r, sigma - epsilon, s_0, k, t, 500);
    double c2 = binomEuroCallA(r, sigma + epsilon, s_0, k, t, 500);
    double result = (c2 - c1)/2/epsilon;
    return result;
}

double theta(double r, double t, double sigma, double s_0, double k){
    double epsilon = 0.0001;
    double c1 = binomEuroCallA(r, sigma, s_0, k, t - epsilon, 500);
    double c2 = binomEuroCallA(r, sigma, s_0, k, t + epsilon, 500);
    double result = (c2 - c1)/(-1.0)/2.0/epsilon;
    return result;
}

double rho(double r, double t, double sigma, double s_0, double k){
    double epsilon = 0.00001;
    double c1 = binomEuroCallA(r - epsilon, sigma, s_0, k, t, 500);
    double c2 = binomEuroCallA(r + epsilon, sigma, s_0, k, t, 500);
    double result = (c2 - c1)/2/epsilon;
    return result;
}

double gamma(double r, double t, double sigma, double s_0, double k){
    double epsilon = 0.0001;
    double delta1 = delta(r, t, sigma, s_0 - epsilon, k);
    double delta2 = delta(r, t, sigma, s_0 + epsilon, k);
    double result = (delta2 - delta1)/2/epsilon;
    return result;
}

double binomEuroPutCRR(double r, double sigma, double s, double k, double t, int n){
    double delta = t/n;
    double d = exp(-1.0*sigma*sqrt(delta));
    double u = exp(1.0*sigma*sqrt(delta));
    double p = 0.5 * (1.0 + (r - sigma * sigma / 2.0)*sqrt(delta)/sigma);
    vector<vector<double>> tree;
    vector<double> endVals;
    for(int i = 0; i < n + 1; i++){
        double val = k - s*power(u, n-i)*power(d, i);
        if(val < 0){
            val = 0;
        }
        endVals.push_back(val);
    }
    tree.push_back(endVals);
    for(int i = 1; i < n + 1; i++){
        vector<double> valHere;
        for(int j = 0; j < n-i+1; j++){
            valHere.push_back((tree[i-1][j] * p + (1-p) * tree[i-1][j+1])*exp(-1*r*delta));
        }
        tree.push_back(valHere);
    }
    double result = tree[n][0];
    return result;
}

double binomAmPutCRR(double r, double sigma, double s, double k, double t, int n){
    double delta = t/n;
    double d = exp(-1.0*sigma*sqrt(delta));
    double u = exp(1.0*sigma*sqrt(delta));
    double p = 0.5 * (1.0 + (r - sigma * sigma / 2.0)*sqrt(delta)/sigma);
    vector<vector<double>> s_tree;
    for(int i = 0; i < n ; i++){
        vector<double> sHere;
        for (int j = 0; j < i + 2; j++) {
            sHere.push_back(s*power(u, i+1-j)*power(d, j));
        }
        s_tree.push_back(sHere);
    }
    vector<vector<double>> payoff_tree;
    vector<double> endPayoff;
    for(int i = 0; i < n + 1; i++){
        double val = k - s_tree[n-1][i];
        if(val < 0){
            val = 0;
        }
        endPayoff.push_back(val);
    }
    payoff_tree.push_back(endPayoff);
    for(int i = 1; i < n; i++){
        vector<double> payoffHere;
        for(int j = 0; j < n-i+1; j++){
            double cv = (payoff_tree[i-1][j] * p + payoff_tree[i-1][j+1]*(1-p))*exp(-1*r*delta);
            double ev =  k - s_tree[n-1-i][j];
            if(ev < 0){
                ev = 0;
            }
            if(cv >= ev){
                payoffHere.push_back(cv);
            }
            else{
                payoffHere.push_back(ev);
            }
            
        }
        payoff_tree.push_back(payoffHere);
    }
    double start = (payoff_tree[n-1][0] * p + payoff_tree[n-1][1]*(1-p))*exp(-1*r*delta);
    double result = start;
    return result;
}

double trinomEuroCall(double r, double sigma, double s, double k, double t, int n){
    double delta = t / n;
    double d = exp(-1.0*sigma*sqrt(3.0*delta));
    double u = 1.0 / d;
    double p_d = (r*delta*(1.0-u)+power(r*delta,2.0)+sigma*sigma*delta)/(u-d)/(1.0-d);
    double p_u = (r*delta*(1.0-d)+power(r*delta,2.0)+sigma*sigma*delta)/(u-d)/(u-1.0);
    double p_m = 1.0-p_d - p_u;
    vector<vector<double>> tree;
    vector<double> endPayoff;
    for(int i = 0; i < n; i++){
        double val = s*power(u, n - i) - k;
        if(val < 0){
            val = 0;
        }
        endPayoff.push_back(val);
    }
    double samepay = s - k;
    if(samepay < 0){
        samepay = 0;
    }
    endPayoff.push_back(samepay);
    for(int i = 1; i < n+1; i++){
        double val = s*power(d, i) - k;
        if(val < 0){
            val = 0;
        }
        endPayoff.push_back(val);
    }
    tree.push_back(endPayoff);
    for(int i = 1; i < n + 1; i++){
        vector<double> valHere;
        for(int j = 0; j < (2*(n-i) + 1); j++){
            valHere.push_back((tree[i-1][j] * p_u + p_d * tree[i-1][j+2] + p_m * tree[i-1][j+1])*exp(-1*r*delta));
        }
        tree.push_back(valHere);
    }
    double result = tree[n][0];
    return result;
}

double trinomEuroCallLog(double r, double sigma, double s, double k, double t, int n){
    double x_0 = log(s);
    double delta = t / n;
    double x_u = sigma*sqrt(3.0*delta);
    double x_d = -1*sigma*sqrt(3.0*delta);
    double p_d = 0.5*(((sigma*sigma*delta + power(r - sigma*sigma/2.0,2.0)*delta*delta) / (x_u*x_u)) - (r - sigma*sigma/2.0)*delta/x_u);
    double p_u = 0.5*(((sigma*sigma*delta + power(r - sigma*sigma/2.0,2.0)*delta*delta) / (x_u*x_u)) + (r - sigma*sigma/2.0)*delta/x_u);
    double p_m = 1.0-p_d - p_u;
    vector<vector<double>> tree;
    vector<double> endPayoff;
    for(int i = 0; i < n; i++){
        double val = x_0 + x_u*(n-i);
        val = exp(val) - k;
        if(val < 0){
            val = 0;
        }
        endPayoff.push_back(val);
    }
    double samepay = exp(x_0) - k;
    if(samepay < 0){
        samepay = 0;
    }
    endPayoff.push_back(samepay);
    for(int i = 1; i < n + 1; i++){
        double val = x_0 + x_d*i;
        val = exp(val) - k;
        if(val < 0){
            val = 0;
        }
        endPayoff.push_back(val);
    }
    tree.push_back(endPayoff);
    for(int i = 1; i < n + 1; i++){
        vector<double> valHere;
        for(int j = 0; j < (2*(n-i) + 1); j++){
            valHere.push_back((tree[i-1][j] * p_u + p_d * tree[i-1][j+2] + p_m * tree[i-1][j+1])*exp(-1*r*delta));
        }
        tree.push_back(valHere);
    }
    double result = tree[n][0];
    return result;
}

//generate halton sequence
//m is base, k is number of values
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

//generates normal values
vector<double> normGenBMHalton(int num, int base1, int base2){
    vector<double> h1 = halton(num, base1);
    vector<double> h2 = halton(num, base2);
    vector<double> result;
    for(int i = 0; i < h1.size(); i++){
        result.push_back(sqrt(-2*log(h1[i]))*cos(2*M_PI*h2[i]));
        result.push_back(sqrt(-2*log(h1[i]))*sin(2*M_PI*h2[i]));
    }
    return result;
}

//calculates dot product of vectors
double dotProd(const vector<double>& a, const vector<double>& b){
    double sum = 0;
    for(int i = 0; i < a.size(); i++){
        sum += a[i]/b[i];
    }
    return sum;
}

//simulates ending brownian motion value at time T
vector<double> brownianSimEnd(int base1, int base2, double t, int n){
    vector<double> z = normGenBMHalton(n, base1, base2);
    vector<double> brownian;
    for(double i : z){
        brownian.push_back(sqrt(t)*i);
    }
    return brownian;
}

//simulates european call option payoff at final time (t)
double euroCall(double r, double t, double sigma, double s_0, double k, int n, int base1, int base2){
    vector<double> result;
    vector<double> w = brownianSimEnd(base1, base2, t, n);
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
    return calcMean(result);
}
