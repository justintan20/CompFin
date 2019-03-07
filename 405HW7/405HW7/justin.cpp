//
//  justin.cpp
//  405HW7
//
//  Created by Justin Tan on 2/26/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#include "justin.h"

double euroPutEFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_x, double k){
    int s0_index = round((log(s_high) - log(s0))/delta_x) - 1;
    int n = (log(s_high) - log(s_low))/delta_x;
    int m = T / delta_t;
    double p_u = delta_t*(sigma*sigma/(2*delta_x*delta_x)+(r - sigma*sigma/2)/(2*delta_x));
    double p_d = delta_t*(sigma*sigma/(2*delta_x*delta_x)-(r - sigma*sigma/2)/(2*delta_x));
    double p_m = 1-delta_t*sigma*sigma/delta_x/delta_x - r*delta_t;

    vector<double> ending_log_prices;
    for(int j = 0; j <= s0_index; j++){
        ending_log_prices.push_back(log(s0) + (s0_index-j) * delta_x);
    }
    for(int j = 0; j < n/2+s0_index; j++){
        ending_log_prices.push_back(log(s0) - (j+1) * delta_x);
    }
    
    vector<double> ending_payoffs;
    for(double num : ending_log_prices){
        double payoff = k-exp(num);
        if(payoff < 0){
            payoff = 0;
        }
        ending_payoffs.push_back(payoff);
    }
    vector<double> here_values;
    
    for(int times = 0; times < m; times++){
        here_values.push_back(0);
        for(int i = 1; i < n; i++){
            double val = p_u*ending_payoffs[i-1]+p_m*ending_payoffs[i]+p_d*ending_payoffs[i+1];
            here_values.push_back(val);
        }
        here_values[0] = here_values[1];
        here_values.push_back(here_values[n-1] - (ending_log_prices[n] - ending_log_prices[n-1]));
        
        ending_payoffs.clear();
        ending_payoffs = here_values;
        here_values.clear();
    }
    
    
    
    
    return ending_payoffs[s0_index];
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

double euroPutBS(double r, double t, double sigma, double s_0, double k){
    double call = euroCallBS(r, t, sigma, s_0, k);
    double k_now = exp(-r*t)*k;
    double put = call + k_now - s_0;
    return put;
}

double euroPutIFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_x, double k){
    int s0_index = round((log(s_high) - log(s0))/delta_x) - 1;
    int n = (log(s_high) - log(s_low))/delta_x;
    int m = T / delta_t;
    double p_u = -0.5*delta_t*(sigma*sigma/(delta_x*delta_x)+(r - sigma*sigma/2)/(delta_x));
    double p_d = -0.5*delta_t*(sigma*sigma/(delta_x*delta_x)-(r - sigma*sigma/2)/(delta_x));
    double p_m = 1+delta_t*sigma*sigma/delta_x/delta_x + r*delta_t;
    
    vector<double> ending_log_prices;
    for(int j = 0; j <= s0_index; j++){
        ending_log_prices.push_back(log(s0) + (s0_index-j) * delta_x);
    }
    for(int j = 0; j < n/2+s0_index; j++){
        ending_log_prices.push_back(log(s0) - (j+1) * delta_x);
    }
    
    vector<double> ending_payoffs;
    for(double num : ending_log_prices){
        double payoff = k-exp(num);
        if(payoff < 0){
            payoff = 0;
        }
        ending_payoffs.push_back(payoff);
    }
    vector<double> here_values;
    vector<vector<double>> A;
    vector<double> temp{1,-1};
    for(int i = 0; i < n-1; i++){
        temp.push_back(0);
    }
    A.push_back(temp);
    temp.clear();
    for(int i = 0; i < n-1; i++){
        int j = 0;
        while(j < n + 1){
            if(j == i){
                temp.push_back(p_u);
                temp.push_back(p_m);
                temp.push_back(p_d);
                j += 3;
            }
            else{
                temp.push_back(0);
                j++;
            }
        }
        A.push_back(temp);
        temp.clear();
    }
    for(int i = 0; i < n-1; i++){
        temp.push_back(0);
    }
    temp.push_back(1);
    temp.push_back(-1);
    A.push_back(temp);
    
    for(int times = 0; times < m; times++){
        vector<double> B;
        B.push_back(0);
        for(int num = 1; num < n; num++){
            B.push_back(ending_payoffs[num]);
        }
        B.push_back(ending_log_prices[ending_log_prices.size()-1]-ending_log_prices[ending_log_prices.size()-2]);
        here_values = solveLinearSystemLU(A, B);
        ending_payoffs.clear();
        ending_payoffs = here_values;
        here_values.clear();
    }

    
    return ending_payoffs[s0_index];
}

vector<double> vectorMultiplication(const vector<double>& v, double times_num){
    vector<double> output;
    for(double i : v){
        output.push_back(i*times_num);
    }
    return output;
}

vector<double> vectorSubtract(const vector<double>& v1, const vector<double>& v2){
    vector<double> output;
    for(int i = 0; i < v1.size(); i++){
        output.push_back(v1[i]-v2[i]);
    }
    return output;
}

vector<double> solveLinearSystemLU(vector<vector<double>>& A, vector<double>& b){
    vector<vector<double>> lu;
    int n = b.size();
    for(int i = 0; i < n; i++){
        vector<double> temp;
        for(int j = 0; j < n; j++){
            temp.push_back(0);
        }
        lu.push_back(temp);
    }
    
    double sum = 0;
    for(int i = 0; i < n; i++){
        for(int j = i; j < n; j++){
            sum = 0;
            for(int k = 0; k < i; k++){
                sum += lu[i][k] * lu[k][j];
            }
            lu[i][j] = A[i][j] - sum;
        }
        for(int j = i+1; j < n; j++){
            sum = 0;
            for(int k = 0; k < i; k++){
                sum += lu[j][k] * lu[k][i];
            }
            lu[j][i] = (1.0 / lu[i][i])*(A[j][i] - sum);
        }
    }
    
    vector<double> y;
    for(int i = 0; i < n; i++){
        sum = 0;
        for(int k = 0; k < i; k++){
            sum += lu[i][k]*y[k];
        }
        y.push_back(b[i]-sum);
    }
    vector<double> x(n);
    for(int i = n-1; i >= 0; i--){
        sum = 0;
        for(int k = i+1; k < n; k++){
            sum += lu[i][k]*x[k];
        }
        x[i] = (1/lu[i][i])*(y[i]-sum);
    }
    return x;
}

double euroPutCNFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_x, double k){
    //set index and range of values, set up grid
    int s0_index = round((log(s_high) - log(s0))/delta_x) - 1;
    int n = (log(s_high) - log(s_low))/delta_x;
    int m = T / delta_t;
    vector<double> ending_log_prices;
    for(int j = 0; j <= s0_index; j++){
        ending_log_prices.push_back(log(s0) + (s0_index-j) * delta_x);
    }
    for(int j = 0; j < n/2+s0_index; j++){
        ending_log_prices.push_back(log(s0) - (j+1) * delta_x);
    }
    //pu, pd, pm
    double p_u = -0.25*delta_t*(sigma*sigma/(delta_x*delta_x)+(r - sigma*sigma/2)/(delta_x));
    double p_d = -0.25*delta_t*(sigma*sigma/(delta_x*delta_x)-(r - sigma*sigma/2)/(delta_x));
    double p_m = 1+delta_t*sigma*sigma/2.0/delta_x/delta_x + r*delta_t/2.0;
    //ending payoffs
    vector<double> ending_payoffs;
    for(double num : ending_log_prices){
        double payoff = k-exp(num);
        if(payoff < 0){
            payoff = 0;
        }
        ending_payoffs.push_back(payoff);
    }
    //here values is holder
    vector<double> here_values;
    //initialize A matrix
    vector<vector<double>> A;
    //first row
    vector<double> temp{1,-1};
    for(int i = 0; i < n-1; i++){
        temp.push_back(0);
    }
    A.push_back(temp);
    temp.clear();
    //middle rows
    for(int i = 0; i < n-1; i++){
        int j = 0;
        while(j < n + 1){
            if(j == i){
                temp.push_back(p_u);
                temp.push_back(p_m);
                temp.push_back(p_d);
                j += 3;
            }
            else{
                temp.push_back(0);
                j++;
            }
        }
        A.push_back(temp);
        temp.clear();
    }
    //last row
    for(int i = 0; i < n-1; i++){
        temp.push_back(0);
    }
    temp.push_back(1);
    temp.push_back(-1);
    A.push_back(temp);
    
    for(int times = 0; times < m; times++){
        vector<double> B;
        B.push_back(0);
        for(int num = 1; num < n; num++){
            double val = -p_u*ending_payoffs[num - 1] - (p_m - 2.0)*ending_payoffs[num] - p_d*ending_payoffs[num + 1];
            B.push_back(val);
        }
        B.push_back(ending_log_prices[ending_log_prices.size()-1]-ending_log_prices[ending_log_prices.size()-2]);
        here_values = solveLinearSystemLU(A, B);
        ending_payoffs.clear();
        ending_payoffs = here_values;
        here_values.clear();
    }
    
    
    return ending_payoffs[s0_index];
}

double amPutEFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_s, double k){
    int s0_index = (s_high - s0)/delta_s;
    int n = (s_high - s_low)/delta_s;
    int m = T / delta_t;
    
    
    vector<double> ending_prices;
    for(int j = 0; j <= n; j++){
        ending_prices.push_back(s_low + (n-j) * delta_s);
    }
   
    
    vector<double> p_u;
    vector<double> p_m;
    vector<double> p_d;
    for(int i = 1; i < n; i++){
        double j = ending_prices[i]/delta_s;
        p_u.push_back(delta_t*(sigma*sigma*j*j/2 + r*j/2));
        p_d.push_back(delta_t*(sigma*sigma*j*j/2 - r*j/2));
        p_m.push_back(1-delta_t*(sigma*sigma*j*j+r));
    }
    
    
    vector<double> ending_payoffs;
    for(double num : ending_prices){
        double payoff = k-num;
        if(payoff < 0){
            payoff = 0;
        }
        ending_payoffs.push_back(payoff);
    }
    vector<double> here_values;

    for(int times = 0; times < m; times++){
        here_values.push_back(0);
        for(int i = 1; i < n; i++){
            double val = p_u[i-1]*ending_payoffs[i-1]+p_m[i-1]*ending_payoffs[i]+p_d[i-1]*ending_payoffs[i+1];
            double exercise = k - ending_prices[i];
            if(val < exercise){
                val = exercise;
            }
            here_values.push_back(val);
        }
        double initial_exercise = k - ending_prices[0];
        if(initial_exercise > here_values[1]){
            here_values[0] = initial_exercise;
        }
        else{
            here_values[0] = here_values[1];
        }
        double last = here_values[n-1] - (ending_prices[n] - ending_prices[n-1]);
        double last_exercise = k - ending_prices[n];
        if(last < last_exercise){
            last = last_exercise;
        }
        here_values.push_back(last);
        ending_payoffs.clear();
        ending_payoffs = here_values;
        here_values.clear();
    }




    return ending_payoffs[s0_index];
}

double amCallEFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_s, double k){
    int s0_index = (s_high - s0)/delta_s;
    int n = (s_high - s_low)/delta_s;
    int m = T / delta_t;
    
    
    vector<double> ending_prices;
    for(int j = 0; j <= n; j++){
        ending_prices.push_back(s_low + (n-j) * delta_s);
    }
    
    
    vector<double> p_u;
    vector<double> p_m;
    vector<double> p_d;
    for(int i = 1; i < n; i++){
        double j = ending_prices[i]/delta_s;
        p_u.push_back(delta_t*(sigma*sigma*j*j/2 + r*j/2));
        p_d.push_back(delta_t*(sigma*sigma*j*j/2 - r*j/2));
        p_m.push_back(1-delta_t*(sigma*sigma*j*j+r));
    }
    
    
    vector<double> ending_payoffs;
    for(double num : ending_prices){
        double payoff = num-k;
        if(payoff < 0){
            payoff = 0;
        }
        ending_payoffs.push_back(payoff);
    }
    vector<double> here_values;
    
    for(int times = 0; times < m; times++){
        here_values.push_back(0);
        for(int i = 1; i < n; i++){
            double val = p_u[i-1]*ending_payoffs[i-1]+p_m[i-1]*ending_payoffs[i]+p_d[i-1]*ending_payoffs[i+1];
            double exercise = ending_prices[i]-k;
            if(val < exercise){
                val = exercise;
            }
            here_values.push_back(val);
        }
        double initial = here_values[1]+(ending_prices[0]-ending_prices[1]);
        double initial_exercise = ending_prices[0]-k;
        if(initial < initial_exercise){
            initial = initial_exercise;
        }
        here_values[0] = initial;
        double last = here_values[n-1];
        double last_exercise = ending_prices[n]-k;
        if(last < last_exercise){
            last = last_exercise;
        }
        here_values.push_back(last);
        ending_payoffs.clear();
        ending_payoffs = here_values;
        here_values.clear();
    }
    return ending_payoffs[s0_index];
}


double amPutIFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_s, double k){
    int s0_index = (s_high - s0)/delta_s;
    int n = (s_high - s_low)/delta_s;
    int m = T / delta_t;
    
    vector<double> ending_prices;
    for(int j = 0; j <= n; j++){
        ending_prices.push_back(s_low + (n-j) * delta_s);
    }
    
    vector<double> ending_payoffs;
    for(double num : ending_prices){
        double payoff = k-num;
        if(payoff < 0){
            payoff = 0;
        }
        ending_payoffs.push_back(payoff);
    }
    vector<double> p_u;
    vector<double> p_m;
    vector<double> p_d;
    for(int i = 1; i < n; i++){
        double j = ending_prices[i]/delta_s;
        p_u.push_back(-0.5*delta_t*(sigma*sigma*j*j + r*j));
        p_d.push_back(-0.5*delta_t*(sigma*sigma*j*j - r*j));
        p_m.push_back(1+delta_t*(sigma*sigma*j*j+r));
    }
    
    vector<double> here_values;
    vector<vector<double>> A;
    vector<double> temp{1,-1};
    for(int i = 0; i < n-1; i++){
        temp.push_back(0);
    }
    A.push_back(temp);
    temp.clear();
    for(int i = 0; i < n-1; i++){
        int j = 0;
        while(j < n + 1){
            if(j == i){
                temp.push_back(p_u[i]);
                temp.push_back(p_m[i]);
                temp.push_back(p_d[i]);
                j += 3;
            }
            else{
                temp.push_back(0);
                j++;
            }
        }
        A.push_back(temp);
        temp.clear();
    }
    for(int i = 0; i < n-1; i++){
        temp.push_back(0);
    }
    temp.push_back(1);
    temp.push_back(-1);
    A.push_back(temp);
    
    for(int times = 0; times < m; times++){
        vector<double> B;
        B.push_back(0);
        for(int num = 1; num < n; num++){
            B.push_back(ending_payoffs[num]);
        }
        B.push_back(ending_prices[ending_prices.size()-1]-ending_prices[ending_prices.size()-2]);
        here_values = solveLinearSystemLU(A, B);
        //compare with exercise value
        for(int i = 0; i < here_values.size(); i++){
            double exercise = k - ending_prices[i];
            if(exercise > here_values[i]){
                here_values[i] = exercise;
            }
        }
        ending_payoffs.clear();
        ending_payoffs = here_values;
        here_values.clear();
    }
    return ending_payoffs[s0_index];
}

double amCallIFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_s, double k){
    int s0_index = (s_high - s0)/delta_s;
    int n = (s_high - s_low)/delta_s;
    int m = T / delta_t;
    
    vector<double> ending_prices;
    for(int j = 0; j <= n; j++){
        ending_prices.push_back(s_low + (n-j) * delta_s);
    }
    
    vector<double> ending_payoffs;
    for(double num : ending_prices){
        double payoff = num-k;
        if(payoff < 0){
            payoff = 0;
        }
        ending_payoffs.push_back(payoff);
    }
    vector<double> p_u;
    vector<double> p_m;
    vector<double> p_d;
    for(int i = 1; i < n; i++){
        double j = ending_prices[i]/delta_s;
        p_u.push_back(-0.5*delta_t*(sigma*sigma*j*j + r*j));
        p_d.push_back(-0.5*delta_t*(sigma*sigma*j*j - r*j));
        p_m.push_back(1+delta_t*(sigma*sigma*j*j+r));
    }
    
    vector<double> here_values;
    vector<vector<double>> A;
    vector<double> temp{1,-1};
    for(int i = 0; i < n-1; i++){
        temp.push_back(0);
    }
    A.push_back(temp);
    temp.clear();
    for(int i = 0; i < n-1; i++){
        int j = 0;
        while(j < n + 1){
            if(j == i){
                temp.push_back(p_u[i]);
                temp.push_back(p_m[i]);
                temp.push_back(p_d[i]);
                j += 3;
            }
            else{
                temp.push_back(0);
                j++;
            }
        }
        A.push_back(temp);
        temp.clear();
    }
    for(int i = 0; i < n-1; i++){
        temp.push_back(0);
    }
    temp.push_back(1);
    temp.push_back(-1);
    A.push_back(temp);
    
    for(int times = 0; times < m; times++){
        vector<double> B;
        B.push_back(ending_prices[0]-ending_prices[1]);
        for(int num = 1; num < n; num++){
            B.push_back(ending_payoffs[num]);
        }
        B.push_back(0);
        here_values = solveLinearSystemLU(A, B);
        //compare with exercise value
        for(int i = 0; i < here_values.size(); i++){
            double exercise = ending_prices[i] - k;
            if(exercise > here_values[i]){
                here_values[i] = exercise;
            }
        }
        ending_payoffs.clear();
        ending_payoffs = here_values;
        here_values.clear();
    }
    return ending_payoffs[s0_index];
}

double amPutCNFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_s, double k){
    //set index and range of values, set up grid
    int s0_index = (s_high - s0)/delta_s;
    int n = (s_high - s_low)/delta_s;
    int m = T / delta_t;
    
    vector<double> ending_prices;
    for(int j = 0; j <= n; j++){
        ending_prices.push_back(s_low + (n-j) * delta_s);
    }
    
    vector<double> ending_payoffs;
    for(double num : ending_prices){
        double payoff = k - num;
        if(payoff < 0){
            payoff = 0;
        }
        ending_payoffs.push_back(payoff);
    }
    
    vector<double> p_u;
    vector<double> p_m;
    vector<double> p_d;
    for(int i = 1; i < n; i++){
        double j = ending_prices[i]/delta_s;
        p_u.push_back(-0.25*delta_t*(sigma*sigma*j*j + r*j));
        p_d.push_back(-0.25*delta_t*(sigma*sigma*j*j - r*j));
        p_m.push_back(1+delta_t*sigma*sigma*j*j/2+r*delta_t/2);
    }
    
    //here values is holder
    vector<double> here_values;
    //initialize A matrix
    vector<vector<double>> A;
    //first row
    vector<double> temp{1,-1};
    for(int i = 0; i < n-1; i++){
        temp.push_back(0);
    }
    A.push_back(temp);
    temp.clear();
    //middle rows
    for(int i = 0; i < n-1; i++){
        int j = 0;
        while(j < n + 1){
            if(j == i){
                temp.push_back(p_u[i]);
                temp.push_back(p_m[i]);
                temp.push_back(p_d[i]);
                j += 3;
            }
            else{
                temp.push_back(0);
                j++;
            }
        }
        A.push_back(temp);
        temp.clear();
    }
    //last row
    for(int i = 0; i < n-1; i++){
        temp.push_back(0);
    }
    temp.push_back(1);
    temp.push_back(-1);
    A.push_back(temp);
    
    for(int times = 0; times < m; times++){
        vector<double> B;
        B.push_back(0);
        for(int num = 1; num < n; num++){
            double val = -p_u[num - 1]*ending_payoffs[num - 1] - (p_m[num-1] - 2.0)*ending_payoffs[num] - p_d[num-1]*ending_payoffs[num + 1];
            B.push_back(val);
        }
        B.push_back(ending_prices[ending_prices.size()-1]-ending_prices[ending_prices.size()-2]);
        here_values = solveLinearSystemLU(A, B);
        for(int i = 0; i < here_values.size(); i++){
            double exercise = k - ending_prices[i];
            if(exercise > here_values[i]){
                here_values[i] = exercise;
            }
        }
        ending_payoffs.clear();
        ending_payoffs = here_values;
        here_values.clear();
    }
    return ending_payoffs[s0_index];
}



double amCallCNFD(double s0, double s_low, double s_high, double r, double sigma, double T, double delta_t, double delta_s, double k){
    //set index and range of values, set up grid
    int s0_index = (s_high - s0)/delta_s;
    int n = (s_high - s_low)/delta_s;
    int m = T / delta_t;
    
    vector<double> ending_prices;
    for(int j = 0; j <= n; j++){
        ending_prices.push_back(s_low + (n-j) * delta_s);
    }
    
    vector<double> ending_payoffs;
    for(double num : ending_prices){
        double payoff = num - k;
        if(payoff < 0){
            payoff = 0;
        }
        ending_payoffs.push_back(payoff);
    }
    
    vector<double> p_u;
    vector<double> p_m;
    vector<double> p_d;
    for(int i = 1; i < n; i++){
        double j = ending_prices[i]/delta_s;
        p_u.push_back(-0.25*delta_t*(sigma*sigma*j*j + r*j));
        p_d.push_back(-0.25*delta_t*(sigma*sigma*j*j - r*j));
        p_m.push_back(1+delta_t*sigma*sigma*j*j/2+r*delta_t/2);
    }
    
    //here values is holder
    vector<double> here_values;
    //initialize A matrix
    vector<vector<double>> A;
    //first row
    vector<double> temp{1,-1};
    for(int i = 0; i < n-1; i++){
        temp.push_back(0);
    }
    A.push_back(temp);
    temp.clear();
    //middle rows
    for(int i = 0; i < n-1; i++){
        int j = 0;
        while(j < n + 1){
            if(j == i){
                temp.push_back(p_u[i]);
                temp.push_back(p_m[i]);
                temp.push_back(p_d[i]);
                j += 3;
            }
            else{
                temp.push_back(0);
                j++;
            }
        }
        A.push_back(temp);
        temp.clear();
    }
    //last row
    for(int i = 0; i < n-1; i++){
        temp.push_back(0);
    }
    temp.push_back(1);
    temp.push_back(-1);
    A.push_back(temp);
    
    for(int times = 0; times < m; times++){
        vector<double> B;
        B.push_back(ending_prices[0]-ending_prices[1]);
        for(int num = 1; num < n; num++){
            double val = -p_u[num - 1]*ending_payoffs[num - 1] - (p_m[num-1] - 2.0)*ending_payoffs[num] - p_d[num-1]*ending_payoffs[num + 1];
            B.push_back(val);
        }
        B.push_back(0);
        here_values = solveLinearSystemLU(A, B);
        for(int i = 0; i < here_values.size(); i++){
            double exercise = ending_prices[i] - k;
            if(exercise > here_values[i]){
                here_values[i] = exercise;
            }
        }
        ending_payoffs.clear();
        ending_payoffs = here_values;
        here_values.clear();
    }
    return ending_payoffs[s0_index];
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
