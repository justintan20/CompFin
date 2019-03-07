//
//  justin.cpp
//  405HW6
//
//  Created by Justin Tan on 2/20/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#include "justin.h"

vector<vector<double>> stockPaths(int m, int n, double t, double s0, double sigma, double r){
    vector<double> z = normGenBM(m*n);
    vector<vector<double>> sPaths;
    double delta = t/n;
    int zIndex = 0;
    for(int i = 0; i < m; i++){
        vector<double> here;
        double s_before = s0;
        here.push_back(s_before);
        for(int j = 0; j < n; j++){
            double s_new = s_before + r*s_before*delta + sigma*s_before*sqrt(delta)*z[zIndex];
            zIndex++;
            here.push_back(s_new);
            s_before = s_new;
        }
        sPaths.push_back(here);
    }
    return sPaths;
}

//generates normal values
vector<double> normGenBM(int num){
    int a = pow(7,5);
    int m = pow(2,31) - 1;
    int b = 0;
    vector<double> uniform;
    long int xNow = 30000;
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

double fixedStrikeLookbackCall(double r, double s0, double k, double sigma, double t, int m, int n){
    vector<vector<double>> stock_paths = stockPaths(m, n, t, s0, sigma, r);
    vector<double> s_max;
    for(vector<double> path : stock_paths){
        s_max.push_back(vec_max(path));
    }
    double sum = 0;
    for (int i = 0; i < s_max.size(); i++) {
        double payoff = s_max[i] - k;
        if (payoff < 0) {
            payoff = 0;
        }
        sum += payoff;
    }
    double result = sum / m;
    return result;
}

double fixedStrikeLookbackPut(double r, double s0, double k, double sigma, double t, int m, int n){
    vector<vector<double>> stock_paths = stockPaths(m, n, t, s0, sigma, r);
    vector<double> s_min;
    for(vector<double> path : stock_paths){
        s_min.push_back(vec_min(path));
    }
    double sum = 0;
    for (int i = 0; i < s_min.size(); i++) {
        double payoff = k - s_min[i];
        if (payoff < 0) {
            payoff = 0;
        }
        sum += payoff;
    }
    double result = sum / m;
    return exp(-1.0*r*t)*result;
}

double vec_max(const vector<double>& v){
    double max = v[0];
    for(int i = 1; i < v.size(); i++){
        if(v[i] > max){
            max = v[i];
        }
    }
    return max;
}

double vec_min(const vector<double>& v){
    double min = v[0];
    for(int i = 1; i < v.size(); i++){
        if(v[i] < min){
            min = v[i];
        }
    }
    return min;
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

vector<double> expGen(double lambda, int num, int seed){
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
    for(double i : uniform){
        result.push_back((-1.0)*lambda*log(1-i));
    }
    return result;
}

vector<vector<double>> stockPathsJumps(int m, int n, double t, double s0, double sigma, double r, double lambda, double gamma){
//    int exp_index = 0;
//    vector<double> exp_nums = expGen(1.0/lambda, m*n*2, 77777);
//    vector<vector<double>> jump_times;
//    for(int i = 0; i < m; i++){
//        vector<double> jumps_here;
//        double exp_val = 0;
//        bool run = true;
//        while(run){
//            exp_val += exp_nums[exp_index];
//            exp_index++;
//            if(exp_val < t){
//                jumps_here.push_back(exp_val);
//            }
//            else{
//                run = false;
//            }
//        }
//        jump_times.push_back(jumps_here);
//    }
//    vector<double> z = normGenBM(m*n);
//    vector<vector<double>> sPaths;
//    double delta = t/n;
//    int zIndex = 0;
//    for(int i = 0; i < m; i++){
//        vector<double> here;
//        double s_before = s0;
//        here.push_back(s_before);
//        vector<double> jumptimes = jump_times[i];
//        if(jumptimes.size() == 0){
//            for(int j = 0; j < n; j++){
//                double s_new = s_before + r*s_before*delta + sigma*s_before*sqrt(delta)*z[zIndex];
//                zIndex++;
//                here.push_back(s_new);
//                s_before = s_new;
//            }
//            sPaths.push_back(here);
//        }
//        else{
//            double next_jump = jumptimes[0];
//            int curr_jump_index = 0;
//            int j = 0;
//            while(j < n){
//                double s_new = s_before + r*s_before*delta + sigma*s_before*sqrt(delta)*z[zIndex];
//                zIndex++;
//                here.push_back(s_new);
//                j++;
//                double time_here = (j+1.0)*delta;
//                if((time_here + delta) > next_jump && j != n){
//                    s_new += gamma*s_new;
//                    curr_jump_index++;
//                    here.push_back(s_new);
//                    j++;
//                    if(curr_jump_index > (jumptimes.size() - 1)){
//                        next_jump = t*2.0;
//                    }
//                    else{
//                        next_jump = jumptimes[curr_jump_index];
//                    }
//                }
//                s_before = s_new;
//            }
//            sPaths.push_back(here);
//        }
//    }
    double delta = t / n;
    vector<double> z = normGenBM(m*n);
    vector<double> N = poisson(lambda*delta, m*n, 12345);
    vector<vector<double>> M;
    double poi_index = 0;
    for(int i = 0; i < m; i++){
        vector<double> m_here;
        for(int j = 0; j < n; j++){
            if(N[poi_index] != 0){
                m_here.push_back(N[poi_index]);
            }
            else{
                m_here.push_back(0);
            }
            poi_index++;
        }
        M.push_back(m_here);
    }
    int zIndex = 0;
    vector<vector<double>> sPaths;
    for(int i = 0; i < m; i++){
        vector<double> here;
        double s_before = s0;
        here.push_back(s_before);
        for(int j = 0; j < n; j++){
            double s_new = s_before + r*s_before*delta + sigma*sqrt(delta)*z[zIndex]*s_before;
            s_new += s_new * M[i][j]*gamma;
            zIndex++;
            here.push_back(s_new);
            s_before = s_new;
        }
        sPaths.push_back(here);
    }
    
    return sPaths;
}

void Proj6_2function(double lambda1, double lambda2, double T, double& value, double& prob, double& Et){
    double v0 = 20000.0;
    double l0 = 22000.0;
    double mu = -0.1;
    double sigma = 0.2;
    double gamma = -0.4;
    double delta = 0.25;
    double r0 = 0.02 + delta * lambda2;
    double alpha = 0.7;
    double epsilon = 0.95;
    double r_month = r0 / 12.0;
    int n = T * 12;
    int m = 10000;
    double pmt = l0*r_month/(1.0-1.0/power(1.0+r_month,n));
    double a = pmt / r_month;
    double b = pmt / (r_month*power(1.0+r_month,n));
    double c = 1.0 + r_month;
    double beta = (epsilon - alpha) / T;
    vector<double> l_t{l0};
    for(int i = 0; i < n; i++){
        l_t.push_back(a-b*power(c, i+1));
    }

    vector<vector<double>> collateral_paths = stockPathsJumps(m, n, T, v0, sigma, mu, lambda1, gamma);
    
    vector<double> q_t;
    for(double t = 0; t <= T; t += T/n){
        q_t.push_back(alpha+beta*t);
    }
    vector<double> poi = poisson(lambda2*delta, m*n, 777777);
    
    vector<vector<double>> s_values;
    double poi_index = 0;
    for(int i = 0; i < m; i++){
        vector<double> temp;
        for(int j = 0; j < n; j++){
            temp.push_back(poi[poi_index]);
            poi_index++;
        }
        s_values.push_back(temp);
    }
    vector<double> payoffs;
    vector<double> expected_times;
    int default_times = 0;
    for(int i = 0; i < m; i++){
        int Q_index = 2*n;
        for(int j = 0; j <= n; j++){
            double v_here = collateral_paths[i][j];
            double l_here = l_t[j];
            double q_here = q_t[j];
            if(v_here <= l_here*q_here){
                Q_index = j;
                break;
            }
        }
        double Q = Q_index * (T / n);
        int S_index = 2*n;
        for(int j = 0; j < n; j++){
            if(s_values[i][j] > 0){
                S_index = j + 1;
                break;
            }
        }
        double S = S_index * (T/n);
        if(S_index < 0 && Q_index < 0){
            cout << S_index << Q_index << endl;
        }
        
        double payoff_here = 0;
        if(Q <= S && Q <= T){
            payoff_here = l_t[Q_index] - epsilon * collateral_paths[i][Q_index];
            if(payoff_here < 0){
                payoff_here = 0;
            }
            default_times++;
            expected_times.push_back(Q);
            payoffs.push_back(payoff_here*exp(-1.0*r_month*Q_index));
            
        }
        else if (S < Q && S <= T){
            payoff_here = l_t[S_index] - epsilon*collateral_paths[i][S_index];
            if(payoff_here < 0){
                payoff_here *= -1.0;
            }
            default_times++;
            expected_times.push_back(S);
            payoffs.push_back(payoff_here*exp(-1.0*r_month*S_index));
        }
        else{
            expected_times.push_back(T);
            payoffs.push_back(0);
        }
        
        
    }
    value = calcMean(payoffs);
    double times = default_times;
    prob = times / m;
    Et = calcMean(expected_times);
    
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

//returns mean of vector of doubles
double calcMean(const vector<double>& v){
    double size = v.size();
    double sum = 0.0;
    for(double i : v){
        sum += i;
    }
    return (sum / size);
}

vector<double> poisson(double lambda, int num, int seed){
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
    for(double u : uniform){
        double z = exp(-1.0*lambda);
        double k = 0;
        double x = z;
        bool run = true;
        while(run){
            if(u < x){
                result.push_back(k);
                run = false;
            }
            else{
                z *= lambda/(k+1.0);
                x += z;
                k++;
            }
        }
    }
    return result;
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
