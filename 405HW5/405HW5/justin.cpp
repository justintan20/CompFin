//
//  justin.cpp
//  405HW5
//
//  Created by Justin Tan on 2/14/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#include "justin.h"

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

vector<double> solveLinearSystem(vector<vector<double>>& A, vector<double>& b){
    vector<double> x;
    unsigned int num = A[0].size();
    for(int m = 0; m < num; m++){
        vector<double> row_1 = A[m];
        double a_1_m = row_1[m];
        double b_1 = b[m];
        for(int i = m+1; i < A.size(); i++){
            vector<double> oldRow = A[i];
            double a_i_m = oldRow[m];
            double old_b = b[i];
            vector<double> newRow = vectorSubtract(oldRow, vectorMultiplication(row_1, a_i_m / a_1_m));
            A[i] = newRow;
            double new_b = old_b - (a_i_m / a_1_m)*b_1;
            b[i] = new_b;
        }
    }
    vector<double> results;
    results.push_back(b[A.size()-1] / A[A.size()-1][A.size()-1]);
    for(int j = A.size()-2; j >=0; j -= 1){
        vector<double> thisRow = A[j];
        double this_a = thisRow[j];
        double bNum = b[j];
        for(int constantIndex = j + 1; constantIndex < A.size(); constantIndex++){
            bNum -= (thisRow[constantIndex] * results[A.size()-constantIndex-1]);
        }
        double this_x = bNum / this_a;
        results.push_back(this_x);
    }
    for(int i = 0; i < results.size(); i++){
        x.push_back(results[results.size()-i-1]);
    }
    return x;
}

double LSMC(double s0, double k, double r, double sigma, double t, int m, int n, int num_func, const string& type){
    double delta = t / n;
    //stock price and index matrices
    vector<vector<double>>sMatrix = stockPaths(m, n, t, s0, sigma, r);
    vector<vector<double>> indexMatrix;
    //initialize index matrix
    for(int i = 0; i < m; i++){
        vector<double> row;
        for(int j = 0; j <= n; j++){
            row.push_back(0);
        }
        indexMatrix.push_back(row);
    }
    //start at t=n
    int currCol = n;
    for(int row = 0; row < m; row++){
        double payoff = k - sMatrix[row][currCol];
        if(payoff > 0){
            indexMatrix[row][currCol] = 1;
        }
    }
    
    //t = n-1
    currCol = n - 1;
    while(currCol >= 0){
        vector<double> x;
        vector<double> y;
        for(int i = 0; i < m; i++){
            double val = 0;
            for(int time = currCol + 1; time <= n; time++){
                if(indexMatrix[i][time] == 1){
                    val = exp(-1.0*r*(time-currCol)*delta)*(k - sMatrix[i][time]);
                }
            }
            //only use non-zero paths
            if(val != 0){
                y.push_back(val);
                x.push_back(sMatrix[i][currCol]);
            }
        }
        
        vector<vector<double>> f;
        vector<double> L1;
        for(int i = 0; i < x.size(); i++){
            if(type == "laguerre"){
                L1.push_back(exp(-1.0*x[i]/2.0));
            }
            else if(type == "hermite" || type == "monomial"){
                L1.push_back(1);
            }
        }
        f.push_back(L1);
        vector<double> L2;
        for(int i = 0; i < x.size(); i++){
            if(type == "laguerre"){
                L2.push_back(exp(-1.0*x[i]/2.0)*(1.0-x[i]));
            }
            else if(type == "hermite"){
                L2.push_back(2.0*x[i]);
            }
            else if(type == "monomial"){
                L2.push_back(x[i]);
            }
        }
        f.push_back(L2);
        if(num_func >= 3){
            vector<double> L3;
            for(int i = 0; i < x.size(); i++){
                if(type == "laguerre"){
                    L3.push_back(exp(-1.0*x[i]/2.0)*(1.0-2.0*x[i]+x[i]*x[i]/2.0));
                }
                else if(type == "hermite"){
                    L3.push_back(4.0*x[i]*x[i]-2.0);
                }
                else if(type == "monomial"){
                    L3.push_back(x[i]*x[i]);
                }
            }
            f.push_back(L3);
        }
        if(num_func >= 4){
            vector<double> L4;
            for(int i = 0; i < x.size(); i++){
                if(type == "laguerre"){
                    L4.push_back(exp(-1.0*x[i]/2.0)*(1.0-3.0*x[i]+3*x[i]*x[i]/2.0 - x[i]*x[i]*x[i]/6.0));
                }
                else if(type == "hermite"){
                    L4.push_back(8.0*x[i]*x[i]*x[i]-12.0*x[i]);
                }
                else if(type == "monomial"){
                    L4.push_back(x[i]*x[i]*x[i]);
                }
            }
            f.push_back(L4);
        }
        
        vector<double> b;
        for(vector<double> f_k : f){
            b.push_back(dotProd(y, f_k));
        }
        
        vector<vector<double>> A;
        unsigned int A_dim = f.size();
        for(int i = 0; i < A_dim; i++){
            vector<double> row;
            for(int j = 0; j < A_dim; j++){
                row.push_back(dotProd(f[j], f[i]));
            }
            A.push_back(row);
        }
        
        vector<double> coef = solveLinearSystemLU(A, b);
        for(int i = 0; i < m; i++){
            double s_here = sMatrix[i][currCol];
            double exercise_value = k - s_here;
            if(exercise_value < 0){
                exercise_value = 0;
            }
            double continuation_value = 0;
            if(type == "laguerre"){
                continuation_value += coef[0]*exp(-1.0*s_here/2.0);
                continuation_value += coef[1]*exp(-1.0*s_here/2.0)*(1.0-s_here);
                if(num_func >= 3){
                    continuation_value += coef[2]*exp(-1.0*s_here/2.0)*(1.0-2*s_here+s_here*s_here/2.0);
                }
                if(num_func >= 4){
                    continuation_value += coef[3]*exp(-1.0*s_here/2.0)*(1.0-3.0*s_here+3.0*s_here*s_here/2.0+s_here*s_here*s_here/6.0);
                }
            }
            else if(type == "hermite"){
                continuation_value += coef[0];
                continuation_value += coef[1]*2.0*s_here;
                if(num_func >= 3){
                    continuation_value += coef[2]*(4.0*s_here*s_here-2.0);
                }
                if(num_func >= 4){
                    continuation_value += coef[3]*(8.0*s_here*s_here*s_here-12.0*s_here);
                }
            }
            else if(type == "monomial"){
                continuation_value += coef[0];
                continuation_value += coef[1]*s_here;
                if(num_func >= 3){
                    continuation_value += coef[2]*s_here*s_here;
                }
                if(num_func >= 4){
                    continuation_value += coef[3]*s_here*s_here*s_here;
                }
            }
            if(continuation_value < exercise_value && continuation_value > 0){
                indexMatrix[i][currCol] = 1;
                for(int time = currCol+1; time <= n; time++){
                    if(indexMatrix[i][time] == 1){
                        indexMatrix[i][time] = 0;
                    }
                }
            }
        }
        currCol -= 1;
    }
    double sum = 0;
    for(int i = 0; i < m; i++){
        for(int j = 0; j <= n; j++){
            if(indexMatrix[i][j] == 1){
                sum += (exp(-1.0*r*j*delta)*(k-sMatrix[i][j]));
            }
        }
    }
    double result = sum / m;
    return result;
}

vector<vector<double>> stockPaths(int m, int n, double t, double s0, double sigma, double r){
    vector<double> z = normGenBM(m*n/2);
    int size = m*n/2;
    for(int j = 0; j < size; j++){
        z.push_back(z[j]*(-1.0));
    }
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

//calculates dot product of vectors
double dotProd(const vector<double>& a, const vector<double>& b){
    double sum = 0;
    for(int i = 0; i < a.size(); i++){
        sum += a[i]*b[i];
    }
    return sum;
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

double forwardStartEuroPut(double s0, double sigma, double r, double t, double T, int n, int m){
    vector<vector<double>> s_paths = stockPaths(m, n, T, s0, sigma, r);
    int t_period = t/(T/n);
    double payoff_sum = 0;
    for(int i = 0; i < m; i++){
        double payoff = s_paths[i][t_period] - s_paths[i][n];
        if(payoff < 0){
            payoff = 0;
        }
        payoff_sum += payoff;
    }
    return exp(-r*T)*payoff_sum / m;
    
}

double forwardStartAmericanPut(double s0, double sigma, double r, double t, double t_ex, int m, int n, int num_func, const string& type){
    double delta = t / n;
    //stock price and index matrices
    vector<vector<double>>sMatrix = stockPaths(m, n, t, s0, sigma, r);
    int t_period = t_ex/delta;
    vector<double> exercise_prices;
    for(int i = 0; i < m; i++){
        exercise_prices.push_back(sMatrix[i][t_period]);
    }
    vector<vector<double>> indexMatrix;
    //initialize index matrix
    for(int i = 0; i < m; i++){
        vector<double> row;
        for(int j = 0; j <= n; j++){
            row.push_back(0);
        }
        indexMatrix.push_back(row);
    }
    //start at t=n
    int currCol = n;
    for(int row = 0; row < m; row++){
        double payoff = exercise_prices[row] - sMatrix[row][currCol];
        if(payoff > 0){
            indexMatrix[row][currCol] = 1;
        }
    }
    
    //t = n-1
    currCol = n - 1;
    while(currCol >= t_period){
        vector<double> x;
        vector<double> y;
        for(int i = 0; i < m; i++){
            double val = 0;
            for(int time = currCol + 1; time <= n; time++){
                if(indexMatrix[i][time] == 1){
                    val = exp(-1.0*r*(time-currCol)*delta)*(exercise_prices[i] - sMatrix[i][time]);
                }
            }
            //only use non-zero paths
            if(val != 0){
                y.push_back(val);
                x.push_back(sMatrix[i][currCol]);
            }
        }
        
        vector<vector<double>> f;
        vector<double> L1;
        for(int i = 0; i < x.size(); i++){
            if(type == "laguerre"){
                L1.push_back(exp(-1.0*x[i]/2.0));
            }
            else if(type == "hermite" || type == "monomial"){
                L1.push_back(1);
            }
        }
        f.push_back(L1);
        vector<double> L2;
        for(int i = 0; i < x.size(); i++){
            if(type == "laguerre"){
                L2.push_back(exp(-1.0*x[i]/2.0)*(1.0-x[i]));
            }
            else if(type == "hermite"){
                L2.push_back(2.0*x[i]);
            }
            else if(type == "monomial"){
                L2.push_back(x[i]);
            }
        }
        f.push_back(L2);
        if(num_func >= 3){
            vector<double> L3;
            for(int i = 0; i < x.size(); i++){
                if(type == "laguerre"){
                    L3.push_back(exp(-1.0*x[i]/2.0)*(1.0-2.0*x[i]+x[i]*x[i]/2.0));
                }
                else if(type == "hermite"){
                    L3.push_back(4.0*x[i]*x[i]-2.0);
                }
                else if(type == "monomial"){
                    L3.push_back(x[i]*x[i]);
                }
            }
            f.push_back(L3);
        }
        if(num_func >= 4){
            vector<double> L4;
            for(int i = 0; i < x.size(); i++){
                if(type == "laguerre"){
                    L4.push_back(exp(-1.0*x[i]/2.0)*(1.0-3.0*x[i]+3*x[i]*x[i]/2.0 - x[i]*x[i]*x[i]/6.0));
                }
                else if(type == "hermite"){
                    L4.push_back(8.0*x[i]*x[i]*x[i]-12.0*x[i]);
                }
                else if(type == "monomial"){
                    L4.push_back(x[i]*x[i]*x[i]);
                }
            }
            f.push_back(L4);
        }
        
        vector<double> b;
        for(vector<double> f_k : f){
            b.push_back(dotProd(y, f_k));
        }
        
        vector<vector<double>> A;
        unsigned int A_dim = f.size();
        for(int i = 0; i < A_dim; i++){
            vector<double> row;
            for(int j = 0; j < A_dim; j++){
                row.push_back(dotProd(f[j], f[i]));
            }
            A.push_back(row);
        }
        
        vector<double> coef = solveLinearSystemLU(A, b);
        for(int i = 0; i < m; i++){
            double s_here = sMatrix[i][currCol];
            double exercise_value = exercise_prices[i] - s_here;
            if(exercise_value < 0){
                exercise_value = 0;
            }
            double continuation_value = 0;
            if(type == "laguerre"){
                continuation_value += coef[0]*exp(-1.0*s_here/2.0);
                continuation_value += coef[1]*exp(-1.0*s_here/2.0)*(1.0-s_here);
                if(num_func >= 3){
                    continuation_value += coef[2]*exp(-1.0*s_here/2.0)*(1.0-2*s_here+s_here*s_here/2.0);
                }
                if(num_func >= 4){
                    continuation_value += coef[3]*exp(-1.0*s_here/2.0)*(1.0-3.0*s_here+3.0*s_here*s_here/2.0+s_here*s_here*s_here/6.0);
                }
            }
            else if(type == "hermite"){
                continuation_value += coef[0];
                continuation_value += coef[1]*2.0*s_here;
                if(num_func >= 3){
                    continuation_value += coef[2]*(4.0*s_here*s_here-2.0);
                }
                if(num_func >= 4){
                    continuation_value += coef[3]*(8.0*s_here*s_here*s_here-12.0*s_here);
                }
            }
            else if(type == "monomial"){
                continuation_value += coef[0];
                continuation_value += coef[1]*s_here;
                if(num_func >= 3){
                    continuation_value += coef[2]*s_here*s_here;
                }
                if(num_func >= 4){
                    continuation_value += coef[3]*s_here*s_here*s_here;
                }
            }
            if(continuation_value < exercise_value && continuation_value > 0){
                indexMatrix[i][currCol] = 1;
                for(int time = currCol+1; time <= n; time++){
                    if(indexMatrix[i][time] == 1){
                        indexMatrix[i][time] = 0;
                    }
                }
            }
        }
        currCol -= 1;
    }
    double sum = 0;
    for(int i = 0; i < m; i++){
        for(int j = 0; j <= n; j++){
            if(indexMatrix[i][j] == 1){
                sum += (exp(-1.0*r*j*delta)*(exercise_prices[i]-sMatrix[i][j]));
            }
        }
    }
    double result = sum / m;
    return result;
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
