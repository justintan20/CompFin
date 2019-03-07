//
//  justin.cpp
//  405HW9
//
//  Created by Justin Tan on 3/6/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#include "justin.h"

double priceMBSNumerix(unsigned int T_years, double loan_amount, double wac, double r0, double kappa, double r_bar, double sigma, int num_paths){
    unsigned int numMonths = T_years * 12;
    double monthlyMortgageRate = wac / 12.0;
    double delta = 1/12.0;
    double pmt = loan_amount*monthlyMortgageRate/(1.0-1.0/(power(1.0+monthlyMortgageRate, numMonths)));
    
    vector<vector<double>> interest_rates = ratesCIR(r0, kappa, r_bar, sigma, num_paths, numMonths, T_years);
    
    long double price_sum = 0;
    for(int j = 0; j < num_paths; j++){
        vector<double> ir_here = interest_rates[j];
        
        vector<double> PV_t{loan_amount};
        vector<double> c_t{0};
        int currMonth = 0;
        
        
        for(int i = 0; i < numMonths; i++){
            currMonth++;
            if(currMonth == 13){
                currMonth = 1;
            }
            if(PV_t[i] <= 0){
                break;
            }
            double IP_t = PV_t[i]*monthlyMortgageRate;
            double SP_t = pmt - IP_t;
            double r_tenyear = 0;
            double bondPrice = ZCBClosedFormCIR(ir_here[i], r_bar, sigma, kappa, 1, 0, 10);
            r_tenyear = -1.0/10.0*log(bondPrice);
            double CPR_t = cprNumerixCIR(r_tenyear, wac, loan_amount, PV_t[i], i+1, currMonth);
            double PP_t = (PV_t[i] - SP_t)*(1.0-power(1.0-CPR_t, 1.0/12.0));
            double remaining = PV_t[i] - SP_t - PP_t;
            if(remaining < 0){
                double extra = 0 - remaining;
                PP_t -= extra;
                remaining = 0;
            }
            double cash = pmt + PP_t;
            PV_t.push_back(remaining);
            c_t.push_back(cash);
        }
        
        long double price = 0;
        int duration = c_t.size() - 1;
        for(int i = 1; i <= duration; i++){
            double cash_here = c_t[i];
            double r_sum = 0;
            for(int index = 1; index <= i; index++){
                r_sum += ir_here[index];
            }
            double d_t = exp(-delta*r_sum);
            double value = cash_here*d_t;
            price += value;
        }
        price_sum += price;
    }
    double result = price_sum / num_paths;
    return result;
}



double cprNumerixCIR(double r, double mort_rate, double pv0, double pv_tMinus1, unsigned int curr_period, unsigned int month){
    double RI_t = 0.28 + 0.14*atan(-8.57+430.0*(mort_rate - r));
    double BU_t = 0.3 + 0.7*pv_tMinus1/pv0;
    double SG_t = curr_period / 30.0;
    if(SG_t > 1){
        SG_t = 1;
    }
    double SY_t = 0;
    switch (month) {
        case 1:
            SY_t = 0.94;
            break;
        
        case 2:
            SY_t = 0.76;
            break;
            
        case 3:
            SY_t = 0.74;
            break;
        
        case 4:
            SY_t = 0.95;
            break;
            
        case 5:
            SY_t = 0.98;
            break;
            
        case 6:
            SY_t = 0.92;
            break;
            
        case 7:
            SY_t = 0.98;
            break;
            
        case 8:
            SY_t = 1.1;
            break;
            
        case 9:
            SY_t = 1.18;
            break;
            
        case 10:
            SY_t = 1.22;
            break;
            
        case 11:
            SY_t = 1.23;
            break;
            
        case 12:
            SY_t = 0.98;
            break;
    }
    double cpr = RI_t*BU_t*SG_t*SY_t;
    return cpr;
}


vector<vector<double>> ratesCIR(double r0, double kappa, double r_bar, double sigma, int num_paths, int num_steps, int T_years){
    vector<double> z = normGenBM(num_paths*num_steps, 12345);
    int zIndex = 0;
    vector<vector<double>> allPaths;
    double delta = T_years / ((double)num_steps);
    for(int i = 0; i < num_paths; i++){
        vector<double> thisPath{r0};
        for(int j = 0; j < num_steps; j++){
            double r_before = thisPath[j];
            double r_here = r_before+kappa*(r_bar-r_before)*delta;
            r_here += sigma*sqrt(r_before)*sqrt(delta)*z[zIndex];
            zIndex++;
            if(r_here < 0){
                r_here *= -1;
            }
            thisPath.push_back(r_here);
        }
        allPaths.push_back(thisPath);
    }
    return allPaths;
}

double ZCBClosedFormCIR(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T){
    double h1 = sqrt(kappa*kappa + 2*sigma*sigma);
    double h2 = (kappa + h1)/2;
    double h3 = 2*kappa*r_bar/sigma/sigma;
    double a = power(h1*exp(h2*(T-t))/(h2*(exp(h1*(T-t))-1)+h1), h3);
    double b = (exp(h1*(T-t))-1)/(h2*(exp(h1*(T-t))-1)+h1);
    double price = a*exp(-b*r0)*faceValue;
    return price;
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
