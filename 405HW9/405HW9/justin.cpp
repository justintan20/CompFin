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
    
    vector<vector<double>> interest_rates = ratesCIR(r0, kappa, r_bar, sigma, num_paths, numMonths, T_years);
    
    long double price_sum = 0;
    for(int j = 0; j < num_paths; j++){
        vector<double> ir_here = interest_rates[j];
        
        vector<double> PV_t{loan_amount};
        int currMonth = 0;
        double price_here = 0;
        double r_sum = 0;
        for(int i = 0; i < numMonths; i++){
            currMonth++;
            if(currMonth == 13){
                currMonth = 1;
            }
            if(PV_t[i] <= 0){
                break;
            }
            r_sum += ir_here[i+1];
            double IP_t = PV_t[i]*monthlyMortgageRate;
            double MP_t = PV_t[i]*monthlyMortgageRate/(1.0-1.0/power(1.0+monthlyMortgageRate, numMonths - i));
            double SP_t = MP_t - IP_t;
            double r_tenyear = 0;
            double bondPrice = ZCBClosedFormCIR(ir_here[i], r_bar, sigma, kappa, 1, 0, 10);
            r_tenyear = -1.0/10.0*log(bondPrice);
            double CPR_t = cprNumerixCIR(r_tenyear, wac, loan_amount, PV_t[i], i+1, currMonth);
            double PP_t = (PV_t[i] - SP_t)*(1.0-power(1.0-CPR_t, 1.0/12.0));
            double remaining = PV_t[i] - SP_t - PP_t;
            double c_t = MP_t + PP_t;
            PV_t.push_back(remaining);
            price_here += c_t*exp(-delta*r_sum);
            
        }
        
        price_sum += price_here;
    }
    double result = price_sum / num_paths;
    return result;
}

double oas(unsigned int T_years, double loan_amount, double wac, double r0, double kappa, double r_bar, double sigma, int num_paths, double marketPrice){
    unsigned int numMonths = T_years * 12;
    double monthlyMortgageRate = wac / 12.0;
    double delta = 1/12.0;
    
    vector<vector<double>> interest_rates = ratesCIR(r0, kappa, r_bar, sigma, num_paths, numMonths, T_years);
    vector<vector<double>> cashFlowMatrix;
    for(int j = 0; j < num_paths; j++){
        vector<double> ir_here = interest_rates[j];
        vector<double> cashflow{0};
        
        vector<double> PV_t{loan_amount};
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
            double MP_t = PV_t[i]*monthlyMortgageRate/(1.0-1.0/power(1.0+monthlyMortgageRate, numMonths - i));
            double SP_t = MP_t - IP_t;
            double r_tenyear = 0;
            double bondPrice = ZCBClosedFormCIR(ir_here[i], r_bar, sigma, kappa, 1, 0, 10);
            r_tenyear = -1.0/10.0*log(bondPrice);
            double CPR_t = cprNumerixCIR(r_tenyear, wac, loan_amount, PV_t[i], i+1, currMonth);
            double PP_t = (PV_t[i] - SP_t)*(1.0-power(1.0-CPR_t, 1.0/12.0));
            double remaining = PV_t[i] - SP_t - PP_t;
            double c_t = MP_t + PP_t;
            PV_t.push_back(remaining);
            cashflow.push_back(c_t);
        }
        cashFlowMatrix.push_back(cashflow);
    }
    //bisection with 0.01 precision on price
    double x_up = 0.1;
    double x_down = -0.1;
    double x = 0;
    bool solved = false;
    while(!solved){
        double payoffSum = 0;
        for(int path = 0; path < cashFlowMatrix.size(); path++){
            double priceHere = 0;
            vector<double> cfHere = cashFlowMatrix[path];
            for(int i = 1; i < cfHere.size(); i++){
                double r_sum = 0;
                for(int j = 1; j < i + 1; j++){
                    r_sum += (interest_rates[path][j]+x);
                }
                double discounted_cf = cfHere[i]*exp(-delta*r_sum);
                priceHere += discounted_cf;
            }
            payoffSum += priceHere;
        }
        double computed_price = payoffSum / num_paths;
        if(computed_price > marketPrice){
            x_down = x;
            x = (x_up - x_down)/2.0 + x_down;
        }
        else if(computed_price < marketPrice){
            x_up = x;
            x = (x_up - x_down)/2.0 + x_down;
        }
        double diff = computed_price - marketPrice;
        if((diff < 0.01 && diff > 0)||(diff > -0.01 && diff < 0)){
            solved = true;
        }
    }
    
    return x;
}

double priceMBSNumerixWithOAS(unsigned int T_years, double loan_amount, double wac, double r0, double kappa, double r_bar, double sigma, int num_paths, double oasSpread){
    unsigned int numMonths = T_years * 12;
    double monthlyMortgageRate = wac / 12.0;
    double delta = 1/12.0;
    
    vector<vector<double>> interest_rates = ratesCIR(r0, kappa, r_bar, sigma, num_paths, numMonths, T_years);
    vector<vector<double>> cashFlowMatrix;
    for(int j = 0; j < num_paths; j++){
        vector<double> ir_here = interest_rates[j];
        vector<double> cashflow{0};
        
        vector<double> PV_t{loan_amount};
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
            double MP_t = PV_t[i]*monthlyMortgageRate/(1.0-1.0/power(1.0+monthlyMortgageRate, numMonths - i));
            double SP_t = MP_t - IP_t;
            double r_tenyear = 0;
            double bondPrice = ZCBClosedFormCIR(ir_here[i], r_bar, sigma, kappa, 1, 0, 10);
            r_tenyear = -1.0/10.0*log(bondPrice);
            double CPR_t = cprNumerixCIR(r_tenyear, wac, loan_amount, PV_t[i], i+1, currMonth);
            double PP_t = (PV_t[i] - SP_t)*(1.0-power(1.0-CPR_t, 1.0/12.0));
            double remaining = PV_t[i] - SP_t - PP_t;
            double c_t = MP_t + PP_t;
            PV_t.push_back(remaining);
            cashflow.push_back(c_t);
        }
        cashFlowMatrix.push_back(cashflow);
    }

    double payoffSum = 0;
    for(int path = 0; path < cashFlowMatrix.size(); path++){
        double priceHere = 0;
        vector<double> cfHere = cashFlowMatrix[path];
        for(int i = 1; i < cfHere.size(); i++){
            double r_sum = 0;
            for(int j = 1; j < i + 1; j++){
                r_sum += (interest_rates[path][j]+oasSpread);
            }
            double discounted_cf = cfHere[i]*exp(-delta*r_sum);
            priceHere += discounted_cf;
        }
        payoffSum += priceHere;
    }
    double computed_price = payoffSum / num_paths;
    
    return computed_price;
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
    vector<double> z = normGenBM(num_paths*num_steps, 111);
    int zIndex = 0;
    vector<vector<double>> allPaths;
    double delta = T_years / ((double)num_steps);
    for(int i = 0; i < num_paths; i++){
        vector<double> thisPath{r0};
        for(int j = 0; j < num_steps; j++){
            double r_before = thisPath[j];
            double r_here = r_before+kappa*(r_bar-r_before)*delta;
            if(r_before < 0){
                r_here += sigma*sqrt(-r_before)*sqrt(delta)*z[zIndex];
            }
            else{
                r_here += sigma*sqrt(r_before)*sqrt(delta)*z[zIndex];
            }
            zIndex++;
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
