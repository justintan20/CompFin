//
//  justin.cpp
//  405HW8
//
//  Created by Justin Tan on 3/3/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#include "justin.h"

double ZCBVasicek(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T, int steps, int paths){
    vector<double> z = normGenBM(steps*paths,54321);
    int zIndex = 0;
    double delta = T / steps;
    double inner_sum = 0;
    for(int i = 0; i < paths; i++){
        double r_before = r0;
        double sum = 0;
        for(int j = t*steps; j < steps; j++){
            double r_after = r_before + kappa*(r_bar - r_before)*delta + sigma*sqrt(delta)*z[zIndex];
            zIndex++;
            sum += r_after;
            r_before = r_after;
        }
        double value = exp(-delta*sum)*faceValue;
        inner_sum += value;
    }
    double result = inner_sum / paths;
    return result;
}

double ZCBClosedFormVasicek(double r_t, double r_bar, double sigma, double kappa, double faceValue, double t, double T){
    double b = 1/kappa*(1-exp(-kappa*(T-t)));
    double a = exp((r_bar - sigma*sigma/2/kappa/kappa)*(b - (T - t)) - sigma*sigma/4/kappa*b*b);
    double price = a*exp(-b*r_t)*faceValue;
    return price;
}

double couponBondVasicek(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T, const vector<double>& c_i, const vector<double>& T_i, int steps, int paths){
    double sum = 0;
    for(int i = 0; i < c_i.size(); i++){
        sum += ZCBVasicek(r0, r_bar, sigma, kappa, c_i[i], t, T_i[i], steps, paths);
    }
    double result = sum;
    return result;
}

double callZCBVasicek(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T, double S, double k, int steps, int paths){
    vector<double> z = normGenBM(steps*paths,77777);
    int zIndex = 0;
    double delta = S / steps;
    double inner_sum = 0;
    for(int i = 0; i < paths; i++){
        double r_before = r0;
        double sum = 0;
        double price_bond = 0;
        for(int j = t*steps; j < steps; j++){
            double curr_time = j * delta;
            if(curr_time == T){
                price_bond = ZCBClosedFormVasicek(r_before, r_bar, sigma, kappa, faceValue, T, S);
                break;
            }
        
            double r_after = r_before + kappa*(r_bar - r_before)*delta + sigma*sqrt(delta)*z[zIndex];
            zIndex++;
            sum += r_after;
            r_before = r_after;
        }
        double call_val = price_bond - k;
        if(call_val < 0){
            call_val = 0;
        }
        double value = exp(-delta*sum)*call_val;
        inner_sum += value;
    }
    double result = inner_sum / paths;
    return result;
}

double callZCBClosedFormVasicek(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T, double S, double k){
    double p_t_s = ZCBClosedFormVasicek(r0, r_bar, sigma, kappa, 1, t, S);
    double p_t_t = ZCBClosedFormVasicek(r0, r_bar, sigma, kappa, 1, t, T);
    double sigma_p = sqrt((1.0-exp(-2.0*kappa*(T-t)))/2/kappa)*(1-exp(-kappa*(S-T)))/kappa*sigma;
    double d1 = log(faceValue*p_t_s / k / p_t_t)/sigma_p + sigma_p/2.0;
    double d2 = d1 - sigma_p;
    double result = faceValue*p_t_s*cdfNorm(d1)-k*p_t_t*cdfNorm(d2);
    return result;
}

double callCouponBondVasicek(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T, double S, double k, const vector<double>& c_i, const vector<double>& T_i, int steps, int paths){
    vector<double> z = normGenBM(steps*paths,19333);
    int zIndex = 0;
    double delta = S / steps;
    double inner_sum = 0;
    for(int i = 0; i < paths; i++){
        double r_before = r0;
        double sum = 0;
        double price_bond = 0;
        for(int j = t*steps; j < steps; j++){
            double curr_time = j * delta;
            if(curr_time == T){
                price_bond = couponBondVasicek(r_before, r_bar, sigma, kappa, faceValue, 0, S - T, c_i, T_i, steps, paths);
                break;
            }
            double r_after = r_before + kappa*(r_bar - r_before)*delta + sigma*sqrt(delta)*z[zIndex];
            zIndex++;
            sum += r_after;
            r_before = r_after;
        }
        double call_val = price_bond - k;
        if(call_val < 0){
            call_val = 0;
        }
        double value = exp(-delta*sum)*call_val;
        inner_sum += value;
    }
    double result = inner_sum / paths;
    return result;
}

double ZCBCIR(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T, int steps, int paths){
    vector<double> z = normGenBM(steps*paths,52675);
    int zIndex = 0;
    double delta = T / steps;
    double inner_sum = 0;
    for(int i = 0; i < paths; i++){
        double r_before = r0;
        double sum = 0;
        for(int j = t*steps; j < steps; j++){
            double r_after = 0;
            if(r_before < 0){
                r_after = r_before + kappa*(r_bar - r_before)*delta + sigma*sqrt(delta)*sqrt(-r_before)*z[zIndex];
            }
            else{
                r_after = r_before + kappa*(r_bar - r_before)*delta + sigma*sqrt(delta)*sqrt(r_before)*z[zIndex];
            }
            zIndex++;
            sum += r_after;
            r_before = r_after;
        }
        double value = exp(-delta*sum)*faceValue;
        inner_sum += value;
    }
    double result = inner_sum / paths;
    return result;
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

//uses explicit bond price
double callZCBCIR(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T, double S, double k, int steps, int paths){
    vector<double> z = normGenBM(steps*paths,79584);
    int zIndex = 0;
    double delta = S / steps;
    double inner_sum = 0;
    for(int i = 0; i < paths; i++){
        double r_before = r0;
        double sum = 0;
        double price_bond = 0;
        for(int j = t*steps; j < steps; j++){
            double curr_time = j * delta;
            if(curr_time == T){
                price_bond = ZCBClosedFormCIR(r_before, r_bar, sigma, kappa, faceValue, T, S);
                break;
            }
            double r_after = 0;
            if(r_before < 0){
                r_after = r_before + kappa*(r_bar - r_before)*delta + sigma*sqrt(delta)*sqrt(-r_before)*z[zIndex];
            }
            else{
                r_after = r_before + kappa*(r_bar - r_before)*delta + sigma*sqrt(delta)*sqrt(r_before)*z[zIndex];
            }
            zIndex++;
            sum += r_after;
            r_before = r_after;
        }
        double call_val = price_bond - k;
        if(call_val < 0){
            call_val = 0;
        }
        double value = exp(-delta*sum)*call_val;
        inner_sum += value;
    }
    double result = inner_sum / paths;
    return result;
}

//uses monte carlo bond price
double callZCBCIR2(double r0, double r_bar, double sigma, double kappa, double faceValue, double t, double T, double S, double k, int steps, int paths){
    vector<double> z = normGenBM(steps*paths,88223);
    int zIndex = 0;
    double delta = S / steps;
    double inner_sum = 0;
    for(int i = 0; i < paths; i++){
        double r_before = r0;
        double sum = 0;
        double price_bond = 0;
        for(int j = t*steps; j < steps; j++){
            double curr_time = j * delta;
            if(curr_time == T){
                price_bond = ZCBCIR(r_before, r_bar, sigma, kappa, faceValue, T, S,steps,paths);
                break;
            }
            double r_after = 0;
            if(r_before < 0){
                r_after = r_before + kappa*(r_bar - r_before)*delta + sigma*sqrt(delta)*sqrt(-r_before)*z[zIndex];
            }
            else{
                r_after = r_before + kappa*(r_bar - r_before)*delta + sigma*sqrt(delta)*sqrt(r_before)*z[zIndex];
            }
            zIndex++;
            sum += r_after;
            r_before = r_after;
        }
        double call_val = price_bond - k;
        if(call_val < 0){
            call_val = 0;
        }
        double value = exp(-delta*sum)*call_val;
        inner_sum += value;
    }
    double result = inner_sum / paths;
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

double ZCBG2(double r0, double sigma, double faceValue, double t, double T, double a, double b, double x0, double y0, double phi0, double eta, double rho, int steps, int paths){
    vector<double> z = normGenBM(steps*paths*2,12345);
    double delta = T / steps;
    int zIndex1 = 0;
    int zIndex2 = steps*paths;
    double tot_sum = 0;
    for(int i = 0; i < paths; i++){
        double r_before = r0;
        double x_before = x0;
        double y_before = y0;
        double sum = 0;
        for(int j = 0; j < steps; j++){
            double x_after = x_before - a*x_before*delta + sigma*sqrt(delta)*z[zIndex1];
            double y_after = y_before - b*y_before*delta + eta*sqrt(delta)*(rho*z[zIndex1] + sqrt(1-rho*rho)*z[zIndex2]);
            zIndex1++;
            zIndex2++;
            double r_after = x_after + y_after + phi0;
            
            sum += r_after;
            r_before = r_after;
            x_before = x_after;
            y_before = y_after;
            

        }
        double value = exp(-delta*sum)*faceValue;
        tot_sum += value;
    }
    double price = tot_sum / paths;
    return price;
}

double putZCBG2(double r0, double sigma, double faceValue, double t, double T, double S, double k, double a, double b, double x0, double y0, double phi0, double eta, double rho, int steps, int paths){
    vector<double> z = normGenBM(steps*paths*2,30000);
    double delta = S / steps;
    int zIndex1 = 0;
    int zIndex2 = steps*paths;
    double tot_sum = 0;
    for(int i = 0; i < paths; i++){
        double r_before = r0;
        double x_before = x0;
        double y_before = y0;
        double sum = 0;
        double price_bond = 0;
        for(int j = 0; j < steps; j++){
            double curr_time = (j+1) * delta;
            double x_after = x_before - a*x_before*delta + sigma*sqrt(delta)*z[zIndex1];
            double y_after = y_before - b*y_before*delta + eta*sqrt(delta)*(rho*z[zIndex1] + sqrt(1-rho*rho)*z[zIndex2]);
            zIndex1++;
            zIndex2++;
            double r_after = x_after + y_after + phi0;
            sum += r_after;
            r_before = r_after;
            x_before = x_after;
            y_before = y_after;
            if(curr_time == T){
                price_bond = ZCBG2(r_before, sigma, faceValue, T, S, a, b, x_before, y_before, phi0, eta, rho,steps,paths);
                break;
            }
        }
        double put = k - price_bond;
        if(put < 0){
            put = 0;
        }
        double value = exp(-delta*sum)*put;
        tot_sum += value;
    }
    double price = tot_sum / paths;
    return price;
}

double ZCBClosedFormG2(double r0, double sigma, double faceValue, double t, double T, double a, double b, double x0, double y0, double phi0, double eta, double rho){
    double v_left = sigma*sigma/a/a*(T-t + 2.0/a*exp(-a*(T-t)) - 1/2.0/a*exp(-2.0*a*(T-t)) - 3.0/2.0/a);
    double v_mid = eta*eta/b/b*(T-t + 2.0/b*exp(-b*(T-t)) - 1/2.0/b*exp(-2.0*b*(T-t)) - 3.0/2.0/b);
    double v_right = 2.0*rho*sigma*eta/a/b*(T-t + (exp(-a*(T-t))-1.0)/a + (exp(-b*(T-t))-1.0)/b - (exp(-(a+b)*(T-t))-1.0)/(a+b));
    double v = v_left+v_mid+v_right;
    double e1 = phi0*(T-t);
    double e2 = (1.0-exp(-a*(T-t)))/a*x0;
    double e3 = (1.0-exp(-b*(T-t)))/b*y0;
    double e4 = 0.5*v;
    double price = exp(-e1-e2-e3+e4)*faceValue;
    return price;
}

double putZCBClosedFormG2(double r0, double sigma, double faceValue, double t, double T, double S, double k, double a, double b, double x0, double y0, double phi0, double eta, double rho){
    double p_t_s = ZCBClosedFormG2(r0, sigma, 1, t, S, a, b, x0, y0, phi0, eta, rho);
    double p_t_t = ZCBClosedFormG2(r0, sigma, 1, t, T, a, b, x0, y0, phi0, eta, rho);
    double sig_left = sigma*sigma/2/a/a/a*(1-exp(-a*(S-T)))*(1-exp(-a*(S-T)))*(1-exp(-2*a*(T-t)));
    double sig_mid = eta*eta/2/b/b/b*(1-exp(-b*(S-T)))*(1-exp(-b*(S-T)))*(1-exp(-2*b*(T-t)));
    double sig_right = 2*rho*sigma*eta/a/b/(a+b)*(1-exp(-a*(S-T)))*(1-exp(-b*(S-T)))*(1-exp(-(a+b)*(T-t)));
    double Sig = sqrt(sig_left+sig_mid+sig_right);
    
    double d1 = log(k*p_t_t/faceValue/p_t_s)/Sig - 0.5*Sig;
    double d2 = log(k*p_t_t/faceValue/p_t_s)/Sig + 0.5*Sig;
    double price = -p_t_s*faceValue*cdfNorm(d1) + p_t_t*k*cdfNorm(d2);
    return price;
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
