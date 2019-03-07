//
//  justinComp.hpp
//  405HW2
//
//  Created by Justin Tan on 1/19/19.
//  Copyright © 2019 UCLA MFE. All rights reserved.
//

#pragma once

#include <stdio.h>
#include <vector>
#include <math.h>

using namespace std;

//hw1 methods
vector<double> randUniformGen(int num, int seed);
vector<double> normGenBM(vector<double> uniform);


double calcMean(vector<double> v);
double calcStdDev(vector<double> v);

