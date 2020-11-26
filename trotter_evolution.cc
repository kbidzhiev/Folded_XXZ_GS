//I would like to give MPS, trotter order and have MPS


//additional parameters:
//args, Ham.ampo, tau

#pragma once
#include "itensor/all.h"
#include <iostream>
#include <complex>

using namespace itensor;
using namespace std;

//exponentiation of H
//double tau = param.val("tau");
//const int o = param.val("TrotterOrder");

//auto args = Args("Method=", "Fit", "Cutoff", param.val("trunc"), "MaxDim",
//param.longval("max_bond"), "Normalize", false); // for FitApplyMPO

void PrintKemal(){
	cout << "Kemal"<< endl;
}
