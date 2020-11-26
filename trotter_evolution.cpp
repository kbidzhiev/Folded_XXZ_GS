//I would like to give MPS, trotter order and have MPS


//additional parameters:
//args, Ham.ampo, tau

#pragma once
#include "itensor/all.h"
#include <iostream>
#include <complex>

using namespace itensor;


//exponentiation of H
double tau = param.val("tau");
const int o = param.val("TrotterOrder");

auto args = Args("Method=", "Fit", "Cutoff", param.val("trunc"), "MaxDim",
		param.longval("max_bond"), "Normalize", false); // for FitApplyMPO



MPS TimeEvol(AutoMPO  Ham, const int o, double tau){
	MPO expH1, expH2, expH3, expH4, expH5, expH6, expH7; //7 possible complex valued time steps
	Cplx t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0, t7 = 0;
	bool ok = false; //flag. is 1 when toExpH worked
	if (o == 1) {
		//Approx. with error O(tau^3)
		t1 = Cplx_i * tau;
		expH1 = toExpH(Ham, t1);
		ok = true;
	}
	if (o == 2) {
		//Approx. with error O(tau^3)
		t1 = 0.5 * (1 + Cplx_i) * tau;
		t2 = 0.5 * (-1 + Cplx_i) * tau;	
		expH1 = toExpH(Ham, t1);
		expH2 = toExpH(Ham, t2);
		ok = true;
	}
	if (o == 3) {
		//Approx. with error O(tau^4) [thanks to Kemal]
		//Tested -it's ok
		t1 = .10566243270259355887 - .39433756729740644113 * Cplx_i;
		t2 = Cplx_i * t1;
		t3 = conj(t2);
		t4 = Cplx_i * t3;
		t1 *= Cplx_i * tau;
		t2 *= Cplx_i * tau;
		t3 *= Cplx_i * tau;
		t4 *= Cplx_i * tau;
		expH1 = toExpH(Ham, t1);
		expH2 = toExpH(Ham, t2);
		expH3 = toExpH(Ham, t3);
		expH4 = toExpH(Ham, t4);
		ok = true;
	}
	if (o == 4) {
		//Approx. with error O(tau^5) [thanks to Kemal twice]
		//Tested -it's ok
		t1 = 0.25885339861091821723 + 0.04475613401114190287 * Cplx_i;
		t2 = -0.03154685814880379274 + 0.24911905427556321757 * Cplx_i;
		t3 = 0.19082905211066719664 - 0.23185374923210605447 * Cplx_i;
		t4 = 0.1637288148544367438753;
		t5 = conj(t3);
		t6 = conj(t2);
		t7 = conj(t1);
		t1 *= Cplx_i * tau;
		t2 *= Cplx_i * tau;
		t3 *= Cplx_i * tau;
		t4 *= Cplx_i * tau;
		t5 *= Cplx_i * tau;
		t6 *= Cplx_i * tau;
		t7 *= Cplx_i * tau;
		expH1 = toExpH(Ham, t1);
		expH2 = toExpH(Ham, t2);
		expH3 = toExpH(Ham, t3);
		expH4 = toExpH(Ham, t4);
		expH5 = toExpH(Ham, t5);
		expH6 = toExpH(Ham, t6);
		expH7 = toExpH(Ham, t7);
		ok = true;
	}
	if (ok == false)
		cout << "Error, TrotterOrder=" << o << " not implemented.\n", exit(
				1);
	cout << "exp(iHt) constructed.\n";
	cout.flush();
}

