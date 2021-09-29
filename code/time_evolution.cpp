#include "itensor/all.h"
#include "time_evolution.hpp"

using namespace itensor;
using namespace std;



//_____________________________________________________

//Trotter Gates
TrotterExp::TrotterExp(const SiteSet &sites, const ThreeSiteParam &param,
		const complex<double> tau) {
	initialize(sites, param, tau);
}
void TrotterExp::initialize(const SiteSet &sites, const ThreeSiteParam &param,
		const complex<double> tau) {
	//const int begin = param.val("begin");
	const int begin = 1;
	const int end = param.val("N");
	const int order = param.val("TrotterOrder");
	if (order == 1) {
		cout << "trotter 1 scheme" << endl;

		TimeGates(begin, end, tau, sites, param);
		TimeGates(begin + 1, end, tau, sites, param);
		TimeGates(begin + 2, end, tau, sites, param);

	} else {
		cout << "trotter 2 scheme" << endl;
		/*
		 double a1 = 1. / 6;		// more precise arrpoximation coefficients
		 double a2 = 1 - 2. * a1;
		 double b1 = (3 - sqrt(3)) / 6.;
		 double b2 = 1. / 2 - b1;
		 double c1 = 1. / 2;
		 */
		double begin0 = begin; //this variable are needed to change operators ABC
		double begin2 = begin + 1;
		double begin4 = begin + 2;
		//Trotter gates from arxiv.org/abs/1901.04974
		// Eq. (38),(47)

		cout << "Time evolutions " << endl;
		TimeGates(begin0, end, 0.5 * tau, sites, param); //A
		TimeGates(begin2, end, 0.5 * tau, sites, param); //B
		TimeGates(begin4, end, tau, sites, param); //C
		TimeGates(begin2, end, 0.5 * tau, sites, param); //B
		TimeGates(begin0, end, 0.5 * tau, sites, param); //A
		/*
		 TimeGates(begin0, end, a1 * tau, sites, param); //A
		 TimeGates(begin2, end, b1 * tau, sites, param); //B
		 TimeGates(begin4, end, c1 * tau, sites, param); //C
		 TimeGates(begin2, end, b2 * tau, sites, param); //B
		 TimeGates(begin0, end, a2 * tau, sites, param); //A
		 TimeGates(begin2, end, b2 * tau, sites, param); //B
		 TimeGates(begin4, end, c1 * tau, sites, param); //C
		 TimeGates(begin2, end, b1 * tau, sites, param); //B
		 TimeGates(begin0, end, a1 * tau, sites, param); //A
		 */
	}
}

void TrotterExp::TimeGates(const int begin, const int end,
		const complex<double> tau, const SiteSet &sites,
		const ThreeSiteParam &param) {
	const int step = 3;
	const double J = param.val("J");
	const double Dhar = param.val("Dhar");
	const double PXXP = param.val("PXXP");
	//cout << "Gates starts from " << begin << endl;
	for (int j = begin; j < end - 1; j += step) {
		//cout << "j = (" << j << ", " << j + 1 << ", " << j + 2 << ")"
		//		<< endl;
		//this part act on real sites
		auto hh = J * 4 * 0.25 * op(sites, "Sp", j) * op(sites, "Id", j + 1)
				* op(sites, "Sm", j + 2);
		hh += J * 4 * 0.25 * op(sites, "Sm", j) * op(sites, "Id", j + 1)
				* op(sites, "Sp", j + 2);
		hh += -J * 8 * 0.25 * op(sites, "Sp", j) * op(sites, "Sz", j + 1)
				* op(sites, "Sm", j + 2);
		hh += -J * 8 * 0.25 * op(sites, "Sm", j) * op(sites, "Sz", j + 1)
				* op(sites, "Sp", j + 2);

		//Deepak Dhar term
		if (Dhar > 0) {
			cout << "Dhar term is included" << endl;
			hh += Dhar * 4 * 0.5 * op(sites, "Sz", j) * op(sites, "Id", j + 1)
					* op(sites, "Sz", j + 2);
			hh += -Dhar * 8 * 0.5 * op(sites, "Sz", j) * op(sites, "Sz", j + 1)
					* op(sites, "Sz", j + 2);
			hh += -Dhar * 1 * 0.5 * op(sites, "Id", j) * op(sites, "Id", j + 1)
					* op(sites, "Id", j + 2);
			hh += Dhar * 2 * 0.5 * op(sites, "Id", j) * op(sites, "Sz", j + 1)
					* op(sites, "Id", j + 2);
		}
		if (PXXP > 0) {
			cout << "PXXP = (1-Sz)(1-Sz)(Sx+Sx)(1-Sz)(1-Sz)"
					<< "term is included " << endl;
			auto P =
					0.25
							* (op(sites, "Id", j) * op(sites, "Id", j + 1)
									+ 2.0 * op(sites, "Sz", j)
											* op(sites, "Id", j + 1))
							* (op(sites, "Id", j) * op(sites, "Id", j + 1)
									+ 2.0 * op(sites, "Id", j)
											* op(sites, "Sz", j + 1));

			hh += PXXP * P * 2.0
					* (op(sites, "Sx", j) * op(sites, "Id", j + 1)
							* op(sites, "Id", j + 2)
							+ op(sites, "Id", j) * op(sites, "Sx", j + 1)
									* op(sites, "Id", j + 2)) * P;
		}

		auto G = expHermitian(hh, tau);
		gates.emplace_back(j, move(G));
	}
}
void TrotterExp::Evolve(MPS &psi, const Args &args) {
	for (auto &gate : gates) {
		auto j = gate.i1;
		auto &G = gate.G;
		psi.position(j);
		auto WF = psi(j) * psi(j + 1) * psi(j + 2);
		WF = G * WF;
		WF /= norm(WF);
		WF.noPrime();
		{
			auto [Uj1, Vj1] = factor(WF,
					{ siteIndex(psi, j), leftLinkIndex(psi, j) }, args);
			auto indR = commonIndex(Uj1, Vj1);
			auto [Uj2, Vj2] = factor(Vj1, { siteIndex(psi, j + 1), indR }, args);
			psi.set(j, Uj1);
			psi.set(j + 1, Uj2);
			psi.set(j + 2, Vj2);

		}
	}
}

vector<MPO> XXZ_time_evol(const SiteSet &sites, const ThreeSiteParam &param) {

	double tau = param.val("tau");
	const int order = param.val("TrotterOrder");

	vector<complex<double>> time_steps;
	time_steps.reserve(order);
	if (order == 1) {
		//Approx. with error O(tau^3)
		time_steps[0] = Cplx_i * tau;
	}
	if (order == 2) {
		//Approx. with error O(tau^3)
		time_steps[0] = 0.5 * ( 1 + Cplx_i) * tau;
		time_steps[1] = 0.5 * (-1 + Cplx_i) * tau;

	}
	if (order == 3) {
		//Approx. with error O(tau^4) [thanks to Kemal]
		//Tested -it's ok
		time_steps[0] = 0.10566243270259355887 - 0.39433756729740644113 * Cplx_i;
		time_steps[1] = Cplx_i * time_steps[0];
		time_steps[2] = conj(time_steps[1]);
		time_steps[3] = Cplx_i * time_steps[2];
		for (auto &time_step : time_steps){
			time_step *= Cplx_i * tau;
		}
	}
	if (order == 4) {
		//Approx. with error O(tau^5) [thanks to Kemal twice]
		//WARNINIG Not fully TESTED
		//Tested -it's ok
		time_steps[0] = 0.25885339861091821723 + 0.04475613401114190287 * Cplx_i;
		time_steps[1] = -0.03154685814880379274 + 0.24911905427556321757 * Cplx_i;
		time_steps[2] = 0.19082905211066719664 - 0.23185374923210605447 * Cplx_i;
		time_steps[3] = 0.1637288148544367438753;
		time_steps[4] = conj(time_steps[2]);
		time_steps[5] = conj(time_steps[1]);
		time_steps[6] = conj(time_steps[0]);
		for (auto &time_step : time_steps){
			time_step *= Cplx_i * tau;
		}

	}

	vector<MPO> Exp_H_vec;
	Exp_H_vec.reserve(order);

	XXZ Init_H(sites, param);
	for (auto time_step : time_steps){
		MPO expH_i = toExpH(Init_H.ampo, time_step);
		Exp_H_vec.push_back(expH_i);
	}
	return Exp_H_vec;
}



