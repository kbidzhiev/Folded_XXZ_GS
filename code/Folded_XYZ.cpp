#include "itensor/all.h"
#include "time_evolution.hpp"

using namespace itensor;
using namespace std;


//https://arxiv.org/abs/2110.11322 in SIGMA repr Eq(1)
HamiltonianFoldedXYZ::HamiltonianFoldedXYZ(const SiteSet &sites,
		const ThreeSiteParam &param)
			: ampo(sites)
			, N(length(sites)) {
	//N = length(sites); // size of Hamiltonian comes with onject SiteSet sites, not just as number N
	init(param);   // initializing the Hamiltonian
	cout << "A Hamiltonian with " << N << " sites was constructed." << endl;
}

//Initialize Hamiltonian parameters
void HamiltonianFoldedXYZ::init(const ThreeSiteParam &param) {
	const double Jx = param.val("J");
	const double Jy = param.val("Jy");
	const double Delta = param.val("Delta");
	const double J2 = param.val("J2");

	for (int j = 1; j < N - 1; ++j) {
		//XYZ in sigma basis
		ampo += Jx * 4, "Sx", j, "Sx", j + 2;
		ampo += -Jy * 8, "Sx", j, "Sz", j + 1, "Sx", j + 2;
		ampo += Delta * 2, "Sz", j;
		// integrability breaking term
		ampo += J2 * 4, "Sz", j, "Sz", j + 1;
	}
	// boundary terms. "for" loop doesn't reach j == N-1 and N
	ampo += Delta * 2, "Sz", N - 1;
	ampo += Delta * 2, "Sz", N;
	ampo += J2 * 4, "Sz", N-1, "Sz", N;
}




