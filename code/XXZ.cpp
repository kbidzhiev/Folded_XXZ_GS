#include "itensor/all.h"
#include "time_evolution.hpp"

using namespace itensor;
using namespace std;



//_____________________________________________________

//
XXZ::XXZ(const SiteSet &sites,
		const ThreeSiteParam &param)
			: ampo(sites)
			, N(length(sites)) {
	//N = length(sites); // size of Hamiltonian comes with onject SiteSet sites, not just as number N
	init(param);   // initializing the Hamiltonian
	cout << "A Hamiltonian with " << N << " sites was constructed." << endl;
}

//Initialize Hamiltonian parameters
void XXZ::init(const ThreeSiteParam &param) {
	const double J = param.val("J");
	double mu = 0;
	const double hL = param.val("hL");
	const double hR = param.val("hR");

	const double Delta = param.val("Delta");
	dot = N / 2 + 1;  //Position of the "dot"
	cout << "The dot is on site #" << dot << endl;
	for (int j = 1; j < N - 1; ++j) {
		//Strange coefficients are needed to match
		// spin matrices and Pauli matrices -> Pauli = 2*Spin, so each matrix gives factor 2
		// one of 1/2 comes from the projector (1-sigma_z)/2
		// and the other is "Jacobian", i.e. 0.5 (SpSm+ SmSp) = SxSx + SySy
		ampo += J * 4 * 0.5, "S+", j, "S-", j + 1; //
		ampo += J * 4 * 0.5, "S-", j, "S+", j + 1;
		ampo += J * 4 * Delta, "Sz", j, "Sz", j + 1;
		if (j <= dot) {
			mu = hL;
		} else {
			mu = hR;
		}
		ampo += mu, "Sz", j;
	}
	// boundary terms. for loop doesn't reach j == N-1 OR N
	ampo += mu, "Sz", N - 1;
	ampo += mu, "Sz", N;
}


//
XP_Hamiltonian::XP_Hamiltonian(const SiteSet &sites,
		const ThreeSiteParam &param)
			: ampo(sites)
			, N(length(sites)) {
	//N = length(sites); // size of Hamiltonian comes with onject SiteSet sites, not just as number N
	init(param);   // initializing the Hamiltonian
	cout << "A Hamiltonian with " << N << " sites was constructed." << endl;
}

//Initialize Hamiltonian parameters
void XP_Hamiltonian::init(const ThreeSiteParam &param) {
	const double J = param.val("J");
	dot = N / 2 + 1;  //Position of the "dot"
	cout << "The dot is on site #" << dot << endl;
	for (int j = 1; j < N - 1; ++j) {
		//Strange coefficients are needed to match
		// spin matrices and Pauli matrices -> Pauli = 2*Spin, so each matrix gives factor 2
		// one of 1/2 comes from the projector (1-sigma_z)/2
		// and the other is "Jacobian", i.e. 0.5 (SpSm+ SmSp) = SxSx + SySy


		ampo += J * 2 * 0.5, "Sx", j;
		ampo += J * 2 * 0.5, "Sx", j + 1;
		ampo += J * 4 * 0.5, "Sx", j, "Sz", j + 1;
		ampo += J * 4 * 0.5, "Sz", j, "Sx", j + 1;


//		ampo += J * 2 * 0.25, "Sx", j; //
//		ampo += J * 4 * 0.25, "Sx", j, "Sz", j + 1;
//		ampo += J * 4 * 0.25, "Sz", j - 1, "Sx", j;
//		ampo += J * 8 * 0.25, "Sz", j - 1, "Sx", j , "Sz", j + 1;
	}
}








