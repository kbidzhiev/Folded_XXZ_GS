#include "itensor/all.h"
#include "time_evolution.hpp"
#include <cmath>

using namespace itensor;
using namespace std;


//
ThreeSiteHamiltonian::ThreeSiteHamiltonian(const SiteSet &sites,
		const ThreeSiteParam &param)
			: ampo(sites)
			, N(length(sites)) {
	//N = length(sites); // size of Hamiltonian comes with onject SiteSet sites, not just as number N
	init(param);   // initializing the Hamiltonian
	cout << "A Hamiltonian with " << N << " sites was constructed." << endl;
}

//Initialize Hamiltonian parameters
void ThreeSiteHamiltonian::init(const ThreeSiteParam &param) {
	const double J = param.val("J");
	double mu = 0;
	const double hL = param.val("hL");
	const double hR = param.val("hR");
	const int PPK = sqrt(2.0)*param.val("PPK");
	const int PPX = sqrt(2.0)*param.val("PPX");
	int shift_of_range = 0;
	if (PPK != 0) shift_of_range = 1;
	dot = N / 2 + 1;  //Position of the "dot"
	cout << "The dot is on site #" << dot << endl;
	for (int j = 1; j < N - 1 - shift_of_range; ++j) {
		//Strange coefficients are needed to match
		// spin matrices and Pauli matrices -> Pauli = 2*Spin, so each matrix gives factor 2
		// one of 1/2 comes from the projector (1-sigma_z)/2
		// and the other is "Jacobian", i.e. 0.5 (SpSm+ SmSp) = SxSx + SySy
		ampo += J *  4 * 0.25, "S+", j, "S-", j + 2; //
		ampo += J *  4 * 0.25, "S-", j, "S+", j + 2;
		ampo += J * -8 * 0.25, "S+", j, "Sz", j + 1, "S-", j + 2;
		ampo += J * -8 * 0.25, "S-", j, "Sz", j + 1, "S+", j + 2;

		if (PPK != 0) {
			ampo +=  PPK * 0.5, "S+", j + 2, "S-", j + 3; // (I)
			ampo +=  PPK * 0.5, "S-", j + 2, "S+", j + 3;
			ampo += -PPK * 2 * 0.5, "Sz", j, "S+", j + 2, "S-", j + 3; // (II)
			ampo += -PPK * 2 * 0.5, "Sz", j, "S-", j + 2, "S+", j + 3;
			ampo += -PPK * 2 * 0.5, "Sz", j + 1, "S+", j + 2, "S-", j + 3; // (III)
			ampo += -PPK * 2 * 0.5, "Sz", j + 1, "S-", j + 2, "S+", j + 3;
			ampo +=  PPK * 4 * 0.5, "Sz", j, "Sz", j + 1, "S+", j + 2, "S-", j + 3; // (IV)
			ampo +=  PPK * 4 * 0.5, "Sz", j, "Sz", j + 1, "S-", j + 2, "S+", j + 3;
		}

		if (PPX != 0) {
			ampo +=  PPX * 0.5, "Sx", j + 2; // (I)
			ampo +=  -PPX , "Sz", j, "Sx", j + 2;
			ampo +=  -PPX , "Sz", j + 1, "Sx", j + 2;
			ampo +=  PPX * 2.0, "Sz", j, "Sz", j + 1, "Sx", j + 2;
		}


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


LadderHamiltonian::LadderHamiltonian(const SiteSet &sites,
		const ThreeSiteParam &param, const string ham_type_)
			: ampo(sites)
			, N(length(sites))
			, ham_type(ham_type_) {
				init(param);   // initializing the Hamiltonian
				cout << "A LADDER Hamiltonian with "
						<< N << " sites was constructed."
						<< endl;
}

void LadderHamiltonian::init(const ThreeSiteParam &param) {    //.init (param)
	const double J = param.val("J");
	const double magnetic_field = param.val("h");
	const double rho = param.val("rho");
	const double m = 2.0;
	dot = N / 2 + 1;  //Position of the central spin

	if (ham_type == "Ladder") {
		for (int j = 1; j <= N - 2; j += 2) {
			//Strange coefficients are needed to match
			// spin matrices and Pauli matrices -> Pauli = 2*Spin, so each matrix gives factor 2
			ampo += -J * m * 2, "Sz", j + 1;		 //Sublattice A even sites
			ampo += -J * 4, "Sz", j, "Sz", j + 2; 	 //Sublattice B
			ampo += -J * (magnetic_field + pow(-1, j / 2) * rho) * 2, "Sx", j; //Sublattice B
		}
		ampo += -J * m * 2, "Sz", N - 1; // Sublattice A even sites

		ampo += -J * (magnetic_field + pow(-1, N / 2) * rho) * 2, "Sx", N;
	} else if (ham_type == "Ising") {
		// to create an initial state as a GroundState of ladder ham,
		// we start with GS of uniform Ising model to introduce initial correlations,
		// otherwise DMRG procedure cannot find a GS
		for (int j = 1; j < N-1; j += 2) {
			ampo += -J, "S+", j, "S-", j + 2;
			ampo += -J, "S-", j, "S+", j + 2;
			//ampo += -J, "Sz", j, "Sz", j + 2;
			ampo += 5.5*J, "Sx", j;

			//ampo += -J * m * 2, "Sz", j;

			if (j == N - 2) ampo += 5.5*J, "Sx", N ;
			if (j == N - 3) ampo += 5.5*J, "Sx", N - 1;
		}


	} else {
		throw invalid_argument(
				"One should choose 'Ladder' or 'Ising' in the LadderHamiltonian ham_type_ parameter");
	}
}


