#include "itensor/all.h"
#include "time_evolution.hpp"

using namespace itensor;
using namespace std;




ThreeSiteHamiltonian::ThreeSiteHamiltonian(const SiteSet &sites,
		const ThreeSiteParam &param) :
		ampo(sites), N(length(sites)) {
	//N = length(sites); // size of Hamiltonian comes with onject SiteSet sites, not just as number N
	init(param);   // initializing the Hamiltonian
	cout << "A Hamiltonian with " << N << " sites was constructed." << endl;
}

void ThreeSiteHamiltonian::init(const ThreeSiteParam &param) {   //.init (param)
	const double J = param.val("J");
	double mu = 0;
	const double hL = param.val("hL");
	const double hR = param.val("hR");
	dot = N / 2 + 1;  //Position of the "dot"
	//cout << "The dot is on site #" << dot << endl;
	//if ((2*N)<=3) cout<<"Error, N="<<N<<" is too small.\n",exit(0);
	for (int j = 1; j < N - 1; ++j) {
		//Strange coefficients are needed to match with
		// Pauli matrices instead of spin Sx Sy
		ampo += J * 4 * 0.25, "S+", j, "S-", j + 2; // 0.5 (SpSm+ SmSp) = SxSx + SySy
		ampo += J * 4 * 0.25, "S-", j, "S+", j + 2;
		ampo += J * -8 * 0.25, "S+", j, "Sz", j + 1, "S-", j + 2;
		ampo += J * -8 * 0.25, "S-", j, "Sz", j + 1, "S+", j + 2;
		//cout << "j = "<< j << "/ " << N-3 << endl;
		//cout << "site (" << j << ", " << j + 1 << ", " << j + 2 << ")"
		//		<< endl;
		if (j <= dot) {
			mu = hL;
		} else {
			mu = hR;
		}
		ampo += mu, "Sz", j;
	}
	ampo += mu, "Sz", N - 1;
	ampo += mu, "Sz", N;
}

LadderHamiltonian::LadderHamiltonian(const SiteSet &sites,
		const ThreeSiteParam &param, const string ham_type_) :
		ampo(sites), N(length(sites)), ham_type(ham_type_) {
	init(param);   // initializing the Hamiltonian
	cout << "A LADDER Hamiltonian with " << N << " sites was constructed."
			<< endl;
}

void LadderHamiltonian::init(const ThreeSiteParam &param) {    //.init (param)
	const double J = param.val("J");
	auto h = [&](double j) {
		return param.val("h");
		//return 2.0*cos(j); // this term add "chaotic" h
	};
	const double rho = param.val("rho");
	double m = 2.0;

	dot = N / 2 + 1;  //Position of the central spin

	if (ham_type == "Ladder") {
		for (int j = 1; j <= N - 2; j += 2) {
			// The coefficients 4 and 2 are needed to match between
			// spin Pauli matrices instead of Sx Sy
			ampo += -J * m * 2, "Sz", j + 1; //Sublattice A even sites

			ampo += -J * 4, "Sz", j, "Sz", j + 2; //Sublattice B
			ampo += -J * (h(j) + pow(-1, j / 2) * rho) * 2, "Sx", j; //Sublattice B

//				cout << "Sz : " << j << "\t SxSx : (" << j + 1 << ", " << j + 3
//						<< ")" << "\t h term ("
//						<< (h(j+1) + pow(-1, (j + 1) / 2) * rho) << "): " << j + 1
//						<< endl;
		}
		cout << "Sz : " << N - 1 << "\t \t \t \t h term ("
				<< (h(N) + pow(-1, N / 2) * rho) << "): " << N << endl;
		ampo += -J * m * 2, "Sz", N - 1; // Sublattice A even sites

		ampo += -J * (h(N) + pow(-1, N / 2) * rho) * 2, "Sx", N;
	} else if (ham_type == "Ising") {
		// to create an initial state as GS of ladder ham, we start with GS
		// of uniform ising model to introduce initial correlations,
		// otherwise DMRG procedure cannot find a GS
		for (int j = 1; j < N; j++) {
			//ampo += -J * 4, "Sx", j, "Sx", j+1;
			ampo += -J, "S+", j, "S-", j + 1;
			ampo += -J, "S-", j, "S+", j + 1;

			//ampo += -J * m * 2, "Sz", j;
		}
		//ampo += -J * m * 2, "Sz", N;

	} else {
		throw invalid_argument(
				"One should choose Ladder or Ising in the LadderHamiltonian initialization");
	}

}


