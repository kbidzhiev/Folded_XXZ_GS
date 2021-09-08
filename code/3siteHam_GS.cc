#include "itensor/all.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
#include <sstream>
#include <complex>
#include <math.h>
#include <chrono>
#include <cmath>
#include <stdexcept>


#include <cstdlib>
#include "observables_GS.h"
#include "quantum_measurement.h"

using namespace itensor;
using namespace std;
using namespace std::chrono;

//------------------------------------------------------------------
//The function below translate numbers (etc.) into character strings
//the second parameter (optional) is the precision (digits)

template<class T>
inline string to_string(const T &t, unsigned int precision = 0) {
	stringstream ss;
	if (precision > 0)
		ss.precision(precision);
	ss << t;
	return ss.str();
}
//____________________________________________________________
double char2double(char *a) {
	char *end_ptr;
	const double x = strtod(a, &end_ptr);
	if (end_ptr == a || ('\0' != *end_ptr))
		cout << endl << "ERROR :" << a << " is not a valid format for a double."
				<< endl, exit(0);
	return x;
}
//____________________________________________________________
class Parameters: public map<string, double> { // class Parameters have all methods from container "map" (string key, double value) 
public:
	double val(string var_name) const { // .val("C") gives value for parmeter C
		map<string, double>::const_iterator it = find(var_name);
		if (it == end()) {
			cout << "Error: Parameter " << var_name << " is not defined.\n", exit(
					0);
			return 0;
		} else
			return it->second;
	}
	long longval(string var_name) const {
		double v = val(var_name);
		if (abs(double(round(v)) - v) < 1e-6) {
			return long(round(v));
		} else {
			cout << "Error, parameter " << var_name << "=" << v
					<< " is not a long" << endl, exit(0);
			return 0;
		}
	}
	void PRint(ostream &o) const {
		for (map<string, double>::const_iterator it = begin(); it != end();
				it++) {
			o << it->first << "=" << it->second << endl;
		}
	}
	void ReadArguments(int argc, char *argv[]) {
		for (int n = 1; n < argc; n++) {
			string var_name(argv[n]);
			map<string, double>::const_iterator it = find(var_name);

			if (it != end()) {
				n++;
				if (n == argc)
					cerr << "Error: missing value after " << var_name << endl, exit(
							0);
				operator[](var_name) = char2double(argv[n]);
			} else {
				cerr << "Syntax error :" << var_name << endl;
				cout << "List of command-line parameters :";
				PRint (cout);
				exit(0);
			}
		}
	}
};
//_____________________________________________________

// class of PUBLIC parameters
class ThreeSiteParam: public Parameters {
public:
	ThreeSiteParam() { //Constructor
		//Specify below all the allowed parameter names,
		//and their default values
		operator[]("N") = 10; //Length of the chain
		operator[]("J") = 1.0;
		operator[]("tau") = 0.02;  //time step for the unitary evolution
		operator[]("T") = 2;  //Total (final) time
		operator[]("Sz") = 0.1;
		operator[]("SVD_spec") = 0; //SVD spectrum
		operator[]("max_bond") = 4000;  //maximum bond dimension
		operator[]("trunc") = 1e-8;  //maximum truncation error
		operator[]("energy") = 1e-13;  //convergence criterium on the energy
		operator[]("sweeps") = 999;  //maximum number of sweeps in the DMRG
		operator[]("TrotterOrder") = 2;
		operator[]("GroundState") = 0;
		operator[]("LadderState") = 0;
		operator[]("JammedImpurity") = 0;
		operator[]("alpha") = 0; // U = exp[ i alpha  n*\sigma]
		operator[]("Theta") = 0; // U = exp[ i theta  n*\sigma] turns only the central site
		operator[]("JammedShift") = 0;
		operator[]("DoubleSlit") = 0;
		operator[]("Sav") = 0;
		operator[]("Lenart") = 0;
		operator[]("begin") = 1;
		//operator[]("end") = 10;
		operator[]("hL") = 0; //initial magnetization_L
		operator[]("hR") = 0; //initial magnetization_R
		operator[]("h") = 2.0;
		operator[]("rho") = 0.0;
		operator[]("n") = 1;
		operator[]("Q1Profile") = 0; // energy and current profile
		operator[]("Q2Profile") = 0;
		operator[]("Current") = 0;
		operator[]("Entropy") = 0; //entanglement entropy p*log*p between left and right parts of system
		operator[]("EntropyProfile") = 0; // Entropy Profile- parameter 0 -> nothing, dt>0 each second=integer parameter
		operator[]("Loschmidt") = 0; // loschmidt echo <psi(t)|psi(0)>
		operator[]("Dhar") = 0; // Deepak Dhar term in hamiltonian (time evolution ONLY)
		operator[]("PXXP") = 0; //Integrability breaking term
		operator[]("Measurement") = 0; // Project the central spin to be |Down> (time evolution ONLY)
		operator[]("AlphaGate") = 0;
	}

};
//_____________________________________________________

//I'm creating a 3 site Hamiltonian for system of N sites

class ThreeSiteHamiltonian {
public:
	int dot;
	AutoMPO ampo;  //variable ampo
	//Define a constructor for this class 'ThreeSiteHamiltonian'
	ThreeSiteHamiltonian(const SiteSet &sites, const ThreeSiteParam &param) :
			ampo(sites), N(length(sites)) {
		//N = length(sites); // size of Hamiltonian comes with onject SiteSet sites, not just as number N
		init(param);   // initializing the Hamiltonian
		cout << "A Hamiltonian with " << N << " sites was constructed." << endl;
	}
private:
	int N;
	void init(const ThreeSiteParam &param) {    //.init (param)
		const double J = param.val("J");
		double mu = 0;
		const double hL = param.val("hL");
		const double hR = param.val("hR");
		dot = N / 2 + 1;  //Position of the "dot"
		//cout << "The dot is on site #" << dot << endl;
		//if ((2*N)<=3) cout<<"Error, N="<<N<<" is too small.\n",exit(0);
		for (int j = 1; j < N-1; ++j) {
			//Strange coefficients are needed to match with
			// Pauli matrices instead of spin Sx Sy
			ampo += J * 4 * 0.25, "S+", j, "S-", j + 2; // 0.5 (SpSm+ SmSp) = SxSx + SySy
			ampo += J * 4 * 0.25, "S-", j, "S+", j + 2;
			ampo += J * -8 * 0.25, "S+", j, "Sz", j + 1, "S-", j + 2;
			ampo += J * -8 * 0.25, "S-", j, "Sz", j + 1, "S+", j + 2;
			//cout << "j = "<< j << "/ " << N-3 << endl;
			//cout << "site (" << j << ", " << j + 1 << ", " << j + 2 << ")"
			//		<< endl;
			if(j <= dot){
				mu = hL;
			}else{
				mu = hR;
			}
			ampo += mu, "Sz", j;
		}
		ampo += mu, "Sz", N-1;
		ampo += mu, "Sz", N;
	}

};

class LadderHamiltonian {
public:
	int dot;
	AutoMPO ampo;

	//Define a constructor for this class 'LadderHamiltonian'
	LadderHamiltonian(const SiteSet &sites, const ThreeSiteParam &param, const string ham_type_)
		: ampo(sites)
		, N(length(sites))
		, ham_type(ham_type_){
		init(param);   // initializing the Hamiltonian
		cout << "A LADDER Hamiltonian with " << N << " sites was constructed." << endl;
	}
private:
	int N;
	string ham_type;
	void init(const ThreeSiteParam &param) {    //.init (param)
		const double J = param.val("J");
		auto h = [&](double j){
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

				ampo += -J * 4, "Sz", j , "Sz", j + 2; //Sublattice B
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
		} else if (ham_type == "Ising"){
			// to create an initial state as GS of ladder ham, we start with GS
			// of uniform ising model to introduce initial correlations,
			// otherwise DMRG procedure cannot find a GS
			for (int j = 1; j < N; j++){
				//ampo += -J * 4, "Sx", j, "Sx", j+1;
				ampo += -J , "S+", j, "S-", j+1;
				ampo += -J , "S-", j, "S+", j+1;


				//ampo += -J * m * 2, "Sz", j;
			}
			//ampo += -J * m * 2, "Sz", N;

		} else {
			throw invalid_argument("One should choose Ladder or Ising in the LadderHamiltonian initialization");
		}

	}
};


//Trotter Gates
class TrotterExp {
public:
	struct TGate {
		int i1 = 0;
		ITensor G;
		TGate(int i1_, ITensor G_) :
				i1(i1_), G(G_) {
		}
	};
	TrotterExp(const SiteSet &sites, const ThreeSiteParam &param,
			const complex<double> tau) {
		initialize(sites, param, tau);
	}
	void initialize(const SiteSet &sites, const ThreeSiteParam &param,
			const complex<double> tau) {
		//const int begin = param.val("begin");
		const int begin = 1;
		const int end = param.val("N");
		const int order = param.val("TrotterOrder");
		if (order == 1) {
			cout << "trotter 1 scheme" << endl;

			TimeGates(begin,     end, tau, sites, param);
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
			TimeGates(begin4, end,       tau, sites, param); //C
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

	void TimeGates(const int begin, const int end, const complex<double> tau,
			const SiteSet &sites, const ThreeSiteParam &param) {
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
					<<	"term is included " << endl;
				auto P = 0.25 * (
							op(sites, "Id", j ) * op(sites, "Id", j + 1)
							+ 2.0 * op(sites, "Sz", j ) * op(sites, "Id", j + 1)
						) * (
							op(sites, "Id", j ) * op(sites, "Id", j + 1)
							+ 2.0 * op(sites, "Id", j ) * op(sites, "Sz", j + 1)
							);

				hh += PXXP * P
						* 2.0 *(
							op(sites, "Sx", j ) * op(sites, "Id", j + 1) * op(sites, "Id", j + 2)
							+ op(sites, "Id", j) * op(sites, "Sx", j + 1) * op(sites, "Id", j + 2)
						)
						* P;
			}


			auto G = expHermitian(hh, tau);
			gates.emplace_back(j, move(G));
		}
	}
	void Evolve(MPS &psi, const Args &args) {
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
				auto [Uj2, Vj2] = factor(Vj1, { siteIndex(psi, j + 1), indR },
						args);
				psi.set(j, Uj1);
				psi.set(j + 1, Uj2);
				psi.set(j + 2, Vj2);

			}
		}
	}

private:
	vector<TGate> gates;
};

//------------------------------------------------------------------
class MyDMRGObserver: public DMRGObserver { // vector<int> ~ DMRGObserver<ITensor>, where ITensor are tensor object type.
	double previous_energy;
	const double precision;
public:
	MyDMRGObserver(const MPS &psi, double prec = 1e-10) :
			DMRGObserver(psi), precision(prec) {
	}
	bool checkDone(const Args &args = Args::global()) { const
		double energy = args.getReal("Energy", 0);
		cout << "    Energy change:" << energy - previous_energy << endl;
		if (abs(energy - previous_energy) < precision) {
			cout << "   Energy has converged -> stop the DMRG.\n";
			return true;
		} else {
			previous_energy = energy;
			return false;
		}
	}
};

// main
int main(int argc, char *argv[]) {

	ThreeSiteParam param;
	param.ReadArguments(argc, argv); //Now param contains the parameters, default values or those provided on the command-line

	param.PRint(cout); // Print parameters
	cout.precision(15);

	const int N = param.longval("N");

	SpinHalf sites(N, { "ConserveQNs=", false }); // Creating a Hilbert space
	MPS psi, psi0;
	auto args = Args("Method=", "DensityMatrix", "Cutoff", param.val("trunc"),
			"MaxDim", param.longval("max_bond"), "Normalize", true); // for FitApplyMPO RENAME


// Preparing an initial state

	ThreeSiteHamiltonian Init_H(sites, param);
	auto H0 = toMPO(Init_H.ampo);
	double energy_initial = 0;

	auto HadamarGate = [&](int i) {
		auto ind = sites(i);
		auto indP = prime(sites(i));
		auto Had = ITensor(ind, indP);
		Had.set(ind(1), indP(1), ISqrt2);
		Had.set(ind(1), indP(2), ISqrt2);
		Had.set(ind(2), indP(1), ISqrt2);
		Had.set(ind(2), indP(2), -ISqrt2);
		psi.setA(i, psi.A(i) * Had);
	};

	auto UnitaryGate = [&](int i, const double alpha) {
		// U = exp[i alpha (n*s)]; |n|^2==1, s = {sx,sy,sz}
		// n = {1,1,1}/sqrt(3);
		const double nx = 1./sqrt(3.0);
		const double ny = nx;
		const double nz = nx;
		auto ind = sites(i);
		auto indP = prime(sites(i));
		auto Had = ITensor(ind, indP);
		Had.set(ind(1), indP(1),  cos(alpha) + Cplx_i * nz * sin(alpha));
		Had.set(ind(1), indP(2), ( ny + Cplx_i * nx) * sin(alpha));
		Had.set(ind(2), indP(1), (-ny + Cplx_i * nx) * sin(alpha));
		Had.set(ind(2), indP(2),  cos(alpha) - Cplx_i * nz * sin(alpha));
		psi.setA(i, psi.A(i) * Had);
	};

	auto AlphaGate = [&](int i, const double alpha) {
		// U = exp[i alpha (n*s)]; |n|^2==1, s = {sx,sy,sz}
		// n = {0,1,0}
		// the spin will be {cos a, sin a}
		auto ind = sites(i);
		auto indP = prime(sites(i));
		auto Had = ITensor(ind, indP);
		Had.set(ind(1), indP(1),  cos(alpha));
		Had.set(ind(1), indP(2),  sin(alpha));
		Had.set(ind(2), indP(1), -sin(alpha));
		Had.set(ind(2), indP(2),  cos(alpha));
		psi.setA(i, psi.A(i) * Had);
	};
	auto SigmaXGate = [&](int i) {
		auto ind = sites(i);
		auto indP = prime(sites(i));
		auto Had = ITensor(ind, indP);
		Had.set(ind(1), indP(1), 0);
		Had.set(ind(1), indP(2), 1);
		Had.set(ind(2), indP(1), 1);
		Had.set(ind(2), indP(2), 0);
		psi.setA(i, psi.A(i) * Had);
	};

	if (param.longval("GroundState") == 1) {
		// GS of the initial Hamiltonian
		cout << "initial state is GS" << endl;
		auto sweeps = Sweeps(999); //number of sweeps is 5
		sweeps.maxdim() = 10, 20, 50, 50, 100, 300, 4000;
		sweeps.cutoff() = 1E-10;

		psi = randomMPS(sites);

		cout << "Before DMRG" << endl;

		cout << "Norm is = " << real(innerC(psi, psi)) << endl;
		energy_initial = real(innerC(psi, H0, psi)); //<psi|H0|psi>
		cout << "1. Initial energy=" << energy_initial << endl;

		MyDMRGObserver obs(psi, param.val("energy"));

		tie(energy_initial, psi) = dmrg(H0, psi, sweeps, obs, "Quiet");

		cout << "After DMRG" << endl;

		int central_site = Init_H.dot;
		//SigmaXGate(central_site );
		//HadamarGate(central_site);
		psi.noPrime();

		//psi0 = psi; //we create two states. Psi for my time evolution, psi0 for standard one

	} else if (param.longval("LadderState") == 1) {
		// GS of the LADDER Hamiltonian
		cout << "initial state is LADDER" << endl;

		LadderHamiltonian Init_H_Ising(sites, param, "Ising");
		auto H_Ising = toMPO(Init_H_Ising.ampo);

		auto sweeps = Sweeps(999); //number of sweeps is 5
		sweeps.maxdim() = 20, 50, 50, 100, 300, 4000;
		sweeps.cutoff() = 1E-10;

//		auto psi0 = randomMPS(sites);
		psi = randomMPS(sites);

	    energy_initial = real(innerC(psi, H_Ising, psi)); //<psi|H0|psi>
		MyDMRGObserver obs(psi, param.val("energy"));
		tie(energy_initial, psi) = dmrg(H_Ising, psi, sweeps, obs, "Quiet");
		cout << "First step is DONE. Ising Ham is ready" << endl;

		LadderHamiltonian Init_H_Ladder(sites, param, "Ladder");
		auto H_Ladder = toMPO(Init_H_Ladder.ampo);
		energy_initial = inner(psi, H_Ladder, psi); //<psi|H0|psi>
		//MyDMRGObserver obs(psi, param.val("energy"));
		tie(energy_initial, psi) = dmrg(H_Ladder, psi, sweeps, obs, "Quiet");

		cout << "Second step is DONE. Ladder Ham is ready" << endl;

		cout << "Norm (before unitary gates )is = " << real(innerC(psi, psi)) << endl;

//		HadamarGate(N/2 - 1);
//		HadamarGate(N/2    );
//		HadamarGate(N/2 + 1);
//		HadamarGate(N/2 + 2);
		int central_site = Init_H_Ladder.dot;
		//SigmaXGate(central_site );
		SigmaXGate(central_site);
		psi.noPrime();


		cout << "Norm (after unitary gates )is = " << real(innerC(psi, psi)) << endl;

	} else if ( param.longval("JammedImpurity") == 1) {
		cout << "initial state is  | Up Left Up Right >  with the flipped spin" << endl;
		auto initState = InitState(sites);
		// Hadamar_2 Hadamar_4 |---+> = |- left - right>
		for (int i = 1; i <= N; ++i){
			if (i % 4 == 0){ // We start counting from 1 ! so the first sites will be |Up Up Up Dn>
				initState.set(i, "Dn");
			} else {
				initState.set(i, "Up");
			}
		}
		//initState.set(N/2 - 1,"Dn");
		//initState.set(N/2,    "Dn");
		//initState.set(N/2 + 1,"Up");
		//initState.set(N/2 + 2,"Up");

		//initState.set(N/2 + 1,"Dn");
		//initState.set(N/2 + 2,"Dn");

		//initState.set(N/2 - 1,"Dn");
		psi = MPS(initState);
		for (int i = 1; i <= N; ++i) {
			if (i % 2 == 0  ) {
				HadamarGate(i);
			}
		}

		//SigmaXGate(N/2+1);
		SigmaXGate(N/2);


//		AlphaGate(N/2-1, 0.5);
//		HadamarGate(N / 2 - 1);
//		HadamarGate(N / 2 	 );
//		HadamarGate(N / 2 + 1);
//		HadamarGate(N / 2 + 2);

//		const double alpha = param.val("alpha");
//		UnitaryGate(N/2 - 1,alpha);
//		UnitaryGate(N/2    ,alpha);
//		UnitaryGate(N/2 + 1,alpha);
//		UnitaryGate(N/2 + 2,alpha);

		psi.noPrime();


	} else if ( param.val("DoubleSlit") > 0) {
		cout << "initial state is  | Up Left Up Right >  with the 2 flipped spins" << endl;
		auto initState = InitState(sites);
		// Hadamar_2 Hadamar_4 |---+> = |- left - right>
		for (int i = 1; i <= N; ++i){
			if (i % 4 == 0){ // We start counting from 1 ! so the first sites will be |Up Up Up Dn>
				initState.set(i, "Dn");
			} else {
				initState.set(i, "Up");
			}
		}

		psi = MPS(initState);
		for (int i = 1; i <= N; ++i) {
			if (i % 2 == 0  ) {
				HadamarGate(i);
			}
		}

		int distance_between_spinflips = param.val("DoubleSlit");

		int to_right = distance_between_spinflips/2 ;
		int to_left = distance_between_spinflips - to_right;
		SigmaXGate(N/2 + to_right);
		SigmaXGate(N/2 - to_left);
		psi.noPrime();

	} else if ( param.val("SingleSlit") > 0) {
		cout << "initial state is  | Up Left Up Right >  with the flipped spin" << endl;
		auto initState = InitState(sites);
		// Hadamar_2 Hadamar_4 |---+> = |- left - right>
		for (int i = 1; i <= N; ++i){
			if (i % 4 == 0){ // We start counting from 1 ! so the first sites will be |Up Up Up Dn>
				initState.set(i, "Dn");
			} else {
				initState.set(i, "Up");
			}
		}

		psi = MPS(initState);
		for (int i = 1; i <= N; ++i) {
			if (i % 2 == 0  ) {
				HadamarGate(i);
			}
		}

		int distance_between_spinflips = param.val("DoubleSlit");

		int to_right = distance_between_spinflips/2 ;
		SigmaXGate(N/2 + to_right);
		psi.noPrime();

	} else if ( param.val("AlphaGate") > 0) {
			auto initState = InitState(sites);
			for (int i = 1; i <= N; ++i){
					initState.set(i, "Up");
			}

			initState.set(N/2 - 1,"Dn");

			psi = MPS(initState);

			const double alpha = param.val("AlphaGate");
			for (int i = 1; i <= N; ++i) {
				if (i % 4 == 0) {
					AlphaGate(i, alpha);
				}else if (i % 4 == 2){
					AlphaGate(i, -alpha);
				}
			}
			psi.noPrime();

	} else if ( param.longval("Sav") == 1) {
			cout << "initial state is  | Up Left Up Left > " << endl;
			auto initState = InitState(sites);

			for (int i = 1; i <= N; ++i){
					initState.set(i, "Up");
			}

			//initState.set(N/2 ,"Dn");
			psi = MPS(initState);
			for (int i = 1; i <= N; ++i) {
				if (i % 2 == 0) {
					HadamarGate(i);
				}
			}
			SigmaXGate(N/2);
			psi.noPrime();
	} else if ( param.val("Theta") > 0 ) {
		cout << "initial state is  | Up Left Up Right > " << endl;
		auto initState = InitState(sites);
		// Hadamar_2 Hadamar_4 |---+> = |- left - right>
		for (int i = 1; i <= N; ++i){
			if (i % 4 == 1){ // We start counting from 1 ! so the first sites will be |Up Up Up Dn>
				initState.set(i, "Dn");
			} else {
				initState.set(i, "Up");
			}
		}
		psi = MPS(initState);
		for (int i = 1; i <= N; ++i) {
			if (i % 2 == 1) {
				HadamarGate(i);
			}
		}
		psi.noPrime();
		//psi0 = psi;

		AlphaGate(N/2+1, (double)param.val("Theta"));
		psi.noPrime();
	}else {
		cout << "Choose: GroundState, Neel, DomainWall,Impurity, Jammed = 1" << endl;
		return 1;
	}
	//cout << "PSI = " << psi << endl;
	psi0 = psi; //here I save the initial WF


//--------------------------------------------------------------

	//Hamiltonian for the dynamics
	ThreeSiteHamiltonian Ham(sites, param);
	const int dot = Ham.dot;
	auto H = toMPO(Ham.ampo);
	cout << "H constructed" << endl;

	energy_initial = real(innerC(psi, H, psi)); //<psi|H0|psi>
	cout << "2. Initial Energy psi =" << energy_initial << endl;

	double norm_initial = real(innerC(psi, psi)); //<psi|H0|psi>
	cout << "2. Initial NORM psi =" << norm_initial << endl;

// Output and observables
// _______
	ofstream ent, spec, entropy_profile, sz, sz_avrg, energy_beta,
				q1_profile, q2_profile, sx, sy, sxsysz,  sx_avrg, loschmidt,
				correlation1, correlation2; //here I'm defining output streams == files
		ios_base::openmode mode;
		mode = std::ofstream::out; //Erase previous file (if present)


	double dt = param.val("Loschmidt");
	if (dt != 0) { //Loschmidt echo
		loschmidt.open("Loschmidt.dat", mode);
		loschmidt.precision(15);
		loschmidt << "#time \t  re<psi(t)|psi(0)> \t re<psi(t)|psi(0)> \n";
	}
	//---------------------
	dt = param.val("Entropy");
	if (dt != 0) { //Entropy in the center of the chain
		ent.open("Entropy_center.dat", mode);
		ent.precision(15);
		ent << "#time \t Entropy(dot) \t BondDim(dot) \t MaxBondDim\n";
	}
	//---------------------
	dt = param.val("SVD_spec");
	if (dt > 0) { //SVD Spectrum on central bond
		spec.open("SVD_spec.dat", mode);
		spec.precision(15);
		spec << "#Position=" << dot << "\t<SVD_spectrum>\t\ttime\n";
	}

	//---------------------
	dt = param.val("EntropyProfile");
	if (dt > 0) { //Full entropy Profile
		entropy_profile.open("Entropy_profile.dat", mode);
		entropy_profile.precision(15);
		entropy_profile << "#Position=i-" << dot << setw(16) << "\t Entropy(i)"
				<< setw(16)
				<< "\t Entropy_sqrt \t Entropy_state1 \t time \t\t Bond.Dim(i)\n";
	}
	//---------------------
	dt = param.val("Sz");
	if (dt > 0) { //Full magnetization Profile
		sx.open("Sx_profile.dat", mode);
		sy.open("Sy_profile.dat", mode);
		sz.open("Sz_profile.dat", mode);
		sxsysz.open("SxSySz_profile.dat", mode);

		correlation1.open("Correlation1.dat", mode);
		correlation2.open("Correlation2.dat", mode);

		sz_avrg.open("Sz_average_profile.dat", mode);
		sx_avrg.open("Sx_average_profile.dat", mode);

		sx.precision(15);
		sy.precision(15);
		sz.precision(15);
		sxsysz.precision(15);


		sz_avrg.precision(15);
		sx_avrg.precision(15);

		sx << "#Position=i-"
				<< "\t <Sx_i>\t"
				<< "\t <Sx_i Sx_{i+1}>\t"
				<< "\t <Sx_i Sx_{i+2}>\t"
				<< "\t <Sx_i Sx_{i+3}>\t"
				<< "\t <Sx_i Sx_{i+4}>\t"
				<< "\t <Sx_i Sz_{i+1}>\t"
				<< "\t\ttime\n";
		sy << "#Position=i-"
				<< "\t <Sy_i>\t"
				<< "\t <Sy_i Sy_{i+1}>\t"
				<< "\t <Sy_i Sy_{i+2}>\t"
				<< "\t <Sy_i Sy_{i+3}>\t"
				<< "\t <Sy_i Sy_{i+4}>\t"
				<< "\t <Sy_i Sz_{i+1}>\t"
				<< "\t\ttime\n";
		sz << "#Position=i-"
				<< "\t <Sz_i>\t"
				<< "\t <Sz_i Sz_{i+1}>\t"
				<< "\t <Sz_i Sz_{i+2}>\t"
				<< "\t <Sz_i Sz_{i+3}>\t"
				<< "\t <Sz_i Sz_{i+4}>\t"
				<< "\t (-1)^i<Sz_i>\t"
				<< "\t\ttime\n";

		sxsysz << "#Position=i-"
				<< "\t <Sx_i>\t"
				<< "\t <Sy_i>\t"
				<< "\t <Sz_i>\t"
				<< "\t\ttime\n";

		correlation1<< "#Position=i-"
				<< "\t <Sx_i Sx_{i+1}>\t"	//2
				<< "\t <Sx_i Sy_{i+1}>\t"	//3
				<< "\t <Sx_i Sz_{i+1}>\t"	//4
				<< "\t <Sy_i Sx_{i+1}>\t"	//5
				<< "\t <Sy_i Sy_{i+1}>\t"	//6
				<< "\t <Sy_i Sz_{i+1}>\t"	//7
				<< "\t <Sz_i Sx_{i+1}>\t"	//8
				<< "\t <Sz_i Sy_{i+1}>\t"	//9
				<< "\t <Sz_i Sz_{i+1}>\t"	//10
				<< "\t\ttime\n";

		correlation2<< "#Position=i-"
				<< "\t <Sx_i Sx_{i+2}>\t" 	//2
				<< "\t <Sx_i Sy_{i+2}>\t"	//3
				<< "\t <Sx_i Sz_{i+2}>\t"	//4
				<< "\t <Sy_i Sx_{i+2}>\t"	//5
				<< "\t <Sy_i Sy_{i+2}>\t"	//6
				<< "\t <Sy_i Sz_{i+2}>\t"	//7
				<< "\t <Sz_i Sx_{i+2}>\t"	//8
				<< "\t <Sz_i Sy_{i+2}>\t"	//9
				<< "\t <Sz_i Sz_{i+2}>\t"	//10
				<< "\t\ttime\n";

		sz_avrg << "#Position=i-"
				<< "\t 0.5*<Sz_i + next site>\t"
				<< "\t 0.5*<Sz_i Sz_{i+1} + next site>\t"
				<< "\t 0.5*<Sz_i Sz_{i+2} + next site>\t"
				<< "\t 0.5*<Sz_i Sz_{i+3} + next site>\t"
				<< "\t 0.5*<Sz_i Sz_{i+4} + next site>\t"
				<< "\t\ttime\n";

		sx_avrg << "#Position=i-"
				<< "\t 0.5*<Sx_i + next site>\t"
				<< "\t 0.5*<Sx_i Sx_{i+1} + next site>\t"
				<< "\t 0.5*<Sx_i Sx_{i+2} + next site>\t"
				<< "\t 0.5*<Sx_i Sx_{i+3} + next site>\t"
				<< "\t 0.5*<Sx_i Sx_{i+4} + next site>\t"
				<< "\t\ttime\n";
	}
	//---------------------
	dt = param.val("Q1Profile");
	if (dt > 0) { //Energy Profile and Current profile
		q1_profile.open("Q1_profile.dat", mode);
		q1_profile.precision(15);
		q1_profile << "#Position=i-" << "\t<Q1minus_i>\t" << dot
				<< "\t\ttime(or beta)\n";

	}
	//---------------------
	dt = param.val("Q2Profile");
	if (dt > 0) { //Full Q2 Profile
		q2_profile.open("Q2_profile.dat", mode);
		q2_profile.precision(15);
		q2_profile << "#Position=i-" << dot << setw(16) << "\t Entropy(i)"
				<< setw(16) << "\t Q2plus \t Q2minus \t time \t \n";
	}

	//exp(Ham) for the dynamics

	double tau = param.val("tau");
	const long int n_steps = param.val("T") / param.val("tau");
	TrotterExp expH(sites, param, -Cplx_i * tau);

// Time evolution
	double sz_total_initial = 9999;
	double sz_tot = 999; // will be used in part with magnetization and then output to cout;
	double sx_tot = 999; // will be used in part with magnetization and then output to cout;
	for (int n = 0; n <= n_steps ; ++n) {

		const double time = n * tau; //+param.val("time_shift");
		cout << "Time step #" << n << "/" << n_steps << "\ttime=" << time
				<< endl;

		vector<double> Myspec; //vector which will be the SVD spectrum

		//  ----- Compute observables ---
		//  ----- Entropy between sites i and i+1
		if (param.val("Entropy") != 0) {
			double entr = Entropy(psi, dot, Myspec, 1);
			ent << time << "\t" << setw(16) << setfill('0') << entr << "\t"
					<< BondDim(psi, dot) << "\t" << maxLinkDim(psi) << endl;

			if (param.val("SVD_spec") > 0) {
				if (n % int(param.val("SVD_spec") / tau) == 0) {
					spec << "\"t=" << time << "\"" << endl;
					int si = Myspec.size();
					for (int i = 0; i < si; i++) {
						spec << i + 1 << "\t" << Myspec[i] << "\t"
								<< sqrt(Myspec[i]) << "\t" << time << endl;
					}
					if (n < n_steps)
						spec << endl << endl;
				}
			}
		}


		// ------- Loschmidt echo -----
		if (param.val("Loschmidt") > 0)
			if (n % int(param.val("Loschmidt") / tau) == 0) {
				complex<double> echo = innerC(psi0,psi);
				loschmidt << time << "\t" << setw(16) << setfill('0')
						<< real(echo) << "\t" << imag(echo)
						<< endl;
			}

		// ------- entanglement Entropy Profile -----
		if (param.val("EntropyProfile") > 0)
			if (n % int(param.val("EntropyProfile") / tau) == 0) {
				entropy_profile << "\"t=" << time << "\"" << endl;
				for (int i = 1; i < N; i++) {
					double entr_std = Entropy(psi, i, Myspec, 1); // p log p
					entropy_profile << i + 0.5 - dot << "\t" << setw(16) << setfill('0')
							<< entr_std << "\t" << setw(16) << setfill('0')
							<< BondDim(psi, i) << "\t" << time << endl;
				}
				if (n < n_steps)
					entropy_profile << endl << endl;
			}
		// ------- Sz profile -------

		if (param.val("Sz") > 0) {
			if (n % int(param.val("Sz") / tau) == 0) {
				sx << "\"t=" << time << "\"" << endl;
				sy << "\"t=" << time << "\"" << endl;
				sz << "\"t=" << time << "\"" << endl;

				sxsysz << "\"t=" << time << "\"" << endl;

				correlation1 << "\"t=" << time << "\"" << endl;
				correlation2 << "\"t=" << time << "\"" << endl;

				sx_avrg << "\"t=" << time << "\"" << endl;
				sz_avrg << "\"t=" << time << "\"" << endl;

				double sz_left = 0, sz_right = 0, sz_dot = 0;
				sz_tot = 0;
				sx_tot = 0;

				double sx_odd = 0,  sz_odd = 0;
				double sxsx1_odd = 0, sxsx2_odd = 0, sxsx3_odd = 0, sxsx4_odd = 0;
				//double sysy1_odd = 0, sysy2_odd = 0, sysy3_odd = 0, sysy4_odd = 0;
				double szsz1_odd = 0, szsz2_odd = 0, szsz3_odd = 0, szsz4_odd = 0;

				double sxsz_odd  = 0; //sysz_odd  = 0;

				for (int i = 1; i <= N; i ++) {
					double s = Sz(psi, sites, i);
					double sx1 = Sx(psi, sites, i);
					double sy1 = Sy(psi, sites, i);
					double sxsx3 = 0, sxsx4 = 0;
					double sysy3 = 0, sysy4 = 0;
					double szsz3 = 0, szsz4 = 0;

					// some correlations have form C(i,i+2)
					// if i > N-2 I define correlation to be 0
					double sxsx1 = 0, sxsy1 = 0, sxsz1 = 0;
					double sysx1 = 0, sysy1 = 0, sysz1 = 0;
					double szsx1 = 0, szsy1 = 0, szsz1 = 0;
					double sxsx2 = 0, sxsy2 = 0, sxsz2 = 0;
					double sysx2 = 0, sysy2 = 0, sysz2 = 0;
					double szsx2 = 0, szsy2 = 0, szsz2 = 0;

					complex<double> kss = 0;

					if (i <= N - 1) {
						sxsx1 = real(Correlation(psi, sites, "Sx", "Sx", i, i+1));
						sxsy1 = real(Correlation(psi, sites, "Sx", "Sy", i, i+1));
						sxsz1 = real(Correlation(psi, sites, "Sx", "Sz", i, i+1));

						sysx1 = real(Correlation(psi, sites, "Sy", "Sx", i, i+1));
						sysy1 = real(Correlation(psi, sites, "Sy", "Sy", i, i+1));
						sysz1 = real(Correlation(psi, sites, "Sy", "Sz", i, i+1));

						szsx1 = real(Correlation(psi, sites, "Sz", "Sx", i, i+1));
						szsy1 = real(Correlation(psi, sites, "Sz", "Sy", i, i+1));
						szsz1 = real(Correlation(psi, sites, "Sz", "Sz", i, i+1));

					}
					if (i <= N - 2) {
						sxsx2 = real(Correlation(psi, sites, "Sx", "Sx", i, i+2));
						sxsy2 = real(Correlation(psi, sites, "Sx", "Sy", i, i+2));
						sxsz2 = real(Correlation(psi, sites, "Sx", "Sz", i, i+2));

						sysx2 = real(Correlation(psi, sites, "Sy", "Sx", i, i+2));
						sysy2 = real(Correlation(psi, sites, "Sy", "Sy", i, i+2));
						sysz2 = real(Correlation(psi, sites, "Sy", "Sz", i, i+2));

						szsx2 = real(Correlation(psi, sites, "Sz", "Sx", i, i+2));
						szsy2 = real(Correlation(psi, sites, "Sz", "Sy", i, i+2));
						szsz2 = real(Correlation(psi, sites, "Sz", "Sz", i, i+2));

						kss = KSS(psi, sites, i);
					}
					if (i <= N - 3) {
						sxsx3 = real(Correlation(psi, sites, "Sx", "Sx", i, i+3));
						sysy3 = real(Correlation(psi, sites, "Sy", "Sy", i, i+3));
						szsz3 = real(Correlation(psi, sites, "Sz", "Sz", i, i+3));
					}
					if (i <= N - 4) {
						sxsx4 = real(Correlation(psi, sites, "Sx", "Sx", i, i+4));
						sysy4 = real(Correlation(psi, sites, "Sy", "Sy", i, i+4));
						szsz4 = real(Correlation(psi, sites, "Sz", "Sz", i, i+4));
					}

					sz_tot += s;
					sx_tot += sx1;
					if (i < dot)
						sz_left += s;
					if (i > dot)
						sz_right += s;
					if (i == dot)
						sz_dot += s;
					sz << i - dot  << "\t"			//column ID below
						<< s << "\t"				//2
						<< szsz1 << "\t"			//3
						<< szsz2 << "\t"			//4
						<< szsz3 << "\t"			//5
						<< szsz4 << "\t"			//6
						<< pow(-1, i ) * s << "\t"	//7
						 << "\t"<< time << endl;
					sx << i - dot  << "\t"			//column ID below
						<< sx1 << "\t"				//2
						<< sxsx1 << "\t"			//3
						<< sxsx2 << "\t"			//4
						<< sxsx3 << "\t"			//5
						<< sxsx4 << "\t"			//6
						<< sxsz1 << "\t"				//7
						<< time << endl;
					sy << i - dot  << "\t"			//column ID below
						<< sy1 << "\t"				//2
						<< sysy1 << "\t"			//3
						<< sysy2 << "\t"			//4
						<< sysy3 << "\t"			//5
						<< sysy4 << "\t"			//6
						<< sysz1 << "\t"			//7
						<< time << endl;
					sxsysz << i - dot  << "\t"		//column ID below
						<< sx1 << "\t"				//2
						<< sy1 << "\t"				//3
						<< s << "\t"				//4
						<< real(kss) << "\t"		//5
						<< time << endl;

					correlation1 << i - dot  << "\t"		//column ID below
						<< sxsx1 << "\t"				//2
						<< sxsy1 << "\t"				//3
						<< sxsz1 << "\t"				//4
						<< sysx1 << "\t"				//5
						<< sysy1 << "\t"				//6
						<< sysz1 << "\t"				//7
						<< szsx1 << "\t"				//8
						<< szsy1 << "\t"				//9
						<< szsz1 << "\t"				//10
						<< time << endl;

					correlation2 << i - dot  << "\t"		//column ID below
						<< sxsx2 << "\t"				//2
						<< sxsy2 << "\t"				//3
						<< sxsz2 << "\t"				//4
						<< sysx2 << "\t"				//5
						<< sysy2 << "\t"				//6
						<< sysz2 << "\t"				//7
						<< szsx2 << "\t"				//8
						<< szsy2 << "\t"				//9
						<< szsz2 << "\t"				//10
						<< time << endl;

					if ( i % 2 == 1) { //odd site
						sz_odd = s;
						sx_odd = sx1;
						szsz1_odd = szsz1;
						szsz2_odd = szsz2;
						szsz3_odd = szsz3;
						szsz4_odd = szsz4;

						sxsx1_odd = sxsx1;
						sxsx2_odd = sxsx2;
						sxsx3_odd = sxsx3;
						sxsx4_odd = sxsx4;
						sxsz_odd = sxsz1;
					} else {
						sz_avrg << i - dot + 1 << "\t"
								<< 0.5 * (s + sz_odd) << "\t"
								<< 0.5 * (szsz1 + szsz1_odd) << "\t"
								<< 0.5 * (szsz2 + szsz2_odd) << "\t"
								<< 0.5 * (szsz3 + szsz3_odd) << "\t"
								<< 0.5 * (szsz4 + szsz4_odd) << "\t"
								<< time << endl;
						sx_avrg << i - dot + 1 << "\t"
								<< 0.5 * (sx1 + sx_odd) << "\t"
								<< 0.5 * (sxsx1 + sxsx1_odd) << "\t"
								<< 0.5 * (sxsx2 + sxsx2_odd) << "\t"
								<< 0.5 * (sxsx3 + sxsx3_odd) << "\t"
								<< 0.5 * (sxsx4 + sxsx4_odd) << "\t"
								<< 0.5 * (sxsz1  + sxsz_odd) << "\t"
								<< time << endl;
					}
				}

				{ //I need this part to separate time steps in *.dat files (for gnuplot)
					sx << "\n\n";
					sx_avrg << "\n\n";

					sy << "\n\n";

					sz << "\n\n";
					sz_avrg << "\n\n";

					sxsysz << "\n\n";

					correlation1 << "\n\n";
					correlation2 << "\n\n";
				}
				//sz_tot += Sz(psi, sites, N-1) + Sz(psi, sites, N);
				if (n == 0) {
					sz_total_initial = sz_tot;
				}
				cout << "\n<Sz_left>=" << sz_left << "\t" << "<Sz_right>="
						<< sz_right << "\t" << "<Sz_DOT>=" << sz_dot << "\t"
						<< "<Sz_tot>=" << sz_tot
						<<endl;
			}
		}
		// ------- Energy Profile -------
		if (param.val("Q1Profile") > 0 ) {
			if (n % int(param.val("Q1Profile") / tau) == 0) {
				q1_profile << "\"t=" << time << "\"" << endl;
				for (int i = 1; i <= N - 3; i++) {
					const complex<double> q1 = Q1(psi, sites, i);
					const double en = real(q1);
					const double q1minus = imag(q1);
					q1_profile << i - dot
							<< "\t" << en
							<< "\t" << q1minus
							<< "\t" << time << endl;
				}
				q1_profile << "\n\n"; //I need this part to separate time steps in *.dat files (for gnuplot)
			}
		}
		// ------- Q2 Profile -------
		if (param.val("Q2Profile") > 0 ) {
			if (n % int(param.val("Q2Profile") / tau) == 0) {
				q2_profile << "\"t=" << time << "\"" << endl;
				for (int i = 1; i <= N - 5; i++) {
					const complex<double> q2 = Q2(psi, sites, i);
					const double q2_plus =  real(q2);
					const double q2_minus = imag(q2);
					q2_profile << i - dot << "\t" << q2_plus << "\t"
							<< q2_minus << "\t" << time << endl;
				}
				q2_profile << "\n\n"; //I need this part to separate time steps in *.dat files (for gnuplot)
			}
		}

		if (n < n_steps) {
			//MPS psi_temp = psi;
			cout << "Time evol" << endl;
			expH.Evolve(psi, args);
			psi.orthogonalize(args);

			// Hadamar gates act at time == 5
			if(  (int)param.val("Measurement") != 0
					&& n > 0
					&& n  % (int)(param.val("Measurement")/param.val("tau")) == 0
					){
				//AlphaGate(N/2-1, -0.5);
				cout << "TIME IS == \"Measurement\". I ACT WITH UpToDown GATES" << endl;
				SigmaXGate(N/2-1);
//				HadamarGate(N/2-1);
//				HadamarGate(N/2  );
//				HadamarGate(N/2+1);
//				HadamarGate(N/2+2);

				psi.noPrime();
			}

			double energy = real(innerC(psi, H, psi));
			cout << "max bond dim = " << maxLinkDim(psi) << endl;
			//cout << "Norm = " << real(innerC(psi, psi)) << endl;
			cout << "Energy = " << energy << endl;
			cout << "E(t) - E(0) = " << (energy - energy_initial) << endl;
			cout << "Sz(t) - Sz(0) = " << (sz_tot - sz_total_initial) << endl;
		}
	}
	cout << "\nTime evolution complete.\n";
	println("Done !");
	return 0;
}

