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

#include <cstdlib>
#include "observables_GS.h"

using namespace itensor;
using namespace std;
using namespace std::chrono;

//------------------------------------------------------------------
//The function below translate numbers (etc.) into charcter strings
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
		operator[]("tau") = 0.1;  //time step for the unitary evolution
		operator[]("T") = 2;  //Total (final) time
		operator[]("Sz") = 0;
		operator[]("SVD_spec") = 0; //SVD spectrum
		operator[]("max_bond") = 4000;  //maximum bond dimension
		operator[]("trunc") = 1e-10;  //maximum truncation error
		operator[]("energy") = 1e-10;  //convergence criterium on the energy
		operator[]("sweeps") = 999;  //maximum number of sweeps in the DMRG
		operator[]("TrotterOrder") = 2;
		operator[]("GroundState") = 0;
		operator[]("DomainWall") = 0;
		operator[]("RandomState") = 0;
		operator[]("Neel") = 0;
		operator[]("begin") = 1;
		//operator[]("end") = 10;
		operator[]("hL") = 0; //chemical potential
		operator[]("hR") = 0;
		operator[]("Q2Profile") = 0;
		operator[]("CurrentProfile") = 0;
		operator[]("Current") = 0;
		operator[]("Entropy") = 0; //entanglement entropy p*log*p between left and right parts of system
		operator[]("EntropyProfile") = 0; // Entropy Profile- parameter 0 -> nothing, dt>0 each second=integer parameter
		operator[]("EnergyProfile") = 0;
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
		dot = (N) / 2;  //Position of the "dot"
		//cout << "The dot is on site #" << dot << endl;
		//if ((2*N)<=3) cout<<"Error, N="<<N<<" is too small.\n",exit(0);

		for (int j = 1; j < N-1; ++j) {
			//Strange coefficients are needed to match with
			// spin Pauli matrices instead of Sx Sy
			ampo += J * 4 * 0.25, "S+", j, "S-", j + 2; // 0.5 (SpSm+ SmSp) = SxSx + SySy
			ampo += J * 4 * 0.25, "S-", j, "S+", j + 2;
			ampo += J * -8 * 0.25, "S+", j, "Sz", j + 1, "S-", j + 2;
			ampo += J * -8 * 0.25, "S-", j, "Sz", j + 1, "S+", j + 2;
			//cout << "j = "<< j << "/ " << N-3 << endl;
			cout << "site (" << j << ", " << j + 1 << ", " << j + 2 << ")"
					<< endl;
			if(j <= dot){
				mu = hL;
			}else{
				mu = hR;
			}
			ampo += mu, "Sz", j;
		}
		ampo += mu, "Sz", N-1;
		ampo += mu, "Sz", N;
		cout << "H = 3site is construcnted" << endl;
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
			TimeGates(begin4, end, tau, sites, param); 	 //C
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
		cout << "Gates starts from " << begin << endl;
		for (int j = begin; j < end - 1; j += step) {
			cout << "j = (" << j << ", " << j + 1 << ", " << j + 2 << ")"
					<< endl;
			//this part act on real sites
			auto hh = J * 4 * 0.25 * op(sites, "Sp", j) * op(sites, "Id", j + 1)
					* op(sites, "Sm", j + 2);
			hh += J * 4 * 0.25 * op(sites, "Sm", j) * op(sites, "Id", j + 1)
					* op(sites, "Sp", j + 2);
			hh += -J * 8 * 0.25 * op(sites, "Sp", j) * op(sites, "Sz", j + 1)
					* op(sites, "Sm", j + 2);
			hh += -J * 8 * 0.25 * op(sites, "Sm", j) * op(sites, "Sz", j + 1)
					* op(sites, "Sp", j + 2);
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

// Making an initial state

	ThreeSiteHamiltonian Init_H(sites, param);
	auto H0 = toMPO(Init_H.ampo);
	double energy_initial = 0;


	if (param.longval("GroundState") == 1) {
		// GS of the initial Hamiltonian
		cout << "initial state is GS" << endl;
		auto sweeps = Sweeps(999); //number of sweeps is 5
		sweeps.maxdim() = 10, 20, 50, 50, 100, 300, 4000;
		sweeps.cutoff() = 1E-10;

		psi = randomMPS(sites);

		cout << "Before DMRG" << endl;

		cout << "Norm is = " << inner(psi, psi) << endl;
		energy_initial = inner(psi, H0, psi); //<psi|H0|psi>
		cout << "1. Initial energy=" << energy_initial << endl;

		MyDMRGObserver obs(psi, param.val("energy"));

		tie(energy_initial, psi) = dmrg(H0, psi, sweeps, obs, "Quiet");

		cout << "After DMRG" << endl;

		//psi0 = psi; //we create to states. Psi for my time evolution, psi0 for standard one

	} else if (param.longval("DomainWall") == 1) {
		cout << "initil state is ++++----" << endl;
		auto initState = InitState(sites);
		for (int i = 1; i <= N / 2; ++i)
			initState.set(i, "Up");
		for (int i = N / 2 + 1; i <= N; ++i)
			initState.set(i, "Dn");
		psi = MPS(initState);

	} else if (param.longval("Neel") == 1) {
		// Initial state: |+- +- +- >
		cout << "initil state is {++++} {--} {Neel}" << endl;
		auto initState = InitState(sites);
		for (int i = 1; i <= N / 2; ++i){
			initState.set(i, "Up");
		}
		initState.set(N/2 + 1, "Dn");
		initState.set(N/2 + 2, "Dn");

		bool flag_up_spin = true;
		for (int i = N / 2 + 3; i <= N; ++i){
			if (flag_up_spin){
				initState.set(i, "Up");
				flag_up_spin = false;
			} else{
				initState.set(i, "Dn");
				flag_up_spin = true;
			}
		}
		psi = MPS(initState);

	} else if (param.longval("RandomState") == 1) {
		cout << "initil state is {RND and ---}" << endl;
		auto initState = InitState(sites);
		std::srand(std::time(0));
		for (int i = 1; i <= N / 2 ; ){

			int rnd = rand();
			if (rnd % 2 == 0){
				initState.set(i, "Up");
				i+=1;
			}else{
				initState.set(i, "Up");
				i+=1;
				initState.set(i, "Dn");
				i+=1;
			}
		}
		for (int i = N / 2 + 1; i <= N; ++i){
			initState.set(i, "Dn");
		}
		psi = MPS(initState);


//		// Random state
//		cout << "initial state  is random state" << endl;
//		psi = randomMPS(sites);
//		psi.normalize();
	} else {
		cout << "Choose: GroundState 1 or Neel 1 or DomainWall 1 or RandomState 1" << endl;
		return 1;
	}
	//cout << "PSI = " << psi << endl;

	energy_initial = inner(psi, H0, psi); //<psi|H0|psi>
	cout << "2. Initial energy psi =" << energy_initial << endl;

//--------------------------------------------------------------

	//Hamiltonian for the dynamics
	ThreeSiteHamiltonian Ham(sites, param);
	const int dot = Ham.dot;
	auto H = toMPO(Ham.ampo);
	cout << "H constructed" << endl;

// Output and observables
// _______
	ofstream ent, spec, entropy_profile, sz, sz_avrg, energy_beta, energy_profile,
				q1minus_profile, q2_profile; //here I'm defining output streams == files
		ios_base::openmode mode;
		mode = std::ofstream::out; //Erase previous file (if present)

	double dt = param.val("Entropy");
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
		sz.open("Sz_profile.dat", mode);
		sz.precision(15);
		sz << "#Position=i-" << "\t<Sz_i>\t" << "\t(-1)^i<Sz_i>\t"
				<< "\t\ttime\n";

		sz_avrg.open("Sz_average_profile.dat", mode);
		sz_avrg.precision(15);
		sz_avrg << "#Position=i-" << "\t0.5<Sz_2+1i> + 0.5<Sz_2+1i>\t"
				<< "\t\ttime\n";

	}
	//---------------------
	dt = param.val("EnergyProfile");
	if (dt > 0) { //Energy Profile
		energy_profile.open("Energy_profile.dat", mode);
		energy_profile.precision(15);
		energy_profile << "#Position=i-" << "\t<Ham_i>\t" << dot
				<< "\t\ttime(or beta)\n";

		//Q1minus Profile is initialized simultaniously with energy Profile
		q1minus_profile.open("Q1minus_profile.dat", mode);
		q1minus_profile.precision(15);
		q1minus_profile << "#Position=i-" << "\t<Q1minus_i>\t" << dot
				<< "\t\ttime(or beta)\n";

	}
	//---------------------
	dt = param.val("Q2Profile");
	if (dt > 0) { //Full entropy Profileile
		q2_profile.open("Q2_profile.dat", mode);
		q2_profile.precision(15);
		q2_profile << "#Position=i-" << dot << setw(16) << "\t Entropy(i)"
				<< setw(16) << "\t Q2plus \t Q2minus \t time \t \n";
	}


	//exp(Ham) for the dynamics
	auto args = Args("Method=", "DensityMatrix", "Cutoff", param.val("trunc"),
			"MaxDim", param.longval("max_bond"), "Normalize", true); // for FitApplyMPO RENAME

	double tau = param.val("tau");
	const long int n_steps = param.val("T") / param.val("tau");
	TrotterExp expH(sites, param, -Cplx_i * tau);

// Time evolution
	bool go_on = true;
	for (int n = 0; n <= n_steps && go_on; ++n) {
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
				sz << "\"t=" << time << "\"" << endl;
				sz_avrg << "\"t=" << time << "\"" << endl;
				double sz_tot = 0, sz_left = 0, sz_right = 0, sz_dot = 0;
				double sz_odd = 0;
				for (int i = 1; i <= N; i ++) {
					double s = Sz(psi, sites, i);
					sz_tot += s;
					if (i < dot)
						sz_left += s;
					if (i > dot)
						sz_right += s;
					if (i == dot)
						sz_dot += s;
					sz << i - dot  << "\t" << s << "\t"
							<< pow(-1, i ) * s << "\t" << time << endl;

					if ( i % 2 == 1) { //odd site
						sz_odd = s;
					} else {
						sz_avrg << i - dot + 1 << "\t"
								<< 0.5 * (s + sz_odd) << "\t" << time << endl;
					}
				}

				{ //I need this part to separate time steps in *.dat files (for gnuplot)
					sz << "\n\n";
					sz_avrg << "\n\n";
				}
				cout << "\n<Sz_left>=" << sz_left << "\t" << "<Sz_right>="
						<< sz_right << "\t" << "<Sz_DOT>=" << sz_dot << "\t"
						<< "<Sz_tot>=" << sz_tot << endl;
			}
		}
		// ------- Energy Profile -------
		// ------- Energy Profile -------
		if (param.val("EnergyProfile") > 0 ) {
			if (n % int(param.val("EnergyProfile") / tau) == 0) {
				energy_profile << "\"t=" << time << "\"" << endl;
				q1minus_profile << "\"t=" << time << "\"" << endl;
				for (int i = 1; i <= N - 3; i++) {
					const complex<double> q1 = Q1(psi, sites, i);
					const double en = real(q1);
					energy_profile << i - dot  << "\t" << en << "\t"
							<< time << endl;
					const double q1minus = imag(q1);
					q1minus_profile << i - dot << "\t" << q1minus
							<< "\t" << time << endl;
				}
				energy_profile << "\n\n"; //I need this part to separate time steps in *.dat files (for gnuplot)
				q1minus_profile << "\n\n"; //I need this part to separate time steps in *.dat files (for gnuplot)
			}
		}
		// ------- Q2 Profile -------
		if (param.val("Q2Profile") > 0 ) {
			if (n % int(param.val("EnergyProfile") / tau) == 0) {
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

			double energy = real(innerC(psi, H, psi));
			cout << "max bond dim = " << maxLinkDim(psi) << endl;
			cout << "Norm = " << real(innerC(psi, psi)) << endl;
			cout << "Energy = " << energy << endl;
			//cout << "dE/E = " << (energy - energy_initial)/energy_initial << endl;

		}
	}
	cout << "\nTime evolution complete.\n";
	println("Done !");
	return 0;
}

