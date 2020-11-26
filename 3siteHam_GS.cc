// C++ code for time averaged density matrix. We use operator-state duality |psi><psi| -> |psi>|psi>
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
#include "observables_GS.h"


using namespace itensor;
using namespace std;
using namespace std::chrono;


// _USE_QN_ may be defined in the compilation command (see Makefile)
// in order to obtain the "IQ" version of the program.

#ifdef _USE_QN_
typedef IQMPS MPS_;
typedef IQMPO MPO_;
typedef IQTensor ITensor_;
#else
typedef MPS MPS_;
typedef MPO MPO_;
typedef ITensor ITensor_;
#endif

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
				cout << "Error: Parameter " << var_name << " is not defined.\n", exit(0);
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
						cerr << "Error: missing value after " << var_name << endl, exit(0);
					operator[](var_name) = char2double(argv[n]);
				} else {
					cerr << "Syntax error :" << var_name << endl;
					cout << "List of command-line parameters :";
					PRint(cout);
					exit(0);
				}
			}
		}
};
//_____________________________________________________

// class of PUBLIC parameters
class TEDM_param: public Parameters {
	public:
		TEDM_param() { //Constructor
			//Specify below all the allowed parameter names,
			//and their default values
			operator[]("N") = 10; //Length of the chain
			operator[]("J") = 1.0;
			operator[]("mu") = 1.0;
			operator[]("tau") = 0.1;  //time step for the unitary evolution
			operator[]("T") = 2;  //Total (final) time
			operator[]("Entropy") = 0; //entanglement entropy p*log*p between left and right parts of system
			operator[]("Eprof") = 0; // Entropy profile - parameter 0 -> nothing, dt>0 each second=integer parameter
			operator[]("EnergyProf") = 0;
			operator[]("Sz") = 0;
			operator[]("Entr_states") = 0;
			operator[]("SVD_spec") = 0; //SVD spectrum
			operator[]("H_spec") = 0; //Ham spectrum
			operator[]("max_bond") = 4000;  //maximum bond dimension
			operator[]("trunc") = 1e-10;  //maximum truncation error
			operator[]("energy") = 1e-10;  //convergence criterium on the energy
			operator[]("sweeps") = 999;  //maximum number of sweeps in the DMRG
			operator[]("TrotterOrder") = 4;
			operator[]("GS") = 1;
			operator[]("antal") = 0;
			operator[]("rnd_state") = 0;			
			operator[]("XXZ") = 0;
			operator[]("PBC") = 0;
			operator[]("beta") = 0;
		}
};
//_____________________________________________________

//I'm creating a 3 site Hamiltonian for system of N sites



class TEDM_Hamiltonian {
	public:
		int dot;
	private:
		int N;
		void init(const TEDM_param &param) {    //.init (param)
			const double J = param.val("J");
			const double mu = param.val("mu");
			dot = (N) / 2;  //Position of the "dot"
			//cout << "The dot is on site #" << dot << endl;
			//if ((2*N)<=3) cout<<"Error, N="<<N<<" is too small.\n",exit(0);
			if ( param.val("XXZ") ){
				for (int j = 1; j  < N; ++j) {
					ampo +=  J * 0.5 , "S+", j, "S-", j + 1;
					ampo +=  J * 0.5 , "S-", j, "S+", j + 1;
					ampo +=  J * 1.0 , "Sz", j, "Sz", j + 1;
				}
				cout << "H = XXZ is construcnted"<< endl;
				if (param.val("PBC")){
					ampo +=  J * 0.5 , "S+", N, "S-", 1;
					ampo +=  J * 0.5 , "S-", N, "S+", 1;
					ampo +=  J * 1.0 , "Sz", N, "Sz", 1;
					cout << "XXZ is periodic"<< endl;
				}

			} else{

				for (int j = 2; j  < N; ++j ) {
					//Strange coefficients are needed to match with 
					// spin Pauli matrices instead of Sx Sy 
					ampo += J *  4 * 0.25 , "S+", j-1, "S-", j + 1; // 0.5 (SpSm+ SmSp) = SxSx + SySy
					ampo += J *  4 * 0.25 , "S-", j-1, "S+", j + 1;
					ampo += J * -8 * 0.25 , "S+", j-1, "Sz", j, "S-", j + 1;
					ampo += J * -8 * 0.25 , "S-", j-1, "Sz", j ,"S+", j + 1;
					//cout << "j = "<< j << "/ " << N-3 << endl;
					cout << "site (" << j-1 << " " << j+1 <<")" << endl;
					ampo += mu, "Sz", j ;  
				}
					ampo += mu, "Sz", 1 ;
					ampo += mu, "Sz", N ;
				cout << "H = 3site is construcnted"<< endl;
				if (param.val("PBC")){
					// This part realizes Periodic Boundary Condition (PBC)
					// term (N-1,N,1)
					ampo += J *  4 * 0.25 , "S+", N-1, "S-", 1;
					ampo += J *  4 * 0.25 , "S-", N-1, "S+", 1;
					ampo += J * -8 * 0.25 , "S+", N-1, "Sz", N, "S-", 1;
					ampo += J * -8 * 0.25 , "S-", N-1, "Sz", N, "S+", 1;
					cout << "PBC; sites (" << N-1 << " " << 1 <<"), ";


					// term (N,1,2)
					ampo += J *  4 * 0.25 , "S+", N, "S-", 2;
					ampo += J *  4 * 0.25 , "S-", N, "S+", 2;
					ampo += J * -8 * 0.25 , "S+", N, "Sz", 1, "S-", 2;
					ampo += J * -8 * 0.25 , "S-", N, "Sz", 1, "S+", 2;
					cout << "(" << N << " " << 2 <<")" << endl;

					cout << "H = 3site is periodic" << endl;
				}

			}

		}
	public:
		AutoMPO ampo;  //variable ampo 
		//Define a constructor for this class 'TEDM_Hamiltonian'
		TEDM_Hamiltonian(const SiteSet &sites, const TEDM_param &param) :
			ampo(sites) {
				N = length(sites); // size of Hamiltonian comes with onject SiteSet sites, not just as number N
				init(param);   // initializing the Hamiltonian
				cout << "A Hamiltonian with " << N << " sites was constructed." << endl;
			}

};


// OBSERVABLES



// Time evolution e^-iHt = 1-iHt
template <class T>
void TimeEvolution (MPS& psi, MPO& Ham, double & tau0, const int& o, T& args){
	complex<double> tau = tau0;
	vector <complex <double>> time_steps;
	Cplx t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0, t7 = 0;
	if (o == 2) {
		//Approx. with error O(tau^3)
		t1 = 0.5 * ( 1 + Cplx_i) * tau;
		t2 = 0.5 * (-1 + Cplx_i) * tau;
		time_steps.push_back(t1);
		time_steps.push_back(t2);
	}
	if (o == 3) {
		//Approx. with error O(tau^4) [thanks to Kemali]	//Tested -it's ok
		t1 = .10566243270259355887 - .39433756729740644113 * Cplx_i;
		t2 = Cplx_i * t1;
		t3 = conj(t2);
		t4 = Cplx_i * t3;
		t1 *= Cplx_i * tau;
		t2 *= Cplx_i * tau;
		t3 *= Cplx_i * tau;
		t4 *= Cplx_i * tau;
		time_steps.push_back(t1);
		time_steps.push_back(t2);
		time_steps.push_back(t3);
		time_steps.push_back(t4);
	}
	if (o == 4) {
		//Approx. with error O(tau^5) [thanks to Kemal twice]		//Tested -it's ok
		t1 =  0.25885339861091821723 + 0.04475613401114190287 * Cplx_i;
		t2 = -0.03154685814880379274 + 0.24911905427556321757 * Cplx_i;
		t3 =  0.19082905211066719664 - 0.23185374923210605447 * Cplx_i;
		t4 =  0.1637288148544367438753;
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
		time_steps.push_back(t1);
		time_steps.push_back(t2);
		time_steps.push_back(t3);
		time_steps.push_back(t4);
		time_steps.push_back(t5);
		time_steps.push_back(t6);
		time_steps.push_back(t7);
	}

	//MPS Hpsi;
	//cout << "Function TimeEvolution" << endl; 
	//double norm = real(innerC(psi,psi));
	//cout << "BEFORE: Energy = " << norm << endl;
	MPS Hpsi;
	for (auto t : time_steps){
		Hpsi = applyMPO(Ham, psi, args);
		Hpsi.noPrime("Site");
		psi = sum (psi , t * Hpsi, args);
		psi.noPrime("Site");
	}

	//norm = real(innerC(psi,psi));
	//cout << "After: Energy = " << norm << endl;
	//psi.normalize();
}

//------------------------------------------------------------------
class MyDMRGObserver:public DMRGObserver { // vector<int> ~ DMRGObserver<ITensor>, where ITensor are tensor object type.
	double previous_energy;
	const double precision;
	public:
	MyDMRGObserver(const MPS & psi,double prec=1e-10)
		: DMRGObserver(psi)
		  , precision(prec) {
		  }
	bool checkDone(const Args& args = Args::global()) {
		const double energy = args.getReal("Energy",0);
		cout<<"    Energy change:"<<energy-previous_energy<<endl;
		if (abs(energy-previous_energy)<precision) {
			cout<<"   Energy has converged -> stop the DMRG.\n";
			return true;
		}
		else {
			previous_energy=energy;
			return false;
		}
	}
};


// main
int main(int argc, char *argv[]) {

	TEDM_param param;
	param.ReadArguments(argc, argv); //Now param contains the parameters, default values or those provided on the command-line

	param.PRint(cout);// Print parameters
	cout.precision(15);

	const int N = param.longval("N");
#ifdef _USE_QN_
	cout<<"IQ version.\n";
#else
	cout << "nonIQ version.\n";
#endif

#define HILBERT_SPACE SpinHalf

	HILBERT_SPACE sites(N, { "ConserveQNs=", false }); //right above
	//auto sites = CustomSpin(N,{"2S=",3,"ConserveQNs=", false});
	MPS_ psi, psi0;
	double J = param.val("J");
	param["J"] = J;

	//--------------------------------------------------------------

	// Making an initial state
	/*
	   sites = HILBERT_SPACE(N, { "ConserveQNs=", false });
	   param["JR"] = 0;
	   TEDM_Hamiltonian Init_H(sites, param);
	   auto H0 = toMPO(Init_H.ampo);

	   param["JR"] = -1;

	   auto initState = InitState(sites, "Up");
	//Initial state is a polarized state: |++++>
	//auto initState = InitState(sites);
	//auto initState = randomMPS(sites); 
	//for (int i = 1; i <= N; ++i) {
	//	//state.set(i,"Up");
	//	if (i % 2 == 1)
	//		initState.set(i, "Up");
	//	else
	//		initState.set(i, "Dn");
	//}

	//for(int i = 1; i <= N; ++i) initState.set(i,"Up");
	psi = MPS_(initState);
	//psi.position(1);
	//psi.normalize();
	psi0 = MPS_(initState);
	complex<Real> energy = innerC(psi, H0, psi); //<psi|H0|psi>
	cout << "Initial energy=" << energy << endl;
	 */
	//--------------------------------------------------------------

	TEDM_Hamiltonian Init_H(sites,param);
	auto H0 = toMPO(Init_H.ampo); 
	double energy = 0;

	if (param.longval("GS") == 1) {
		// GS of the initial Hamiltonian
		cout << "initial state is GS" << endl;
		auto sweeps = Sweeps(999); //number of sweeps is 5
		sweeps.maxdim() = 10,20,50,50,100,300,4000;
		sweeps.cutoff() = 1E-10;

		psi0 = randomMPS(sites);

		cout << "Before DMRG" << endl;
		psi = psi0;
		cout << "Norm is = " << inner(psi,psi) << endl;
		energy = inner(psi,H0,psi); //<psi|H0|psi>
		cout<<"1. Initial energy="<<energy<<endl;

		MyDMRGObserver obs(psi0, param.val("energy"));
		auto [en0, psi_t] = dmrg(H0,psi0,sweeps,obs,"Quiet");
		psi = psi_t;
		energy = en0;
		cout << "After DMRG" << endl;

		cout << "Norm is = " << inner(psi,psi) << endl;
		energy = inner(psi,H0,psi); //<psi|H0|psi>
		cout<<"1. Initial energy="<<energy<<endl;

		psi0 = psi; //we create to states. Psi for my teme evolution, psi0 for standard one

	} else if (param.longval("antal") == 1) {
		// Initial state: |+- +- +- >
		cout << "initil state is +-+-+-+-" << endl;
		auto initState = InitState(sites);
		for(int i = 1; i <= N/2; ++i) initState.set(i,"Up");
		for(int i = N/2+1; i <= N; ++i) initState.set(i,"Dn");
		psi = MPS_(initState);
		energy = inner(psi,H0,psi); //<psi|H0|psi>
		cout<<"Initial energy="<<energy<<endl;

		psi0 = psi;
	} else if (param.longval("rnd_state") == 1){
		// Random state
		cout << "initial state  is random state" << endl;
		psi0 = randomMPS(sites);

		psi = psi0;		
	}
	cout << "PSI = " << psi << endl;

	energy = inner(psi,H0,psi); //<psi|H0|psi>
	cout<<"2. Initial energy psi ="<<energy<<endl;
	energy = inner(psi0,H0,psi0); //<psi0|H0|psi0>
	cout<<"2. Initial energy psi0 ="<<energy<<endl;
	//--------------------------------------------------------------

	//Hamiltonian for the dyamics

	//this part will be used to quench
	//param["h"]=h0;
	//TEDM_Hamiltonian Ham1(sites,param);
	//auto H1 = MPO_(Ham1.ampo);

	param["J"] = 1;
	TEDM_Hamiltonian Ham(sites, param);
	const int dot = Ham.dot;
	auto H = toMPO(Ham.ampo);
	cout << "H constructed.\n";
	cout.flush();

	//exponentiation of H
	double tau = param.val("tau");
	const int o = param.val("TrotterOrder");
	MPO expH1, expH2, expH3, expH4, expH5, expH6, expH7;

	auto args = Args("Method=", "DensityMatrix", "Cutoff", param.val("trunc"), "MaxDim",
			param.longval("max_bond"), "Normalize", true);	// for FitApplyMPO
	const long int n_steps = param.val("T") / param.val("tau");

	/*
	if (n_steps > 0) {
		Cplx t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0, t7 = 0;
		bool ok = false;
		if (o == 1) {
			cout << "constructing toExpH" << endl;
			//Approx. with error O(tau^3)
			t1 = Cplx_i * tau;
			t2 = 0;
			t3 = 0;
			t4 = 0;
			ok = true;
			expH1 = toExpH(Ham.ampo, t1);
		}
		if (o == 2) {
			//Approx. with error O(tau^3)
			cout << "constructing toExpH" << endl;
			t1 = 0.5 * (1 + Cplx_i) * tau;
			t2 = 0.5 * (-1 + Cplx_i) * tau;
			t3 = 0;
			t4 = 0;
			ok = true;
			expH1 = toExpH(Ham.ampo, t1);
			expH2 = toExpH(Ham.ampo, t2);
		}
		if (o == 3) {
			//Approx. with error O(tau^4) [thanks to Kemal]
			//WARNINIG more less TESTED
			//Tested -it's ok
			t1 = .10566243270259355887 - .39433756729740644113 * Cplx_i;
			t2 = Cplx_i * t1;
			t3 = conj(t2);
			t4 = Cplx_i * t3;
			t1 *= Cplx_i * tau;
			t2 *= Cplx_i * tau;
			t3 *= Cplx_i * tau;
			t4 *= Cplx_i * tau;
			ok = true;
			cout << "constructing toExpH" << endl;
			expH1 = toExpH(Ham.ampo, t1);
			expH2 = toExpH(Ham.ampo, t2);
			expH3 = toExpH(Ham.ampo, t3);
			expH4 = toExpH(Ham.ampo, t4);
		}
		if (o == 4) {
			//Approx. with error O(tau^5) [thanks to Kemal twice]
			//WARNINIG Not fully TESTED
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
			ok = true;
			cout << "constructing toExpH" << endl;
			expH1 = toExpH(Ham.ampo, t1 );
			expH2 = toExpH(Ham.ampo, t2 );
			expH3 = toExpH(Ham.ampo, t3 );
			expH4 = toExpH(Ham.ampo, t4 );
			expH5 = toExpH(Ham.ampo, t5 );
			expH6 = toExpH(Ham.ampo, t6 );
			expH7 = toExpH(Ham.ampo, t7 );
		}
		if (ok == false)
			cout << "Error, TrotterOrder=" << o << " not implemented.\n", exit(
					0);
		cout << "exp constructed.\n";
		cout.flush();
	}
	*/


	// Output and observables
	// _______
	ofstream ent, spec, eprof, sz, h_spec, h_spec0, entr_states, energy_prof; //here I'm defining output
	ios_base::openmode mode;
	mode = std::ofstream::out; //Erase previous file (if present)

	double dt = param.val("Entropy");
	if (dt != 0) { //Entropy in the center of the chain
		ent.open("Entropy_center.dat", mode);
		ent.precision(15);
		ent
			<< "#time\tEntropy(dot)\tEntropy_SQRT_p_i(dot)\tEntr_labda_1_state\tBondDim(dot)\tMaxBondDim\n";
	}
	//---------------------
	dt = param.val("SVD_spec");
	if (dt > 0) { //SVD Spectrum on central bond
		spec.open("SVD_spec.dat", mode);
		spec.precision(15);
		spec << "#Position=" << dot << "\t<SVD_spectrum>\t\ttime\n";
	}
	//---------------------

	dt = param.val("Eprof");
	if (dt > 0) { //Full entropy profile
		eprof.open("Entropy_profile.dat", mode);
		eprof.precision(15);
		eprof << "#Position=i-" << dot << setw(16) << "\t Entropy(i)"
			<< setw(16)
			<< "\t Entropy_sqrt \t Entropy_state1 \t time \t\t Bond.Dim(i)\n";
	}
	//---------------------

	dt = param.val("Sz");
	if (dt > 0) { //Full magnetization profile
		sz.open("Sz_profile.dat", mode);
		sz.precision(15);
		sz << "#Position=i-" << "\t<Sz_i>\t" << dot << "\t\ttime\n";

	}
	//---------------------
	dt = param.val("H_spec");
	if (dt > 0) { //Hamiltonian spectum
		h_spec.open("H_spec.dat", mode);
		h_spec.precision(15);
		h_spec << "#Position=" << dot << "\t<H_spec>\t\ttime\n";

		h_spec0.open("h_spec0.dat", mode);
		h_spec0.precision(15);
		h_spec0 <<  "time and " << "\t<H_spec>\n";
	}
	//---------------------
	dt = param.val("EnergyProf");
	if (dt > 0) { //Energy profile
		energy_prof.open("Energy_profile.dat", mode);
		energy_prof.precision(15);
		energy_prof << "#Position=i-" << "\t<Ham_i>\t" << dot << "\t\ttime(or beta)\n";
	}
	//---------------------
	dt = param.val("Entr_states");
	if (dt > 0) { //Entr of the first largest states
		entr_states.open("Entropy_states.dat", mode);
		entr_states.precision(15);
		entr_states
			<< "#time\t Entr_1 \t BondDim_1 \t Entr_2 \t BondDim_2 \t  Entr_3 \t BondDim_3 \n";

	}

	// Here we define the initial time step, which will play a role of R_1

	// Time evolution
	bool go_on = true;
	for (int n = 0; n <= n_steps && go_on; ++n) {
		const double time = n * tau; //+param.val("time_shift");
		cout << "Time step #" << n << "/" << n_steps << "\ttime=" << time
			<< endl;
		;
		cout.flush();
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



		// ------- entanglement Entropy profile -----
		if (param.val("Eprof") > 0)
			if (n % int(param.val("Eprof") / tau) == 0) {
				eprof << "\"t=" << time << "\"" << endl;
				for (int i = 1; i < N; i++) {
					double entr_std = Entropy(psi, i, Myspec, 1); // p log p
					eprof << i + 0.5 - dot << "\t" << setw(16) << setfill('0')
						<< entr_std << "\t" << setw(16) << setfill('0')
						<< BondDim(psi, i) << "\t" << time << endl;
				}
				if (n < n_steps)
					eprof << endl << endl;
			}

		// ------- Sz profile -------
		if (param.val("Sz") > 0)
			if (n % int(param.val("Sz") / tau) == 0) {
				sz << "\"t=" << time << "\"" << endl;
				double sz_tot = 0, sz_left = 0, sz_right = 0, sz_dot = 0;
				for (int i = 1; i <= N; i++) {
					const double s = Sz(psi, sites, i);
					sz_tot += s;
					if (i < dot)
						sz_left += s;
					if (i > dot)
						sz_right += s;
					if (i == dot)
						sz_dot += s;
					sz << i - dot << "\t" << s << "\t" << time << endl;
				}
			}
		// ------- Energy profile -------
		if (param.val("EnergyProf") > 0){
			if (n % int(param.val("EnergyProf") / tau) == 0) {
				energy_prof << "\"t=" << time << "\"" << endl;
				for (int i = 1; i <= N-2; ++i) {
					const double en = Energy(psi, sites, i);
					energy_prof << i - dot + 1  << "\t" << en << "\t" << time << endl;
				}
				//I need this part to separate time steps in *.dat files (for gnuplot)
				if (n<= n_steps) energy_prof << "\n\n"; 
			}
		}

		if (n < n_steps) {
			//MPS psi_temp = psi;
			cout << "Time evol" << endl;
			TimeEvolution (psi, H, tau, o, args);
			psi.noPrime("Site");
			//psi.normalize();
			cout << "max bond dim = " << maxLinkDim(psi) << endl;
			cout << "U = 1+H. Norm = " << real(innerC(psi,psi))<< endl;
			cout << "U = 1+H. Energy = " << real(innerC(psi,H,psi))<< endl;

		}
		/*
		if (n < n_steps) {
			psi0 = applyMPO(expH1, psi0, args);
			if (o > 1)
				psi0 = applyMPO(expH2, psi0, args);
			if (o == 3) {
				psi0 = applyMPO(expH3, psi0, args);
				psi0 = applyMPO(expH4, psi0, args);
			}
			if (o == 4) {
				psi0 = applyMPO(expH3, psi0, args);
				psi0 = applyMPO(expH4, psi0, args);
				psi0 = applyMPO(expH5, psi0, args);
				psi0 = applyMPO(expH6, psi0, args);
				psi0 = applyMPO(expH7, psi0, args);
			}
			psi0.noPrime("Site");
			cout << "max bond dim = " << maxLinkDim(psi0) << endl;
			cout << "U = WII. Norm = "<<real(innerC(psi0,psi0))<< endl;
			cout << "U = WII. Energy = "<<real(innerC(psi0,H,psi0))<< endl;
		}
		*/
		cout << "overlap <psi|psi0> = " << (innerC(psi,psi0)) << endl;

	}
	cout << "\nTime evolution complete.\n";
	println("Done !");
	return 0;
}

