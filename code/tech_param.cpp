#include "itensor/all.h"
#include "time_evolution.hpp"

using namespace itensor;
using namespace std;




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
double Parameters::val(string var_name) const { // .val("C") gives value for parameter C
	map<string, double>::const_iterator it = find(var_name);
	if (it == end()) {
		cout << "Error: Parameter " << var_name << " is not defined.\n", exit(0);
		return 0;
	} else
		return it->second;
}
long Parameters::longval(string var_name) const {
	double v = val(var_name);
	if (abs(double(round(v)) - v) < 1e-6) {
		return long(round(v));
	} else {
		cout << "Error, parameter " << var_name << "=" << v << " is not a long"
				<< endl, exit(0);
		return 0;
	}
}
void Parameters::PRint(ostream &o) const {
	for (map<string, double>::const_iterator it = begin(); it != end(); ++it) {
		o << it->first << "=" << it->second << endl;
	}
}

void Parameters::ReadArguments(int argc, char *argv[]) {
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

//_____________________________________________________
// class of PUBLIC parameters
ThreeSiteParam::ThreeSiteParam() { //Constructor
	//Specify below all the allowed parameter names,
	//and their default values
	operator[]("N") = 10; //Length of the chain
	operator[]("J") = 1.0; // coupling
	operator[]("Jy") = 1.0; // coupling in the FoldedXYZ
	operator[]("Delta") = 0; // Delta in XYZ and XXZ
	operator[]("J2") = 0; // Integrability breaking term

	operator[]("tau") = 0.02;  //time step for the unitary evolution
	operator[]("T") = 2;  //Total (final) time
	operator[]("Sz") = 0.1;
	operator[]("SVD_spec") = 0; //SVD spectrum
	operator[]("max_bond") = 4000;  //maximum bond dimension
	operator[]("trunc") = 1e-8;  //maximum truncation error
	operator[]("energy") = 1e-13;  //convergence criterium on the energy
	operator[]("sweeps") = 999;  //maximum number of sweeps in the DMRG
	operator[]("TrotterOrder") = 2;

	//Initial state protocols. One should chose ONE of them to be > 0
	operator[]("GroundState") = 0;
	operator[]("LadderState") = 0;
	operator[]("Neel") = 0;
	operator[]("UUD") = 0;
	operator[]("JammedImpurity") = 0;
	operator[]("PPK") = 0;
	operator[]("PPX") = 0;
	operator[]("XXZ") = 0;
	operator[]("XXZGlobal") = 0;
	operator[]("XXZDW") = 0;
	operator[]("XP") = 0;
	operator[]("FoldedXYZ") = 0;


	operator[]("alpha") = 0; // U = exp[ i alpha  n*\sigma]
	operator[]("begin") = 1;
	operator[]("hL") = 0; //initial magnetization_L
	operator[]("hR") = 0; //initial magnetization_R
	operator[]("h") = 0.0;
	operator[]("rho") = 0.0;
	operator[]("n") = 1;

	// Observables to be measured
	operator[]("Q1Profile") = 0; // energy and current profile
	operator[]("Q2Profile") = 0;
	operator[]("Current") = 0;
	operator[]("Entropy") = 0; //entanglement entropy p*log*p between left and right parts of system
	operator[]("EntropyProfile") = 0; // Entropy Profile- parameter 0 -> nothing, dt>0 each second=integer parameter
	operator[]("Loschmidt") = 0; // loschmidt echo <psi(t)|psi(0)>


	operator[]("Dhar") = 0; // Deepak Dhar term in the Hamiltonian (time evolution ONLY)
	operator[]("PXXP") = 0; //Integrability breaking term
	operator[]("Measurement") = 0; // Project the central spin to be |Down> (time evolution ONLY)
	operator[]("TrotterOrderXXZ") = 3;



	operator[]("Distance") = 5;
}



