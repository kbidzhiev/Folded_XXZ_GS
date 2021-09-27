#pragma once

#include "itensor/all.h"
#include <vector>
#include <map>
#include <string>
#include <complex>

using namespace itensor;
using namespace std;


//------------------------------------------------------------------
//The function below translate numbers (etc.) into character strings
//the second parameter (optional) is the precision (digits)

template<class T>
string to_string(const T &t, unsigned int precision = 0) {
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
	double val(string var_name) const ;
	long longval(string var_name) const ;
	void PRint(ostream &o) const ;
	void ReadArguments(int argc, char *argv[]);
};
//_____________________________________________________
// class of PUBLIC parameters
class ThreeSiteParam: public Parameters {
public:
	ThreeSiteParam();
};
//_____________________________________________________


//I'm creating a 3 site Hamiltonian for system of N sites

class ThreeSiteHamiltonian {
public:
	int dot;
	AutoMPO ampo;
	ThreeSiteHamiltonian(const SiteSet &sites, const ThreeSiteParam &param);
private:
	int N;
	void init(const ThreeSiteParam &param);
};
