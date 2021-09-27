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
double char2double(char *a) ;

class Parameters: public map<string, double> { // class Parameters have all methods from container map<str,double>
public:
	double val(string var_name) const ;
	long longval(string var_name) const ;
	void PRint(ostream &o) const ;
	void ReadArguments(int argc, char *argv[]);
};

// class of PUBLIC parameters
class ThreeSiteParam: public Parameters {
public:
	ThreeSiteParam();
};


//I'm creating a Folded XXZ 3-site Hamiltonian
class ThreeSiteHamiltonian {
public:
	int dot;
	AutoMPO ampo;
	ThreeSiteHamiltonian(const SiteSet &sites, const ThreeSiteParam &param);
private:
	int N;
	void init(const ThreeSiteParam &param);
};



class LadderHamiltonian {
public:
	int dot;
	AutoMPO ampo;
	//Define a constructor for this class 'LadderHamiltonian'
	LadderHamiltonian(const SiteSet &sites, const ThreeSiteParam &param, const string ham_type_);
private:
	int N;
	string ham_type;
	void init(const ThreeSiteParam &param);
};


