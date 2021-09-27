#include "itensor/all.h"
#include "observables_GS.hpp"

using namespace itensor;
using namespace std;








MPS Measure(MPS& psi, const SiteSet &sites, const string op_name, const int j, Args args){

	AutoMPO Op_ampo(sites);

	if (op_name == "Energy"){
		Op_ampo +=  4 * 0.25, "S+", j, "S-", j + 2; // 0.5 (SpSm+ SmSp) = SxSx + SySy
		Op_ampo +=  4 * 0.25, "S-", j, "S+", j + 2;
		Op_ampo += -8 * 0.25, "S+", j, "Sz", j + 1, "S-", j + 2;
		Op_ampo += -8 * 0.25, "S-", j, "Sz", j + 1, "S+", j + 2;

		Op_ampo +=  4 * 0.25, "S+", j + 1, "S-", j + 3; // 0.5 (SpSm+ SmSp) = SxSx + SySy
		Op_ampo +=  4 * 0.25, "S-", j + 1, "S+", j + 3;
		Op_ampo += -8 * 0.25, "S+", j + 1, "Sz", j + 2, "S-", j + 3;
		Op_ampo += -8 * 0.25, "S-", j + 1, "Sz", j + 2, "S+", j + 3;
	} else if (op_name == "Sz"){
		Op_ampo += 1.0, "Sz", j;
		Op_ampo += 1.0, "Sz", j + 1;
		Op_ampo += 1.0, "Sz", j + 2;
		Op_ampo += 1.0, "Sz", j + 3;
	} else if (op_name == "Staggered_Sz"){
		Op_ampo += 1.0, "Sz", j;
		Op_ampo += -1.0, "Sz", j + 1;
		Op_ampo += 1.0, "Sz", j + 2;
		Op_ampo += -1.0, "Sz", j + 3;
	}

	auto Op = toMPO(Op_ampo);
	psi = applyMPO(Op,psi,args);
	return psi;

}
