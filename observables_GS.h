#pragma once
#include "itensor/all.h"


using namespace itensor;
using namespace std;

double Entropy(MPS& psi, int i, vector<double> &sing_vals, double r){ //returns the von Neumann entropy on some bond (i,i+1)
	//index linking j to j+1
	auto bond_index = rightLinkIndex(psi, i); //=commonIndex(psi.A(i),psi.A(i+1),Link);
	int bond_dim = dim(bond_index);

	psi.position(i);
	ITensor wf = psi(i) * psi(i + 1);
	auto U = psi(i);
	ITensor S, V;
	//Remark: We know that the rank of wf is at most bond_dim, so we specify
	//this value to the SVD routine, in order to avoid many spurious small singular values (like ~ 1e-30)
	//which should in fact be exaclty zero.
	auto spectrum = svd(wf, U, S, V, { "MaxDim", bond_dim }); // Todo: we should use min(bond_dim, bond_dim_left*2, bond_dim_right*2)
	Real SvN = 0.;
	Real sum = 0;
	//cout<<"\tSingular value decomp.:"<<endl;
	//cout<<"\t\tdim="<<spectrum.numEigsKept()<<endl;
	sing_vals.resize(spectrum.numEigsKept());
	//cout<<"\t\tLargest sing. val:"<<spectrum.eig(1);
	//cout<<"\tSmallest sing. val:"<<spectrum.eig(spectrum.numEigsKept())<<endl;
	int j = 0;
	for (double p : spectrum.eigs()) { // auto p give small values like 1e-322

		sing_vals[j] = p;
		j++;
		sum += p;
		SvN += -pow(p, r) * r * log(p);

	}
	//cout<<"\t\tsum ="<<sum<<endl;
	return SvN;
}
// Bond Dim
int BondDim(const MPS& psi, int i) {
	auto bond_index = rightLinkIndex(psi, i);
	return (dim(bond_index));
}

// Sz on a site i

double Sz(MPS& psi, const auto& sites, const int i){ //<psi|Sz|psi> at site i
	psi.position(i);
	ITensor ket = psi(i); // read only access
	//ITensor bra = dag(prime(ket, "Site"));
	auto Sz = op(sites,"Sz", i);
	ket *= Sz;
	ket *= dag(prime(psi(i), "Site")); //multipuing by bra
	//ITensor B = ket * bra;
	double sz = real(eltC(ket)); // 2 here is "sigma_z = 2* s_z"
	return sz;
}

complex<double> Correlation(MPS& psi, const auto& sites, const string op_name1, const string op_name2, const int i, const int j) {//< Sp_i Sm_i+2 >	
	ITensor ket = psi(i);
	auto Sp = op(sites, op_name1, i);
	auto Sm = op(sites, op_name2, j);
	ket *= Sp;
	auto ir = commonIndex(psi(i),psi(i+1),"Link");
	ket *= dag(prime(prime(psi(i), "Site"), ir));
	for(int k = i + 1; k < j; ++k){
		ket *= psi(k);
		auto right = commonIndex(psi(k),psi(k+1),"Link");
		auto left  = commonIndex(psi(k-1),psi(k),"Link");
		ket *= dag(prime(prime(psi(k), right), left));
	}
	ket *= psi(j);
	ket *= Sm;
	auto il = commonIndex(psi(j),psi(j-1),"Link");
	ket *= dag(prime(prime(psi(j),"Site"),il));
	complex<double> correlation = eltC(ket);
	return correlation;
}

complex<double> SzCorrelation (MPS& psi, const auto& sites, const string op_name1, const string op_name2, const int i ) {//< Sp_i Sz_i+1 Sm_i+2 > 
	// (i,i+1,i+2)
	ITensor ket = psi(i);
	auto Sp = op(sites,op_name1, i);
	auto Sz = op(sites,"Sz", i+1);
	auto Sm = op(sites,op_name2, i+2);
	//psi(i)*Sp
	ket *= Sp;
	auto ir = commonIndex(psi(i),psi(i+1),"Link"); // this link exist
	// psi(i)*Sp*psi(i)
	ket *= dag(prime(prime(psi(i), "Site"), ir));
	// psi(i)*Sp*psi(i) * psi(i+1)
	ket *= psi(i+1);
	// psi(i)*Sp*psi(i) * psi(i+1)*Sz
	ket *= Sz;
	// psi(i)*Sp*psi(i) * psi(i+1)*Sz*psi(i+1)
	ket *= dag(prime(prime(psi(i+1),"Site"),"Link"));
	// psi(i)*Sp*psi(i) * psi(i+1)*Sz*psi(i+1) * psi(i+2)
	ket *= psi(i+2);
	// psi(i)*Sp*psi(i) * psi(i+1)*Sz*psi(i+1) * psi(i+2)*Sm
	ket *= Sm;
	auto il = commonIndex(psi(i+1),psi(i+2),"Link");
	// psi(i)*Sp*psi(i) * psi(i+1)*Sz*psi(i+1) * psi(i+2)*Sm*psi(i+2) 
	ket *= dag(prime(prime(psi(i+2),"Site"),il));
	complex<double> correlation = eltC(ket);
	return correlation;
}

double Energy(MPS& psi, const auto& sites, const int i) { //EgergyKin + EnergyPot at site i (i,i+2,i+4)
	//4 is here because 2S_+ = \sigma_+, 2S_- = \sigma_-
	psi.position(i);
	complex <double> energy_kin =  4 * 0.25 * ( Correlation(psi,sites, "S+", "S-", i, i+2) +
			Correlation(psi, sites, "S-", "S+", i, i+2) );
	complex <double> energy_pot = -8 * 0.25 * ( SzCorrelation(psi,sites, "S+", "S-", i ) +
			SzCorrelation(psi, sites, "S-", "S+", i ));
	double energy = real(energy_kin + energy_pot);
	return energy;
}


