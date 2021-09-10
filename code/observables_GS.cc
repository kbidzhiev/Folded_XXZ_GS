#include "itensor/all.h"
#include "observables_GS.h"

using namespace itensor;
using namespace std;



double Entropy(MPS&  psi, const int i, vector<double> &sing_vals, const double r){ //returns the von Neumann entropy on some bond (i,i+1)
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
		if(p > 1e-15){
			sing_vals[j] = p;
			j++;
			sum += p;
			SvN += -pow(p, r) * r * log(p);
		}
	}
	//cout<<"Tr(probabilities) at site(" << i << ") = "<< sum << endl;
	return SvN;
}
// Bond Dim
int BondDim(const MPS& psi, const int i) {
	auto bond_index = rightLinkIndex(psi, i);
	return (dim(bond_index));
}

// < Sz_i >
double Sz(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const int i){ //<psi|Sz|psi> at site i
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

double Sx(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const int i){ //<psi|Sz|psi> at site i
	psi.position(i);
	ITensor ket = psi(i); // read only access
	//ITensor bra = dag(prime(ket, "Site"));
	auto Sx = op(sites,"Sx", i);
	ket *= Sx;
	ket *= dag(prime(psi(i), "Site")); //multipuing by bra
	//ITensor B = ket * bra;
	double sx = real(eltC(ket)); // 2 here is "sigma_z = 2* s_z"
	return sx;
}

double Sy(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const int i){ //<psi|Sz|psi> at site i
	psi.position(i);
	ITensor ket = psi(i); // read only access
	//ITensor bra = dag(prime(ket, "Site"));
	auto Sy = op(sites,"Sy", i);
	ket *= Sy;
	ket *= dag(prime(psi(i), "Site")); //multipuing by bra
	//ITensor B = ket * bra;
	double sy = real(eltC(ket)); // 2 here is "sigma_z = 2* s_z"
	return sy;
}

double SxSx1(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const int i){
	return real(Correlation(psi, sites, "Sx", "Sx", i, i + 1));
}
double SxSx2(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const int i){
	return real(Correlation(psi, sites, "Sx", "Sx", i, i + 2));
}
double SxSz(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const int i){
	return real(Correlation(psi, sites, "Sx", "Sz", i, i + 1));
}

double SySy1(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const int i){
	return real(Correlation(psi, sites, "Sy", "Sy", i, i + 1));
}
double SySy2(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const int i){
	return real(Correlation(psi, sites, "Sy", "Sy", i, i + 2));
}
double SySz(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const int i){
	return real(Correlation(psi, sites, "Sy", "Sz", i, i + 1));
}


//< Sp_i Sm_i+4 >
complex<double> Correlation(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const string op_name1, const string op_name2, const int i, const int j) {	
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

complex<double> SzCorrelation (MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const string op_name1, const string op_name2, const int i ) {//< Sp_i Sz_i+2 Sm_i+4|> 
	// (i,i+1) (i+2,i+3) (i+4,i+5)
	auto Sp = op(sites,op_name1, i);
	auto Sz = op(sites,"Sz", i+1);
	auto Sm = op(sites,op_name2, i+2);
	ITensor ket = psi(i);
	ket *= Sp;
	auto ir = commonIndex(psi(i),psi(i+1),"Link"); // this link exist
	ket *= dag(prime(prime(psi(i), "Site"), ir));
	ket *= psi(i+1);
	ket *= Sz;
	ket *= dag(prime(prime(psi(i+1),"Site"),"Link"));
	ket *= psi(i+2);
	ket *= Sm;
	auto il = commonIndex(psi(i+2),psi(i+1),"Link");
	ket *= dag(prime(prime(psi(i+2),"Site"),il));
	complex<double> correlation = eltC(ket);
	return correlation;
}

//EgergyKin + EnergyPot at site i (i,i+2,i+4)
double Energy(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const int i) { 
	psi.position(i);
	/*
	   the "strange" coefficients 4 and 8 here appeared cause we use 
	   sigma matrices in the paper and spin-1/2 in the code.
	   sigma_j = 2 * spin_matrix_j, so
	   4 comes from XX+YY term
	   8 comes from (XX+YY)*Z

	   coeff 0.25 (0.5 in the beginning and 0.5 right before return)
	   appears cause the Hamiltonian H = 1/2 * (XX+YY)(1-Z)
	   then (XX+YY) = 0.5(SpSm+SmSp), so
	   H = 0.25(SpSm+SmSp)(1-Sz)


	 */
	double energy_kin = 2 * real (4 * 0.5 * Correlation(psi,sites, "S+", "S-", i, i+2) );
	double energy_pot = 2 * real(-8 * 0.5 * SzCorrelation(psi,sites, "S+", "S-", i ) );	
	double energy = 0.5*(energy_kin + energy_pot);
	return energy;
}


double SpSp_plus_SmSm(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const int i){
	psi.position(i);

	// SpSp + conj(SpSp) == 2 real(SpSp);
	complex<double> spsp = Correlation(psi, sites, "S+", "S+", i, i + 1);
	complex<double> smsm = Correlation(psi, sites, "S-", "S-", i, i + 1);
	return real(4 * 0.5 * (spsp + smsm));
}

// (SxSy-SySx)/2   - (SxSzSy-SySzSx)/2
double Q1minus(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const int i) { 
	psi.position(i);
	// (SxSy-SySx)/2 = (SmSp-SpSm)/4
	double q_kin = 2 * imag(4 * 0.5 * Correlation(psi,sites, "S-", "S+", i, i+2) );
	double q_pot = 2 * imag(-8 * 0.5 * SzCorrelation(psi,sites, "S-", "S+", i ) );	

	//double q_kin = real(4 * 0.5 * (
	//	Correlation(psi,sites, "Sx", "Sy", i, i+4)-
	//	Correlation(psi,sites, "Sy", "Sx", i, i+4) ));
	//double q_pot = real(-8 * 0.5 * (
	//	SzCorrelation(psi,sites, "Sx", "Sy", i )-
	//	SzCorrelation(psi,sites, "Sy", "Sx", i ) ));

	double conserved_charge_minus = -0.5*(q_kin + q_pot); // '-' is in the definition of the q1minus = -(SxSy-SySx)/2
	return conserved_charge_minus;
}

complex<double> Q1(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite> & sites, const int i) {
	psi.position(i);
	// (SxSy-SySx)/2 = (SmSp-SpSm)/4
	complex<double> q_kin = 2 * 4 * 0.5 * Correlation(psi,sites, "S-", "S+", i, i+2) ;
	complex<double> q_pot = 2 * (-8) * 0.5 * SzCorrelation(psi,sites, "S-", "S+", i ) ;

	//double q_kin = real(4 * 0.5 * (
	//	Correlation(psi,sites, "Sx", "Sy", i, i+4)-
	//	Correlation(psi,sites, "Sy", "Sx", i, i+4) ));
	//double q_pot = real(-8 * 0.5 * (
	//	SzCorrelation(psi,sites, "Sx", "Sy", i )-
	//	SzCorrelation(psi,sites, "Sy", "Sx", i ) ));

	double q1_plus  = 0.5* real(q_kin + q_pot); // '-' is in the definition of the q1minus = -(SxSy-SySx)/2
	double q1_minus = -0.5*imag(q_kin + q_pot);

	return {q1_plus, q1_minus};
}

// K = SpSm + SmSp,   D = SmSp - SpSm
complex<double> KKDD(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const int i) {
	auto Sp1 = op(sites,"S-", i+1);	//i+1 on the physical sites
	auto Sp2 = op(sites,"S-", i+2);	//i+2
	auto Sm1 = op(sites,"S+", i+3); //i+3
	auto Sm2 = op(sites,"S+", i+4); //i+4

	psi.position(i+1);
	ITensor ket = psi(i+1);
	ket *= Sp1;
	auto ir = commonIndex(psi(i+1),psi(i+2),"Link"); // this link exist
	ket *= dag(prime(prime(psi(i+1), "Site"), ir));
	ket *= psi(i+2);
	ket *= Sp2;
	ket *= dag(prime(prime(psi(i+2),"Site"),"Link"));
	ket *= psi(i+3);
	ket *= Sm1;
	ket *= dag(prime(prime(psi(i+3),"Site"),"Link"));
	ket *= psi(i+4);
	ket *= Sm2;
	auto il = commonIndex(psi(i+3),psi(i+4),"Link");
	ket *= dag(prime(prime(psi(i+4),"Site"),il));
	complex<double> kkdd = 0.25*16*2*eltC(ket); //4 is needed to convert four Spin matrices to Pauili
	return kkdd;
}

complex<double> Correlations5site(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites,
		const string op_name1,
		const string op_name2,
		const string op_name3,
		const string op_name4,
		const string op_name5,
		const int i){
	auto Op0 = op(sites, op_name1, i);
	auto Op2 = op(sites, op_name2, i+1); //i+1 on the physical sites
	auto Op4 = op(sites, op_name3, i+2); //i+2
	auto Op6 = op(sites, op_name4, i+3); //i+3
	auto Op8 = op(sites, op_name5, i+4); //i+4
	
	ITensor ket = psi(i);	
	ket *= Op0;
	auto ir = commonIndex(psi(i),psi(i+1),"Link"); // this link exist
	ket *= dag(prime(prime(psi(i), "Site"), ir));
	
	ket *= psi(i+1);
	ket *= Op2;
	ket *= dag(prime(prime(psi(i+1),"Site"),"Link"));
	
	ket *= psi(i+2);
	ket *= Op4;
	ket *= dag(prime(prime(psi(i+2),"Site"),"Link"));

	ket *= psi(i+3);
	ket *= Op6;
	ket *= dag(prime(prime(psi(i+3),"Site"),"Link"));
	
	ket *= psi(i+4);
	ket *= Op8;
	auto il = commonIndex(psi(i+3),psi(i+4),"Link");
	ket *= dag(prime(prime(psi(i+4),"Site"),il));
	
	complex<double> correlation = eltC(ket);
	return correlation;
}

complex<double> Q2(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const int i) {
	complex<double> kkdd = KKDD(psi, sites, i);
	psi.position(i);
	complex<double> zk = 0.5*8*(Correlations5site(psi,sites, "S-", "Id", "Sz", "Id", "S+",  i));// 2comes from REAL,
	// 0.5 from (spsm+smsp)/2 = sxsx+sysy,
	// 8 is 3spin to 3sigma
	complex<double> zzzk = 0.5*32*(Correlations5site(psi,sites, "S-", "Sz", "Sz", "Sz", "S+",  i));
	complex<double> zzk1 = 0.5*16*(Correlations5site(psi,sites, "S-", "Sz", "Sz", "Id", "S+",  i));
	complex<double> zzk2 = 0.5*16*(Correlations5site(psi,sites, "S-", "Id", "Sz", "Sz", "S+",  i));
	complex<double> q2 = 0.25*(kkdd + zk + zzzk - zzk1 - zzk2);
	double q2plus = -2 * real(q2);
	double q2minus = 2 * imag(q2);
	return {q2plus, q2minus};
}

double Q2plus(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const int i){
	return -real(Q2(psi, sites, i));
}

double Q2minus(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const int i){
	return imag(Q2(psi, sites, i));
}

complex<double> KSS(MPS& psi, const itensor::BasicSiteSet<itensor::SpinHalfSite>& sites, const int i){
	// KSS means kinetic(i,i+2) - sz(i) - sz(i+2);
	complex <double> kss = 0;
	kss += 4.0 * Correlation(psi, sites, "Sx", "Sx", i, i+2);
	kss += 4.0 * Correlation(psi, sites, "Sy", "Sy", i, i+2);
	kss += - 2.0 * (Sz(psi, sites, i) + Sz(psi, sites, i + 2));
	kss += 2.0;
	return kss;
}



