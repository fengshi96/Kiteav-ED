/* ***************************************************
 *  Authors:Chris Bishop & Nirav D. Patel
 *	Date:January 14, 2016
 *	Platform: linux
 *  Engine Modified by: Shi Feng
 *  	Date:Mar 24 2020 
 *  Language: C++
 *	Model: 3OrMC_MF-Test
 * ***************************************************
*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <omp.h>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <random>
#include <complex>
#include <cmath>
#include <cassert>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

using namespace std;
//#include "Vector.h"
#include "Memusage.h"
#include "Matrix.h"
#include "ParametersEngine.h"
#include "Lattice.h"
#include "QuantumBasis.h"
#include "Hamiltonian.h"
#include "Exactdiag.h"
#include "DynLanczosSolver.h"
#include "Observables.h"

typedef std::size_t SizeType;
typedef Eigen::VectorXf VectorXf;
typedef ConstVariables ConstVariablestype;
typedef QBasis Basistype;
typedef Hamiltonian HamilType;
typedef DynLanczosSolver DynLancType;
typedef Observables ObservType;

typedef std::complex<double> dcomplex;    // your typedef
void PrintGSSz(Eigen::VectorXcd& Psi, Basistype& Basis, int& Nsite_);
void PrintGSConf(Eigen::VectorXcd& Psi, Basistype& Basis, int& Nsite_);
//void PrintGSConf(VectorXf& Psi, Basistype& Basis, int& Nsite_, int& basisLabel1);

int main(int argc, char *argv[]) {
	if (argc<2) {
		throw std::invalid_argument("USE:: executable inputfile");
	}

	string inputfile = argv[1];
	double pi=acos(-1.0);
    double vm, rss;
    dcomplex zero(0.0,0.0);
    dcomplex onei(0.0,1.0);

	//std::cout.unsetf ( std::ios::floatfield );
	//std::cout.precision(6);

	// Read from the inputfile
	ConstVariablestype Parameters(inputfile);
        omp_set_dynamic(0);
        omp_set_num_threads(Parameters.Threads);
	Lattice Lat(Parameters);
	int Nsite_ = Parameters.NumberofSites;

//	// Build binary basis
    Basistype Basis(Parameters);
    int nn=Basis.basis.size();
    process_mem_usage(vm, rss);
    std::cout << "       VM (MB): " << int(vm/1024) << "       RSS (MB): " << int(rss/1024) << std::endl;


//	// Build non-zero Diagonal and Tight-binding part of Hamiltonian
    HamilType Hamiltonian(Parameters,Lat,Basis);
    Hamiltonian.TightB_Ham();
    process_mem_usage(vm, rss);
    std::cout << "       VM (MB): " << int(vm/1024) << "       RSS (MB): " << int(rss/1024) << std::endl;

    Hamiltonian.Diagonal();
    process_mem_usage(vm, rss);
    std::cout << "       VM (MB): " << int(vm/1024) << "       RSS (MB): " << int(rss/1024) << std::endl;


   ObservType Observ(Parameters, Lat, Basis, Hamiltonian);
   DynLancType Lanczos(Parameters, Basis, Hamiltonian);
   Eigen::VectorXcd Psi;
   if(Parameters.Solver=="ED") {
       Psi = ExactDiag(Parameters, Basis, Hamiltonian, Observ);
   } else {
       Psi = Lanczos.Lanczos_Nirav(Parameters.LancType);
       cout << setprecision(6) << fixed;
       cout << " Ritz values: ";
       for(int i=0;i<Lanczos.TriDeval.size();i++) cout << Lanczos.TriDeval[i] << " ";
       cout << endl;
   }

//    cout << setprecision(6) << fixed;
//    cout << " Ritz values: ";
//    Eigen::VectorXcd Sxi, Syi, Szi, Spi, Smi;
//    Observ.measureLocalS(Psi, Sxi, Syi, Szi, Spi, Smi);

//    std::vector<Eigen::MatrixXcd> A(9);
//    A[0] = Observ.TwoPointCorr("Sx","Sx",Psi); cout << A[0] << flush << endl << endl;
//    A[1] = Observ.TwoPointCorr("Sx","Sy",Psi); cout << A[1] << flush << endl << endl;
//    A[2] = Observ.TwoPointCorr("Sx","Sz",Psi); cout << A[2] << flush << endl << endl;

//    A[3] = Observ.TwoPointCorr("Sy","Sx",Psi); cout << A[3] << flush << endl << endl;
//    A[4] = Observ.TwoPointCorr("Sy","Sy",Psi); cout << A[4] << flush << endl << endl;
//    A[5] = Observ.TwoPointCorr("Sy","Sz",Psi); cout << A[5] << flush << endl << endl;

//    A[6] = Observ.TwoPointCorr("Sz","Sx",Psi); cout << A[6] << flush << endl << endl;
//    A[7] = Observ.TwoPointCorr("Sz","Sy",Psi); cout << A[7] << flush << endl << endl;
//    A[8] = Observ.TwoPointCorr("Sz","Sz",Psi); cout << A[8] << flush << endl << endl;

//    cout << " --> Spin total = ";
//    dcomplex Stotal = A[0].sum() + A[4].sum() + A[8].sum();
//    cout << Stotal << endl;

//    cout << " --> Spin trace = ";
//    Stotal = A[0].trace() + A[4].trace() + A[8].trace();
//    cout << Stotal << endl << endl;


//    cout << " --> Pairs = SxSx - SySy + iSxSy + iSySx \n";
//    //Eigen::MatrixXcd AdagAdag = - A[0] + A[4] + dcomplex(0,1)*A[1] + dcomplex(0,1)*A[3];
//    Eigen::MatrixXcd AdagAdag = A[0] - A[4] + dcomplex(0,1)*A[1] + dcomplex(0,1)*A[3];
//    cout << AdagAdag << endl << endl;

//    Eigen::VectorXcd pairx(Nsite_), pairy(Nsite_), pairz(Nsite_);
//    for(int i=0; i<Nsite_;i++) {
//        int n1x = Lat.N1neigh_(i,1);
//        int n1y = Lat.N1neigh_(i,2);
//        int n1z = Lat.N1neigh_(i,0);

//        pairx(i) = AdagAdag(i,n1x) - Spi[i] *Spi[n1x];
//        pairy(i) = AdagAdag(i,n1y) - Spi[i] *Spi[n1y];
//        pairz(i) = AdagAdag(i,n1z) - Spi[i] *Spi[n1z];

//        cout << i << " " << pairx(i) << " " << pairy(i) << " " << pairz(i) << " \n";
//    }
//    cout << endl << endl;
//    cout << "total pair x,y,z = " << real(pairx.sum())/Nsite_ << " " << imag(pairx.sum())/Nsite_ << " "
//                                  << real(pairy.sum())/Nsite_ << " " << imag(pairy.sum())/Nsite_ << " "
//                                  << real(pairz.sum())/Nsite_ << " " << imag(pairz.sum())/Nsite_ << endl << endl;





//    // -------------------------------
//    // --- Topological Ent Entropy ---
//    // -------------------------------.
//    cout << endl << endl;
//    cout << "Calculating the topological entanglement entropy " << endl;


//    // -- Define a cut;
//    Eigen::VectorXi S1sites(12);  // typedef Matrix< int , Dynamic , 1> Eigen::VectorXi
//    Eigen::VectorXi S2sites(12);
//    S1sites << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11;
//    S2sites << 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23;
//    //S1sites << 0, 7;
//    //S2sites << 1, 2, 3, 4, 5, 6;
//    Eigen::MatrixXi ABCutLabels = Basis.ABSplit(S1sites, S2sites);

//    //cout << ABCutLabels << endl;
//    cout << " --> Print ABCutLabels into ABCutLabels.dat";
//    ofstream outfile;
//    string outfilename1 = "ABCutLabels.dat";
//    outfile.open(outfilename1);
//    outfile << ABCutLabels;
//    outfile.close();

//    int S1Hil = pow(2,S1sites.size());
//    int S2Hil = pow(2,S2sites.size());

//    Eigen::MatrixXcd S1RDM(S1Hil,S1Hil);
//    S1RDM.setZero();
//    for (int i=0; i<S1Hil; i++) {
//        for (int j=0; j<S1Hil; j++) {
//            for (int k=0; k<S2Hil; k++) {
//                SizeType ikL = ABCutLabels(i,k);
//                SizeType jkL = ABCutLabels(j,k);
//                S1RDM(i,j) += conj(Psi[ikL])*Psi[jkL];
//            }
//        }
//    }
//    cout << " --> Print reduced density matrix into S1RDM.dat\n";
//    //cout << S1RDM << endl << endl;
//    string outfilename2 = "S1RDM.dat";
//    outfile.open(outfilename2);
//    outfile << S1RDM;
//    outfile.close();

//    cout << " --> Trace of the reduced density matrix " << S1RDM.trace() << endl;



//    Eigen::MatrixXcd tmp = S1RDM;
//    Eigen::VectorXd D;
//    Diagonalize('V',tmp,D);
//    cout << " --> eigen values of the RDM = ";
//    for(int i=0;i<D.size();i++) cout << D[i] << ", ";
//    cout << endl;


//    double entropy=0.0;
//    for(int i=0;i<D.size();i++) {
//        if(abs(D[i])>1e-6) {
//            entropy += -D[i]*log(D[i]);
//        }
//    }
//    cout << " --> vN Entropy = " << entropy << endl << endl;



//	cout << "--------THE END--------" << endl;
}


/*=======================================================================
 * ======================================================================
 * ======================================================================
*/

void PrintGSSz(Eigen::VectorXcd& Psi, Basistype& Basis, int& Nsite_) {
    int k_maxup = Basis.basis.size();
    int n=k_maxup;

    Eigen::VectorXcd V(Nsite_+1);
    for(int i0=0;i0<n;i0++){
        int ket = Basis.basis[i0];
        dcomplex coef1 = Psi[ket];
        int SzT = Basis.NPartInState(ket,Nsite_);
        V[SzT] += coef1; //*conj(coef1);
    }
    cout << V << endl;
}


void PrintGSConf(Eigen::VectorXcd& Psi, Basistype& Basis, int& Nsite_) {
	int k_maxup = Basis.basis.size();
	int n=k_maxup;

	vector<pair<double,int> >V;
	for(int i=0;i<n;i++){
		pair<double,int>P=make_pair((Psi[i]*conj(Psi[i])).real(),i);
		V.push_back(P);
	}
	sort(V.begin(),V.end());

	for(int i0=n-1;i0>=n-10;i0--){
		int i = V[i0].second;
		int basisLabel1 = i; // i1*k_maxdn + i2;
		int ket = Basis.basis[i];
		double coef = V[i0].first; //Psi[basisLabel1];
		dcomplex coef1 = Psi[basisLabel1];

		cout << basisLabel1 << " - ";
		//if(coef*coef>0.02) {
		for(int s=0;s<Nsite_;s++){
			cout << " \\underline{";
			if(Basis.findBit(ket,s) == 1) {
				cout << " \\uparrow ";
			} else {
				cout << " \\downarrow ";
			}
			cout << "} \\ \\ ";
			if(s==Nsite_-1) cout << " \\ \\ ";
		}
		cout << "& \\ \\ ";
		cout << coef << " \\\\ " << coef1 << endl;
	}
	cout << "   " << endl;

}

/*

void PrintGSConf(VectorXf& Psi, Basistype& Basis, int& Nsite_, int& basisLabel1) {
	int k_maxup = Basis.upbasis.size();
	int k_maxdn = Basis.downbasis.size();
	int n=k_maxup*k_maxdn;

	for(int i1=0; i1<k_maxup; i1++) {
		int up_ket1 = Basis.upbasis[i1];
		for(int i2=0; i2<k_maxdn; i2++) {
			int dn_ket1 = Basis.downbasis[i2];
			int basisLabel = i1*k_maxdn + i2;

			if(basisLabel==basisLabel1) {
				double coef = Psi[basisLabel1];

				cout << basisLabel1 << " - ";
				//if(coef*coef>0.02) {
				for(int s=0;s<Nsite_;s++){
					cout << " \\underline{";
					if(Basis.findBit(up_ket1,s) == 1) {
						cout << " \\uparrow ";
					} else {
						cout << " \\ ";
					}
					if(Basis.findBit(dn_ket1,s) == 1) {
						cout << " \\downarrow ";
					} else {
						cout << " \\ ";
					}
					cout << "} \\ \\ ";
					if(s==Nsite_-1) cout << " \\ \\ ";
				}
				cout << "& \\ \\ ";

//				if(coef<0) {
//					cout << " (-) \\ ";
//				} else {
//					cout << " (+) \\ ";
//				}
				cout << coef << " \\\\ " << endl;
			}
		}
	}

	//}
	//		}}
	cout << "   " << endl;

}




void LaTexPrint(int& upstate, int& dnstate, double& coef, int& Nsite_, Basistype& Basis) {
	for(int s=0;s<Nsite_;s++){
		cout << " \\underline{";
		if(Basis.findBit(upstate,s) == 1) {
			cout << " \\uparrow ";
		} else {
			cout << " \\ ";
		}
		if(Basis.findBit(dnstate,s) == 1) {
			cout << " \\downarrow ";
		} else {
			cout << " \\ ";
		}
		cout << "} \\ \\ ";
		if(s==Nsite_-1) cout << " \\ \\ ";
	}
	cout << "& \\ \\ ";

	if(coef<0) {
		cout << " (-) \\ ";
	} else {
		cout << " (+) \\ ";
	}

	cout << coef*coef << " \\\\ " << endl;

}



*/





/*
-
-
-	{
-		cout << " ==================================================================== " << endl;
-		cout << "                            Finite Temperature						   " << endl;
-		cout << " ==================================================================== " << endl;
-		int Nsite_ = Parameters.NumberofSites;
-		Eigen::VectorXcd phi0 = Lanczos.Vorig;
-
-		double EvTemp_Scale=1; //11605;
-		int ExtState = Lanczos.PsiAll.size();
-		int Nst = int(ExtState*1.0-1);
-		int M = Nst;
-		cout << Nst << " \t " << endl;
-		cout << " energy span: " << Lanczos.TriDeval[M] - Lanczos.TriDeval[0] << endl;
-
-
-		ofstream myfile;
-		myfile.open ("EvsTemp1.dat");
-		// -- partition function --- J. Jaklic and P. Prelovsek - page 14
-
-		Eigen::VectorXcd Hoverlap(M), Hsqoverlap(M), jstOlap(M);
-		Eigen::VectorXcd SxToverlap(M), SyToverlap(M), SzToverlap(M);
-		Eigen::VectorXcd SxTToverlap(M), SyTToverlap(M), SzTToverlap(M);
-
-		Eigen::VectorXcd Hphi0 = Hamiltonian.MatrixVectorMult(phi0);
-		Eigen::VectorXcd HHphi0 = Hamiltonian.MatrixVectorMult(Hphi0);
-		Eigen::VectorXcd Sxphi0 = Observ.ApplySx(phi0,Parameters.SpecialSite);
-		Eigen::VectorXcd Syphi0 = Observ.ApplySy(phi0,Parameters.SpecialSite);
-		Eigen::VectorXcd Szphi0 = Observ.ApplySz(phi0,Parameters.SpecialSite);
-
-		Eigen::VectorXcd SxTotalp0 = Observ.ApplySx(phi0,0);
-		Eigen::VectorXcd SyTotalp0 = Observ.ApplySy(phi0,0);
-		Eigen::VectorXcd SzTotalp0 = Observ.ApplySz(phi0,0);
-		for(int i=1; i<Nsite_;i++) {
-			SxTotalp0 = SxTotalp0 + Observ.ApplySx(phi0,i);
-			SyTotalp0 = SyTotalp0 + Observ.ApplySy(phi0,i);
-			SzTotalp0 = SzTotalp0 + Observ.ApplySz(phi0,i);
-		}
-
-		//Observ.ApplySzTotal(phi0,false);
-		//Szphi0 += Observ.ApplySz(phi0,Parameters.SpecialSite+1);
-		//Eigen::VectorXcd SzT_phi0 = Observ.ApplySzTotal(Szphi0,false);
-
-#pragma omp parallel for
-		for (int j=0; j<M; j++){
-			Hoverlap[j] = Lanczos.PsiAll[j].dot(Hphi0);
-			Hsqoverlap[j] = Lanczos.PsiAll[j].dot(HHphi0);
-
-			SxToverlap[j] = Lanczos.PsiAll[j].dot(Sxphi0);
-			SyToverlap[j] = Lanczos.PsiAll[j].dot(Syphi0);
-			SzToverlap[j] = Lanczos.PsiAll[j].dot(Szphi0);
-
-			SxTToverlap[j] = Lanczos.PsiAll[j].dot(SxTotalp0);
-			SyTToverlap[j] = Lanczos.PsiAll[j].dot(SyTotalp0);
-			SzTToverlap[j] = Lanczos.PsiAll[j].dot(SzTotalp0);
-
-			jstOlap[j] = phi0.dot(Lanczos.PsiAll[j]);
-			cout << j << "    ";
-			cout.flush();
-			//cout << r_jOlap << " \t " << J_HrOlap << " \t " << J_HrOlap*r_jOlap << endl;
-		}
-
-
-		int temax=1000000;
-		double dT=0.00002;
-		double Temp=dT;
-		for(int te=0;te<=temax;te++) {
-			//double Temp = Parameters.Temperature/11605;
-			//double Temp = 0.000002 + te*0.000002;
-			if(Temp>0.05 && Temp<=1) dT=0.0005;
-			if(Temp>1 && Temp<=5) dT=0.005;
-			if(Temp>5) dT=0.05;
-			Temp = Temp + dT;
-			double beta = 1.0/Temp;
-
-			double Energy=0,  Partition=0, Ensq=0, SxT=0, SyT=0, SzT=0;
-			double SxTT=0, SyTT=0, SzTT=0;
-			double expon;
-#pragma omp parallel for private(expon) reduction(+:Partition,Energy,Ensq,SxT,SyT,SzT,SxTT,SyTT,SzTT)
-			for(int j=0; j<M; j++) {
-				expon = exp(-1.0*beta*(Lanczos.TriDeval[j] - Lanczos.TriDeval[0]));
-				Partition+= (expon*jstOlap[j]*jstOlap[j]).real();
-				Energy+= (expon*jstOlap[j]*Hoverlap[j]).real();
-				Ensq+= (expon*jstOlap[j]*Hsqoverlap[j]).real();
-
-				SxT+= (expon*jstOlap[j]*SxToverlap[j]).real();
-				SyT+= (expon*jstOlap[j]*SyToverlap[j]).real();
-				SzT+= (expon*jstOlap[j]*SzToverlap[j]).real();
-
-				SxTT+= (expon*jstOlap[j]*SxTToverlap[j]).real();
-				SyTT+= (expon*jstOlap[j]*SyTToverlap[j]).real();
-				SzTT+= (expon*jstOlap[j]*SzTToverlap[j]).real();
-			}
-
-			double Cv = (Ensq/Partition - (Energy/Partition)*(Energy/Partition))/(Temp*Temp);
-
-			if(te%100==0) {
-				cout << Temp*EvTemp_Scale << " \t " << Partition
-				     << " \t " << Energy
-				     << " \t " << Energy/(Partition*Nsite_)
-				     << " \t " << Ensq/Partition
-				     << " \t " << Cv/Nsite_
-				     << " \t " << (SxT/Partition)
-				     << " \t " << (SyT/Partition)
-				     << " \t " << (SzT/Partition)
-				     << " \t " << (SxTT/Partition)
-				     << " \t " << (SyTT/Partition)
-				     << " \t " << (SzTT/Partition)
-				     << endl;
-			}
-
-			myfile << Temp*EvTemp_Scale << " \t " << Partition
-			       << " \t " << Energy
-			       << " \t " << Energy/(Partition*Nsite_)
-			       << " \t " << Ensq/Partition
-			       << " \t " << Cv/Nsite_
-			       << " \t " << (SxT/Partition)
-			       << " \t " << (SyT/Partition)
-			       << " \t " << (SzT/Partition)
-			       << " \t " << (SxTT/Partition)
-			       << " \t " << (SyTT/Partition)
-			       << " \t " << (SzTT/Partition)
-			       << endl;
-
-			if(Temp*EvTemp_Scale>100) break;
-		}
-		myfile.close();
-		cout << endl;
-
-		// =============== spin - susceptibility ======================
-		cout << " ==================================================================== " << endl;
-		cout << "                            spin - susceptibility					   " << endl;
-		cout << " ==================================================================== " << endl;
-
-		Eigen::MatrixXcd BSz(M,M), BSx(M,M), BSy(M,M);
-		BSz.setZero();
-		BSx.setZero();
-		BSy.setZero();
-
-		for (int i=0; i<M; i++){
-#pragma omp parallel for
-			for (int j=0; j<M; j++) {
-				BSx(i,j) = Lanczos.PsiAll[i].dot(Observ.ApplySx(Lanczos.PsiAll[j],0));
-				BSy(i,j) = Lanczos.PsiAll[i].dot(Observ.ApplySy(Lanczos.PsiAll[j],0));
-				BSz(i,j) = Lanczos.PsiAll[i].dot(Observ.ApplySz(Lanczos.PsiAll[j],0));
-			}
-			cout << i << "    ";
-			cout.flush();
-			//cout << r_jOlap << " \t " << J_HrOlap << " \t " << J_HrOlap*r_jOlap << endl;
-		}
-
-
-		Eigen::MatrixXcd ASz(M,Nsite_), ASx(M,Nsite_), ASy(M,Nsite_);
-		ASz.setZero();
-		ASx.setZero();
-		ASy.setZero();
-
-		for(int s=0; s<Nsite_; s++) {
-			VectorXcd tmpx = Observ.ApplySx(phi0,s);
-			VectorXcd tmpy = Observ.ApplySy(phi0,s);
-			VectorXcd tmpz = Observ.ApplySz(phi0,s);
-
-#pragma omp parallel for
-			for (int j=0; j<M; j++) {
-				ASx(j,s) = Lanczos.PsiAll[j].dot(tmpx);
-				ASy(j,s) = Lanczos.PsiAll[j].dot(tmpy);
-				ASz(j,s) = Lanczos.PsiAll[j].dot(tmpz);
-			}
-
-		}
-
-
-		myfile.open ("SzSzvsTemp.dat");
-		temax=1000000;
-		dT=0.00002;
-		Temp=dT;
-		for(int te=0;te<=temax;te++) {
-			if(Temp>0.05 && Temp<=1) dT=0.0005;
-			if(Temp>1 && Temp<=5) dT=0.005;
-			if(Temp>5) dT=0.05;
-			Temp = Temp + dT;
-			double beta = 1.0/Temp;
-
-
-
-			if(te%10==0) {
-				double SxT=0, SyT=0, SzT=0;
-				double expon;
-
-#pragma omp parallel for private(expon) reduction(+:SxT,SyT,SzT)
-				for(int s=0; s<Nsite_;s++) {
-					for(int i=0; i<M; i++) {
-						expon = exp(-1.0*beta*(Lanczos.TriDeval[i] - Lanczos.TriDeval[0]));
-						for(int j=0; j<M; j++) {
-							SxT+= (expon*jstOlap[j]*BSx(i,j)*ASx(j,s)).real();
-							SyT+= (expon*jstOlap[j]*BSy(i,j)*ASy(j,s)).real();
-							SzT+= (expon*jstOlap[j]*BSz(i,j)*ASz(j,s)).real();
-						}
-					}
-				}
-
-				double Partition=0;
-#pragma omp parallel for private(expon) reduction(+:Partition)
-				for(int j=0; j<M; j++) {
-					expon = exp(-1.0*beta*(Lanczos.TriDeval[j] - Lanczos.TriDeval[0]));
-					Partition+= (expon*jstOlap[j]*jstOlap[j]).real();
-				}
-
-				if(te%100==0) {
-					cout << Temp*EvTemp_Scale << " \t " << Partition
-					     << " \t " << (SxT/Partition)
-					     << " \t " << (SyT/Partition)
-					     << " \t " << (SzT/Partition)
-					     << " \t " << ((SxT+SyT+SzT)/Partition)
-					     << endl;
-				}
-
-				myfile << Temp*EvTemp_Scale << " \t " << Partition
-				       << " \t " << (SxT/Partition)
-				       << " \t " << (SyT/Partition)
-				       << " \t " << (SzT/Partition)
-				       << " \t " << ((SxT+SyT+SzT)/Partition)
-				       << endl;
-
-
-			}
-
-			if(Temp*EvTemp_Scale>100.0) break;
-		}
-		myfile.close();
-
-
-
-
-
-
-		exit(1);

*/







