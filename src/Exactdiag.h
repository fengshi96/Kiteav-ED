#ifndef EXACTDIAG_H
#define EXACTDIAG_H
#include "Observables.h"

typedef Observables ObservType;
typedef std::complex<double> dcomplex;
typedef Eigen::VectorXf VectorXf;
typedef Eigen::VectorXcd VectorXcd;

extern "C" void zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
                       std::complex<double> *,int *, double *, int *);

// ===== define functions ==========
Eigen::VectorXcd ReturnRow(Matrix<dcomplex>& mat, int row);
Eigen::VectorXcd ReturnCol(Matrix<dcomplex>& mat, int col);
VectorXf ReturnRow(Matrix<double>& mat, int row);
VectorXf ReturnCol(Matrix<double>& mat, int col);

Matrix<dcomplex> SpinSpinCorr(ObservType& Observ_, Eigen::VectorXcd& Psi, int& Nsite_);
Eigen::VectorXcd ExactDiag(ConstVariables& variables_, QBasis& Basis_, Hamiltonian& Hamil_, ObservType& Observ_);
void Diagonalize(char option, Matrix<dcomplex>& A, vector<double>& eval);
void FiniteTemperature(ConstVariables& variables_, QBasis& Basis_,
                       Hamiltonian& Hamil_, ObservType& Observ_,
                       Matrix<dcomplex>& Ham,
                       Matrix<dcomplex>& evec,
                       vector<double>& eval);

// ================================================
// --------- distribute dynamics ------------------
Eigen::VectorXcd ExactDiag(ConstVariables& variables_, QBasis& Basis_,
                           Hamiltonian& Hamil_, ObservType& Observ_) {

	std::cout << "Begin Exact Diagonalization! --- "<<std::endl;
	int N = Basis_.basis.size();
	Matrix<dcomplex> Ham(N,N);
	vector<double> D(N);
	//make the Hamiltonian
	for (int i=0; i<N; i++) {
		dcomplex Hij=Hamil_.HDiag[i]; 
		Ham(i,i)+=Hij;
	}

	int hilbert_t=Hamil_.HTSxx.size();
	// Kxx Kitaev - Sector -------
	for(int i=0;i<hilbert_t;i++){
		int Hi=Hamil_.Sxx_init[i];
		int Hj=Hamil_.Sxx_final[i];
		assert(Hi<N && Hj<N);
		dcomplex Hij=Hamil_.HTSxx[i];
		Ham(Hi,Hj) += Hij;
	}

	hilbert_t=Hamil_.HTSyy.size();
	// Kxx Kitaev - Sector -------
	for(int i=0;i<hilbert_t;i++){
		int Hi=Hamil_.Syy_init[i];
		int Hj=Hamil_.Syy_final[i];
		assert(Hi<N && Hj<N);
		dcomplex Hij=Hamil_.HTSyy[i];
		Ham(Hi,Hj) += Hij;
	}
	assert(Ham.IsHermitian());
    //Ham.print();

    Matrix<dcomplex> Evec = Ham;
	Diagonalize('V',Evec,D);
    for(int i=0;i<D.size();i++) cout << D[i] << ", ";


	cout << " Ground-state energy:= " << setprecision(16) << D[0] << endl;
	//for (int i=0; i<D.size(); i++) cout << D[i] << " \t " << endl;
    return ReturnCol(Evec, 0);
}

void FiniteTemperature(ConstVariables& variables_, QBasis& Basis_,
                       Hamiltonian& Hamil_, ObservType& Observ_, Matrix<dcomplex>& Ham,
                       Matrix<dcomplex>& Evec,
                       vector<double>& eval) {

	cout << " ==================================================================== " << endl;
	cout << "                      ED - Finite Temperature						   " << endl;
	cout << " ==================================================================== " << endl;
	double pi=acos(-1.0);
	int Nsite_ = variables_.NumberofSites;

	double EvTemp_Scale=1; //11605;
	int N = Ham.n_row();
	cout << " energy span: " << eval[N-1] - eval[0] << endl;
	Eigen::VectorXcd Psi = ReturnCol(Evec, 0);

	ofstream myfile;
	myfile.open ("ED_EvsTemp1.dat");

	Eigen::VectorXd Hoverlap(N), Hsqoverlap(N);
	Eigen::VectorXcd jstOlap(N), SzToverlap(N), SxToverlap(N), SyToverlap(N), Skx0(N);
	vector<Matrix<dcomplex>> SiSjN(N);

	for (int j=0; j<N; j++){
		Eigen::VectorXcd jstate = ReturnCol(Evec, j);
		Hoverlap(j) = eval[j];
		Hsqoverlap(j) = eval[j]*eval[j];
		Eigen::VectorXcd SzT = Observ_.ApplySzTotal(jstate,false);
		SzToverlap(j) = jstate.dot(SzT);

		Eigen::VectorXcd SxT = Observ_.ApplySxTotal(jstate,false);
		SxToverlap(j) = jstate.dot(SxT);

		Eigen::VectorXcd SyT = Observ_.ApplySyTotal(jstate,false);
		SyToverlap(j) = jstate.dot(SyT);
		//cout << r_jOlap << " \t " << J_HrOlap << " \t " << J_HrOlap*r_jOlap << endl;

//		SiSjN[j] = SpinSpinCorr(Observ_,jstate, variables_.NumberofSites);
//		Skx0(j) = 0.0;
//		for (int ii=0; ii<Nsite_; ii++){
//			for (int jj=ii; jj<Nsite_; jj++) {
//					int R = ii - jj;
//					dcomplex sum=0;
//					if(ii==jj) sum += SiSjN[j](ii,jj);
//					if(ii!=jj) sum += 2.0*SiSjN[j](ii,jj);
//					Skx0(j) += sum; //SiSjN(j)(ii,jj);  //*cos(0*R) * (1.0/(Nsite_*2));
//			}
//		}
		//exit(1);
	}

	int temax=1000000;
	double dT=0.00005;
	double Temp=dT;
	for(int te=0;te<=temax;te++) {
		//double Temp = Parameters.Temperature/11605;
		//double Temp = 0.000002 + te*0.000002;

		if(Temp>0.2) dT=0.05;
		Temp = Temp + dT;
		double beta = 1.0/Temp;

		double Energy=0,  Partition=0, nsq=0, Ensq=0;
		dcomplex SzT=0, SxT=0, SyT=0, SSk0=0;
		for(int j=0; j<N; j++) {
			double expon = exp(-1.0*beta*(eval[j] -eval[0] ));
			Partition+=expon;
			Energy+= expon*Hoverlap(j);
			Ensq+= expon*Hsqoverlap(j);
			SzT+= expon*SzToverlap(j);
			SxT+= expon*SxToverlap(j);
			SyT+= expon*SyToverlap(j);
			SSk0 += expon*Skx0(j);
//			for(int si=0; si<Nsite_; si++)
//				szkpi+=SiSjN(j)(0,si)
		}
		double Cv = (Ensq/Partition - (Energy/Partition)*(Energy/Partition))/(Temp*Temp);
		double SxavgT = (SxT.real()/Partition)/(Temp);
		double SyavgT = (SyT.real()/Partition)/(Temp);
		double SzavgT = (SzT.real()/Partition)/(Temp);
		double magSus = (SzT.real()/Partition)/(Temp);
		double Sk0 = (SSk0.real()/Partition)/(Temp);

		cout << Temp*EvTemp_Scale << " \t " << Partition
		     << " \t " << Energy
		     << " \t " << Energy/(Partition*Nsite_)
		     << " \t " << Ensq/Partition
		     << " \t " << Cv/Nsite_
		     << " \t " << magSus/Nsite_
		     << " \t " << Sk0/(Nsite_*3)
		     << " \t " << SxavgT/(Partition*Nsite_) << " \t " << SyavgT/(Partition*Nsite_)  << " \t " << SzavgT/(Partition*Nsite_)
		     << endl;

		myfile << Temp*EvTemp_Scale << " \t " << Partition
		       << " \t " << Energy
		       << " \t " << Energy/(Partition*Nsite_)
		       << " \t " << Ensq/Partition
		       << " \t " << Cv/Nsite_
		       << " \t " << magSus/Nsite_
		       << " \t " << Sk0/(Nsite_*3)
		       << " \t " << SxavgT/(Partition*Nsite_) << " \t " << SyavgT/(Partition*Nsite_)  << " \t " << SzavgT/(Partition*Nsite_)
		       << endl;

		if(Temp*EvTemp_Scale>1000) break;
	}
	myfile.close();




}


void Diagonalize(char option, Matrix<dcomplex>& A, vector<double>& eval){
	char jobz=option;
	char uplo='U';
	int n=A.n_row();
	int lda=A.n_col();
	assert(n==lda);

	eval.resize(n);
	vector<dcomplex> work(3);
	vector<double> rwork(3*n);
	int info,lwork= -1;

	// ---- spin up part of the Hamiltonian --------------------
	fill(eval.begin(),eval.end(),0);
	// query:
	zheev_(&jobz,&uplo,&n,&(A(0,0)),&lda,&(eval[0]),&(work[0]),&lwork,&(rwork[0]),&info);
	lwork = int(real(work[0]))+1;
	work.resize(lwork+1);
	// real work:
	zheev_(&jobz,&uplo,&n,&(A(0,0)),&lda,&(eval[0]),&(work[0]),&lwork,&(rwork[0]),&info);
	if (info!=0) {
		std::cerr<<"info="<<info<<"\n";
		perror("diag: zheev: failed with info!=0.\n");
	}
	//for(int i=0;i<n;i++) { cout << eval[i] << " \t ";} cout << endl;
	//sort(eigs_.begin(),eigs_.end()); // sort Eigenvalues and Hamiltonian
}




void Diagonalize(char option, Eigen::MatrixXcd& A, Eigen::VectorXd& eval){
    char jobz=option;
    char uplo='U';
    int n=A.rows();
    int lda=A.cols();
    assert(n==lda);

    eval.resize(n);
    eval.fill(0);
    Eigen::VectorXcd work(3);
    Eigen::VectorXd rwork(3*n);
    int info,lwork= -1;

    // ---- spin up part of the Hamiltonian --------------------
    // query:
    zheev_(&jobz,&uplo,&n,&(A(0,0)),&lda,&(eval[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    lwork = int(real(work[0]))+1;
    work.resize(lwork+1);
    // real work:
    zheev_(&jobz,&uplo,&n,&(A(0,0)),&lda,&(eval[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }
    //for(int i=0;i<n;i++) { cout << eval[i] << " \t ";} cout << endl;
    //sort(eigs_.begin(),eigs_.end()); // sort Eigenvalues and Hamiltonian
}

// ================================================
// --- return ith row of the matrix;

Eigen::VectorXcd ReturnRow(Matrix<dcomplex>& mat, int row) {
	int N = mat.n_col();
	Eigen::VectorXcd out(N); out.setZero();
	for(int i=0; i<N;i++) out(i) = mat(row,i);
	//VectorXcd out1=out;
	return out;
}

VectorXf ReturnRow(Matrix<double>& mat, int row) {
	int N = mat.n_col();
	VectorXf out(N); out.setZero();
	for(int i=0; i<N;i++) out[i] = mat(row,i);
	return out;
}

Eigen::VectorXcd ReturnCol(Matrix<dcomplex>& mat, int col) {
	int N = mat.n_col();
	Eigen::VectorXcd out(N); out.setZero();
	for(int i=0; i<N;i++) out(i) = mat(i,col);
	return out;
}

VectorXf ReturnCol(Matrix<double>& mat, int col) {
	int N = mat.n_col();
	VectorXf out(N); out.setZero();
	for(int i=0; i<N;i++) out[i] = mat(i,col);
	return out;
}




#endif // EXACTDIAG_H

