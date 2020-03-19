#ifndef GEOMETRYHONEYCOMB_H
#define GEOMETRYHONEYCOMB_H

class Geometry {
public:
	Hamiltonian(ConstVariables& variables, QBasis& Basis):
	    variables_(variables), Basis_(Basis),
	    Nsite_(variables.NumberofSites)
	{ }

	std::vector<int> indx_,indy_;
	Matrix<int> Nc_;
	Matrix<int> Kxx_pair, Kyy_pair, Kzz_pair;
	std::vector<double> Kxx_Conn,Kyy_Conn,Kzz_Conn;
	std::vector<int> Sxx_init, Sxx_final;
	std::vector<int> Syy_init, Syy_final;
	std::vector<int> Szz_init, Szz_final;
	std::vector<dcomplex> HTSxx,HTSyy,HDiag;

	// -------------------------------------------------------
	void TightB_Ham(){
		if (variables_.Model=="HeisenbergChain") {
			HeisenbergChain_Connectors();
		} else if (variables_.Model=="Kitaev") {
			Connectors(); // Build Connectors in Single-particle states
		} else {
			std::cerr<<"Model="<<variables_.Model<<"\n";
	        throw std::string("Unknown Model parameter \n");
		}
		std::cout << "create Tight-Binding Hamiltonian" << std::endl;
		std::cout << "create Off-Diag matrix elements" << std::endl;
		int Hsize = Basis_.basis.size();
		std::cout << "	hilbert-space size:  "<< Hsize << std::endl;

		// allocate diagonal term!
		dcomplex izero=(0.0,0.0);
		std::fill(HTSxx.begin(),HTSxx.end(),izero);
		std::fill(HTSyy.begin(),HTSyy.end(),izero);

		// Sxi Sx_j |Psi> = distrop down and create up
		// (create down remove up)_i (create up remove down)_j
		for(int bond=0;bond<Kxx_Conn.size();bond++) {
			int si=Kxx_pair(bond,0);
			int sj=Kxx_pair(bond,1);
			if(Kxx_Conn.size()==1) break;
			assert(si!=sj);
			for(int i1=0; i1<Hsize; i1++) {
				int ket = Basis_.basis[i1];

				int keto1,keto2;
				std::complex<double> Sxi,Sxj;
				Basis_.ApplySx(si,ket,keto1,Sxi);
				Basis_.ApplySx(sj,keto1,keto2,Sxj);

				//cout << Sxi << " \t " << Sxj << " \t " << Sxi*Sxj << endl;
				Sxx_init.push_back(ket);
				Sxx_final.push_back(keto2);
				HTSxx.push_back(Sxi*Sxj*Kxx_Conn[bond]);
			}
		}

		// Syi Sy_j |Psi> = distrop down and create up
		// (create down remove up)_i (create up remove down)_j
		for(int bond=0;bond<Kyy_Conn.size();bond++) {
			int si=Kyy_pair(bond,0);
			int sj=Kyy_pair(bond,1);
			if(Kyy_Conn.size()==1) break;
			assert(si!=sj);
			for(int i1=0; i1<Hsize; i1++) {
				int ket = Basis_.basis[i1];

				int keto1,keto2;
				std::complex<double> Syi,Syj;
				Basis_.ApplySy(si,ket,keto1,Syi);
				Basis_.ApplySy(sj,keto1,keto2,Syj);

				//cout << Syi << " \t " << Syj << " \t " << Syi*Syj << endl;
				Syy_init.push_back(ket);
				Syy_final.push_back(keto2);
				HTSyy.push_back(Syi*Syj*Kyy_Conn[bond]);
			}
		}



		// --- ADD Magnetic Field in spin-x and spin-y direction
		// x-direction
		for(int i1=0; i1<Hsize; i1++) {
			int ket = Basis_.basis[i1];
			for(int site=0;site<Nsite_;site++) {
				int keto1=-1;
				//int bit = int(findBit(ketin,site));
				keto1 = Basis_.flipBit(ket,site);
				dcomplex coef=0.5+0i;
				Sxx_init.push_back(keto1);
				Sxx_final.push_back(ket);
				HTSxx.push_back(coef*variables_.Bxx);
			}
		}

		// y-direction
		for(int i1=0; i1<Hsize; i1++) {
			int ket = Basis_.basis[i1];
			for(int site=0;site<Nsite_;site++) {
				int keto1=-1;
				dcomplex coef=0.+0.i;
				Basis_.ApplySy(site,ket,keto1,coef);
				dcomplex magy = variables_.Byy*coef;
				Syy_init.push_back(keto1);
				Syy_final.push_back(ket);
				HTSyy.push_back(magy);
			}
		}

	} // end function

	// -------------------------------------------------------
	double findsign(int newstate,int hopp_Fr,int hopp_To){
		int particles=0;

		if(abs(hopp_Fr-hopp_To)==0 || abs(hopp_Fr-hopp_To)==1) {
			return 1.0;
		}

		if(hopp_To<hopp_Fr){
			for(int pos=hopp_To+1;pos<hopp_Fr;pos++){
				if(Basis_.findBit(newstate,pos)==1) particles++;
			}
		}
		else {
			for(int pos=hopp_Fr+1;pos<hopp_To;pos++){
				if(Basis_.findBit(newstate,pos)==1) particles++;
			}
		}
		return pow(-1.0,particles);
	} // end function

	// -------------------------------------------------------
	void Connectors(){
		std::cout << "creating lattice and Tight-Binding Connectors:" << std::endl;

		int LLX = variables_.Lx;
		int LLY = variables_.Ly;
		int SqLLX = LLX*2;
		int SqLLY = LLY;
		int SqNsite_=SqLLX*SqLLY;

		Matrix<int> SqNc(SqLLX,SqLLY);
		vector<int> Sqindx(SqNsite_),Sqindy(SqNsite_);
		SqNc.fill(-1);

		int counter=0;
		for(int ix=0;ix<SqLLX;ix++){
			for(int iy=0;iy<SqLLY;iy++){
				if(SqNc(ix,iy)==-1){
					if(ix%2==1 and iy%4==0) {
						SqNc(ix,iy) = 1;
						Sqindx[counter] = ix;
						Sqindy[counter] = iy;
						counter+=1;
						//print ix, iy

						if(iy+3<SqLLY){
							SqNc(ix,iy+3) = 1;
							Sqindx[counter] = ix;
							Sqindy[counter] = iy+3;
							counter+=1;
							//print ix, iy+3
						}
					}

					if(ix%2==0 and iy%4==0) {

						if(iy+1<SqLLY){
							SqNc(ix,iy+1) = 1;
							Sqindx[counter] = ix;
							Sqindy[counter] = iy+1;
							counter+=1;
							//print ix, iy+1
						}

						if(iy+2<SqLLY){
							SqNc(ix,iy+2) = 1;
							Sqindx[counter] = ix;
							Sqindy[counter] = iy+2;
							counter+=1;
							//print ix, iy+2
						}
					}

				}
			}
		}


		// Site labeling
		indx_.resize(Nsite_);   indy_.resize(Nsite_);
		Nc_.resize(SqLLX,SqLLY);


		counter = 0;
		for(int ix=0;ix<SqLLX;ix++){
			for(int iy=0;iy<SqLLY;iy++){
				Nc_(ix,iy) = -1;
				if(SqNc(ix,iy)!=-1) {
					Nc_(ix,iy) = counter;
					indx_[counter]=ix;
					indy_[counter]=iy;
					//print ix, iy, counter
					counter+=1;
				}
			}
		}


		Matrix<double> SinglePart_Kxx, SinglePart_Kyy; // local
		Matrix<double> SinglePart_Kzz;
		SinglePart_Kxx.resize(Nsite_,Nsite_);
		SinglePart_Kyy.resize(Nsite_,Nsite_);
		SinglePart_Kzz.resize(Nsite_,Nsite_);

		for(int i=0;i<Nsite_;i++){ 	// ith site
			int ix = indx_[i];
			int iy = indy_[i];

			if(iy+1<SqLLY)  {
				int jy = iy+1;
				int jx = ix;
				if(Nc_(jx,jy)!=-1){
					int j = Nc_(jx,jy);
					SinglePart_Kzz(i,j) = variables_.Kzz;
					//SinglePart_Kzz(j,i) = variables_.Kzz;
				}
			}


			if(ix+1<SqLLX && iy-1>-1) {
				int jy = iy-1;
				int jx = ix+1;
				if(Nc_(jx,jy)!=-1) {
					int j = Nc_(jx,jy);
					SinglePart_Kxx(i,j) = variables_.Kxx;
					//SinglePart_Kxx(j,i) = variables_.Kxx;
				}
			}



			if(ix-1>-1 && iy-1>-1) {
				int jy = iy-1;
				int jx = ix-1;
				if(Nc_(jx,jy)!=-1) {
					int j = Nc_(jx,jy);
					SinglePart_Kyy(i,j) = variables_.Kyy;
					//SinglePart_Kyy(j,i) = variables_.Kyy;
				}
			}

			// - Ly periodic -
			if(iy==LLY-1 && variables_.IsPeriodicY==true) {
				int jy = 0;
				int jx = ix;
				if(Nc_(jx,jy)!=-1){
					int j = Nc_(jx,jy);
					SinglePart_Kzz(i,j) = variables_.Kzz;
					//SinglePart_Kzz(j,i) = variables_.Kzz;
				} else if (ix-1>-1) {
					jx = ix-1;
					int j = Nc_(jx,jy);
					SinglePart_Kzz(i,j) = variables_.Kzz;
					//SinglePart_Kzz(j,i) = variables_.Kzz;
				}
			}

			// - Lx periodic -
			if(variables_.IsPeriodicX==true) {

				// - Sx bond -
				if(iy-1>-1 && ix==SqLLX-1) {
					int jy = iy-1;
					int jx = 0;
					if(Nc_(jx,jy)!=-1) {
						int j = Nc_(jx,jy);
						SinglePart_Kxx(i,j) = variables_.Kxx;
						//SinglePart_Kxx(j,i) = variables_.Kxx;
					}
				}

				// - Sy bond -
				if(ix==0 && iy-1>-1) {
					int jy = iy-1;
					int jx = SqLLX-1;
					if(Nc_(jx,jy)!=-1) {
						int j = Nc_(jx,jy);
						SinglePart_Kyy(i,j) = variables_.Kyy;
						//SinglePart_Kyy(j,i) = variables_.Kyy;
					}
				}


			}





		}


		int Kxxbonds=SinglePart_Kxx.numNonZeros();
		int Kyybonds=SinglePart_Kyy.numNonZeros();
		int Kzzbonds=SinglePart_Kzz.numNonZeros();

		cout << " Single particle Kitaev - xx " << endl;
		SinglePart_Kxx.print();

		cout << " Single particle Kitaev - yy " << endl;
		SinglePart_Kyy.print();

		cout << " Single particle Kitaev - zz " << endl;
		SinglePart_Kzz.print();


		if(Kxxbonds==0) Kxxbonds=1;
		if(Kyybonds==0) Kyybonds=1;
		if(Kzzbonds==0) Kzzbonds=1;
		Kxx_pair.resize(Kxxbonds,2); Kyy_pair.resize(Kyybonds,2); Kzz_pair.resize(Kzzbonds,2);
		Kxx_Conn.resize(Kxxbonds);   Kyy_Conn.resize(Kyybonds);   Kzz_Conn.resize(Kzzbonds);

		counter=0;
		for(int i=0;i<Nsite_;i++) {
			for(int j=0;j<Nsite_;j++) {
				if(SinglePart_Kxx(i,j)!=0.0) {
					//assert(j>i);
					Kxx_pair(counter,0) = i;
					Kxx_pair(counter,1) = j;
					Kxx_Conn[counter] = SinglePart_Kxx(i,j);
					counter++;
				}
			}
		}

		counter=0;
		for(int i=0;i<Nsite_;i++) {
			for(int j=0;j<Nsite_;j++) {
				if(SinglePart_Kyy(i,j)!=0.0) {
					//assert(j>i);
					Kyy_pair(counter,0) = i;
					Kyy_pair(counter,1) = j;
					Kyy_Conn[counter] = SinglePart_Kyy(i,j);
					counter++;
				}
			}
		}

		counter=0;
		for(int i=0;i<Nsite_;i++) {
			for(int j=0;j<Nsite_;j++) {
				if(SinglePart_Kzz(i,j)!=0.0) {
					//assert(j>i);
					Kzz_pair(counter,0) = i;
					Kzz_pair(counter,1) = j;
					Kzz_Conn[counter] = SinglePart_Kzz(i,j);
					counter++;
				}
			}
		}



	} // end function


	// -------------------------------------------------------
	void Diagonal(){
		std::cout << "create Diagonal matrix elements" << std::endl;
		int Hsize = Basis_.basis.size();
		std::cout << "	hilbert-space size:  "<< Hsize << std::endl;

		// allocate diagonal term!
		HDiag.resize(Hsize);
		dcomplex izero=(0.0,0.0);
		std::fill(HDiag.begin(),HDiag.end(),izero);

		// Szi Sz_j |Psi> = distrop down and create up
		// (create down remove up)_i (create up remove down)_j
		for(int bond=0;bond<Kzz_Conn.size();bond++) {
			int si=Kzz_pair(bond,0);
			int sj=Kzz_pair(bond,1);
			if(Kzz_Conn.size()==1) break;
			assert(si!=sj);
			for(int i1=0; i1<Hsize; i1++) {
					int ket = Basis_.basis[i1];
					int keto1,keto2;
					dcomplex Szi,Szj;

					Basis_.ApplySz(si,ket,keto1,Szi); assert(ket==keto1);
					Basis_.ApplySz(sj,keto1,keto2,Szj); assert(keto1==keto2);

					HDiag[i1] += Kzz_Conn[bond]*Szi*Szj;
			}
		}


		// --- ADD Magnetic Field in spin-z-direction
		for(int i1=0; i1<Hsize; i1++) {
			int ket = Basis_.basis[i1];
			dcomplex SzT=0;

			for(int site=0;site<Nsite_;site++) {
				dcomplex Szi; int keto1;
				Basis_.ApplySz(site,ket,keto1,Szi); assert(ket==keto1);
				SzT+=Szi;
			}
			HDiag[i1] += variables_.Bzz*SzT;
		}

	} // end function



	// -------------------------------------------------------
	void HeisenbergChain_Connectors(){
		std::cout << "creating lattice and Tight-Binding Connectors:" << std::endl;

		int LLX,LLY;

		// ** FIXME: Testing on chain first **
		LLX = Nsite_; //variables_.Lx;
		LLY = 1; //variables_.Ly;

		indx_.resize(Nsite_);   indy_.resize(Nsite_);
		Nc_.resize(LLX,LLY);
		Matrix<double> SinglePart_Kxx, SinglePart_Kyy; // local
		Matrix<double> SinglePart_Kzz;
		SinglePart_Kxx.resize(Nsite_,Nsite_);
		SinglePart_Kyy.resize(Nsite_,Nsite_);
		SinglePart_Kzz.resize(Nsite_,Nsite_);

		// Site labeling
		int site=0;
		for(int ix=0;ix<LLX;ix++){
			for(int iy=0;iy<LLY;iy++){
				indy_[site]=iy;
				indx_[site]=ix;
				Nc_(ix,iy)=site;
				site++;
			}
		}


		// Neighbors for each site
		int bonds=0, Kxxbonds=0, Kyybonds=0, Kzzbonds=0;
		for(site=0;site<Nsite_;site++){ 	// ith site
			int lx=indx_[site];
			int ly=indy_[site];

			if(lx<LLX-1) {
				int mx = lx + 1; // +x hopping
				int my = ly;
				int m = Nc_(mx,my);

				SinglePart_Kxx(site,m) = variables_.Kxx;
				if(variables_.Kxx!=0.0) Kxxbonds++;

				SinglePart_Kyy(site,m) = variables_.Kyy;
				if(variables_.Kyy!=0.0) Kyybonds++;

				SinglePart_Kzz(site,m) = variables_.Kzz;
				if(variables_.Kzz!=0.0) Kzzbonds++;
			}


			// periodic boundary conditions in x direction
			if(lx==LLX-1 && variables_.IsPeriodicX==true) {
				int mx = 0; // +x hopping
				int my = ly;
				int m = Nc_(mx,my);

				SinglePart_Kxx(m,site) = variables_.Kxx;
				if(variables_.Kxx!=0.0) Kxxbonds++;

				SinglePart_Kyy(m,site) = variables_.Kyy;
				if(variables_.Kyy!=0.0) Kyybonds++;

				SinglePart_Kzz(m,site) = variables_.Kzz;
				if(variables_.Kzz!=0.0) Kzzbonds++;
			}
		}


		cout << " Single particle Kitaev - xx " << endl;
		SinglePart_Kxx.print();

		cout << " Single particle Kitaev - yy " << endl;
		SinglePart_Kyy.print();

		cout << " Single particle Kitaev - zz " << endl;
		SinglePart_Kzz.print();


		if(Kxxbonds==0) Kxxbonds=1;
		if(Kyybonds==0) Kyybonds=1;
		if(Kzzbonds==0) Kzzbonds=1;
		Kxx_pair.resize(Kxxbonds,2); Kyy_pair.resize(Kyybonds,2); Kzz_pair.resize(Kzzbonds,2);
		Kxx_Conn.resize(Kxxbonds);   Kyy_Conn.resize(Kyybonds);   Kzz_Conn.resize(Kzzbonds);

		int counter=0;
		for(int i=0;i<Nsite_;i++) {
			for(int j=0;j<Nsite_;j++) {
				if(SinglePart_Kxx(i,j)!=0.0) {
					Kxx_pair(counter,0) = i;
					Kxx_pair(counter,1) = j;
					Kxx_Conn[counter] = SinglePart_Kxx(i,j);
					counter++;
				}
			}
		}

		counter=0;
		for(int i=0;i<Nsite_;i++) {
			for(int j=0;j<Nsite_;j++) {
				if(SinglePart_Kyy(i,j)!=0.0) {
					Kyy_pair(counter,0) = i;
					Kyy_pair(counter,1) = j;
					Kyy_Conn[counter] = SinglePart_Kyy(i,j);
					counter++;
				}
			}
		}

		counter=0;
		for(int i=0;i<Nsite_;i++) {
			for(int j=0;j<Nsite_;j++) {
				if(SinglePart_Kzz(i,j)!=0.0) {
					Kzz_pair(counter,0) = i;
					Kzz_pair(counter,1) = j;
					Kzz_Conn[counter] = SinglePart_Kzz(i,j);
					counter++;
				}
			}
		}



	} // end function







	// -------------------------------------------------------
	// Matrix Vector Multiplication H|vin> = |vout>
	Vector<dcomplex> MatrixVectorMult(Vector<dcomplex>& Vin) {
		if(Vin.size()!=HDiag.size()) {
			string errorout ="LanczosSolver.h::Dot-Product size mismatch!";
			throw std::invalid_argument(errorout);
		}

		int hilbert=Vin.size();
		Vector<dcomplex> Vout(hilbert);
		Vout.clear();
		//Vout = VectorXcd::Zero(hilbert); // "Just incase" - set all elements to zero

		// Diagonal Hamiltonian terms!
		for(int i=0;i<hilbert;i++) {
			dcomplex Hij=HDiag[i];
			//assert(Hij.imag()==0);
			Vout(i) += Hij*Vin(i);
		}

		int hilbert_t=HTSxx.size();
		// Kxx Kitaev - Sector -------
		for(int i=0;i<hilbert_t;i++){
			int Hi=Sxx_init[i];
			int Hj=Sxx_final[i];
			dcomplex Hij=HTSxx[i];
			assert(Hi<hilbert && Hj<hilbert);
			//assert(Hij.imag()==0);
			Vout(Hi) += Hij*Vin(Hj);
//			Vout[Hi] += Hij.real()*Vin[Hj];
		}

		hilbert_t=HTSyy.size();
		// Kyy Kitaev - Sector -------
		for(int i=0;i<hilbert_t;i++){
			int Hi=Syy_init[i];
			int Hj=Syy_final[i];
			dcomplex Hij=HTSyy[i];
			assert(Hi<hilbert && Hj<hilbert);
			//assert(Hij.imag()==0);
			Vout(Hi) += Hij*Vin(Hj);
			//cout << Hij << " \t " << Vin(Hj) << " \t " << Hij*Vin(Hj) << endl;
		}

		return Vout;
	}

	// -------------------------------------------------------
	// Matrix Vector Multiplication H|vin> = |vout>
	Vector<double> MatrixVectorMult(Vector<double>& Vin) {
		if(Vin.size()!=HDiag.size()) {
			string errorout ="LanczosSolver.h::Dot-Product size mismatch!";
			throw std::invalid_argument(errorout);
		}

		int hilbert=Vin.size();
		Vector<double> Vout(hilbert);
		Vout.clear();
		//Vout = VectorXcd::Zero(hilbert); // "Just incase" - set all elements to zero

		// Diagonal Hamiltonian terms!
		for(int i=0;i<hilbert;i++) {
			dcomplex Hij=HDiag[i];
			assert(Hij.imag()==0);
			Vout(i) += Hij.real()*Vin(i);
		}

		int hilbert_t=HTSxx.size();
		// Kxx Kitaev - Sector -------
		for(int i=0;i<hilbert_t;i++){
			int Hi=Sxx_init[i];
			int Hj=Sxx_final[i];
			dcomplex Hij=HTSxx[i];
			assert(Hi<hilbert && Hj<hilbert);
			assert(Hij.imag()==0);
			Vout(Hi) += Hij.real()*Vin(Hj);
//			Vout[Hi] += Hij.real()*Vin[Hj];
		}

		hilbert_t=HTSyy.size();
		// Kyy Kitaev - Sector -------
		for(int i=0;i<hilbert_t;i++){
			int Hi=Syy_init[i];
			int Hj=Syy_final[i];
			dcomplex Hij=HTSyy[i];
			assert(Hi<hilbert && Hj<hilbert);
			assert(Hij.imag()==0);
			Vout(Hi) += Hij.real()*Vin(Hj);
			//cout << Hij << " \t " << Vin(Hj) << " \t " << Hij*Vin(Hj) << endl;
		}

		return Vout;
	}





private:
	ConstVariables& variables_;
	QBasis& Basis_;
	std::vector<double> Conn, TJH_Conn, Pair_Conn, TJHSz_Conn;
	int Nsite_,Nup_,Ndn_;
};


#endif // GEOMETRYHONEYCOMB_H

