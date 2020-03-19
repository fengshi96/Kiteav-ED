#ifndef OBSERVABLES_H
#define OBSERVABLES_H

typedef std::complex<double> dcomplex;
typedef Eigen::VectorXf VectorXf;

extern "C" void   dstev_(char *,int *, double *, double *, double *, int *, double *, int *);
class Observables {
public:
    Observables(ConstVariables& variables, Lattice& Lat, QBasis& Basis, Hamiltonian& Hamil):
        variables_(variables), Lat_(Lat), Basis_(Basis), Hamil_(Hamil)
    {

    }

    // -----------------------------------
    Eigen::VectorXcd ApplyOp(Eigen::VectorXcd& inV, const string OpA, const int& site) {
        int Nsite_ = variables_.NumberofSites;
        int Hsize = Basis_.basis.size();
        Eigen::VectorXcd out;

        if (OpA == "Sz") {
           return ApplySz(inV, site);
        } else if (OpA == "Sy") {
           return ApplySy(inV, site);
        } else if (OpA == "Sx") {
            return ApplySx(inV, site);
        } else if (OpA == "I") {
            return inV;
        } else {
            std::cerr<<"Operator="<<OpA<<" does not exist \n";
            throw std::string("Unknown Operator label \n");
        }
    }

    // -----------------------------------
    void measureLocalS(Eigen::VectorXcd& GS_, Eigen::VectorXcd& Sxi,
                       Eigen::VectorXcd& Syi, Eigen::VectorXcd& Szi,
                       Eigen::VectorXcd& Spi, Eigen::VectorXcd& Smi) {
        int Nsite_ = variables_.NumberofSites;
        Sxi.resize(Nsite_); Syi.resize(Nsite_); Szi.resize(Nsite_);
        Spi.resize(Nsite_); Smi.resize(Nsite_);

        cout << "\n-------- measuring <sx>, <sy>, <sz>, <sp>, <sm> --------" << endl;
        for(int si=0; si<Nsite_; si++) {
            Eigen::VectorXcd tmpx = ApplySx(GS_,si);
            Eigen::VectorXcd tmpy = ApplySy(GS_,si);
            Eigen::VectorXcd tmpz = ApplySz(GS_,si);
            Sxi[si] = GS_.dot(tmpx);
            Syi[si] = GS_.dot(tmpy);
            Szi[si] = GS_.dot(tmpz);
            Spi[si] = Sxi[si] + dcomplex(0.0,1.0)*Syi[si];
            Smi[si] = Sxi[si] - dcomplex(0.0,1.0)*Syi[si];
            cout << si << " \t " << Sxi[si] << " \t " << Syi[si] << " \t " << Szi[si]
                 << " \t " << Spi[si] << " \t " << Smi[si] << endl;
        }
    }

    // -----------------------------------
    Eigen::MatrixXcd TwoPointCorr(const string OpA, const string OpB, Eigen::VectorXcd& GS_) {
        //cout << " -------- " << endl;
        std::cout << "Calculating " << OpA << "." << OpB << " correlations " << std::endl;
        int Nsite_ = variables_.NumberofSites;
        int Hsize = Basis_.basis.size();
        Eigen::MatrixXcd Corr(Nsite_,Nsite_);
        Corr.setZero();
        dcomplex dczero(0.0,0.0);

        // #pragma omp parallel for
        for(int sj=0; sj<Nsite_; sj++) {
            Eigen::VectorXcd out1(Hsize);
            out1.setZero();
            out1 = ApplyOp(GS_, OpB, sj);
            for(int si=0; si<Nsite_; si++) {
                if(sj>=si) {
                    Eigen::VectorXcd out2(Hsize);
                    out2.setZero();
                    out2 = ApplyOp(out1, OpA, si);
                    Corr(si,sj) = GS_.dot(out2);
                    Corr(sj,si) = Corr(si,sj);
                } else {
                    Corr(si,sj) = dczero;
                }
            }
        }
        //Corr.print();
        return Corr;
    }

    // -----------------------------------
    Eigen::VectorXcd ApplySz(Eigen::VectorXcd& GS_, const int& i) {
        int Nsite_ = variables_.NumberofSites;
        assert(i<Nsite_);
        int Hsize = Basis_.basis.size();
        Eigen::VectorXcd out(Hsize); out.setZero();

        for(int i1=0; i1<Hsize; i1++) {
            int ket = Basis_.basis[i1];
            int keto1;
            dcomplex Szi;

            Basis_.ApplySz(i,ket,keto1,Szi);
            assert(ket==keto1);
            out(keto1) += Szi*GS_(ket);
        }
        return out;
    }

    // -----------------------------------
    Eigen::VectorXcd ApplySy(Eigen::VectorXcd& GS_, const int& i) {
        int Nsite_ = variables_.NumberofSites;
        assert(i<Nsite_);
        int Hsize = Basis_.basis.size();
        Eigen::VectorXcd out(Hsize); out.setZero();

        for(int i1=0; i1<Hsize; i1++) {
            int ket = Basis_.basis[i1];
            int keto1;
            dcomplex Szi;

            Basis_.ApplySy(i,ket,keto1,Szi);
            out(keto1) += Szi*GS_(ket);
        }
        return out;
    }

    // -----------------------------------
    Eigen::VectorXcd ApplySx(Eigen::VectorXcd& GS_, const int& i) {
        int Nsite_ = variables_.NumberofSites;
        assert(i<Nsite_);
        int Hsize = Basis_.basis.size();
        Eigen::VectorXcd out(Hsize); out.setZero();

        for(int i1=0; i1<Hsize; i1++) {
            int ket = Basis_.basis[i1];
            int keto1;
            dcomplex Szi;

            Basis_.ApplySx(i,ket,keto1,Szi);
            out(keto1) += Szi*GS_(ket);
        }
        return out;
    }

    // -----------------------------------
    Eigen::VectorXcd ApplySzTotal(Eigen::VectorXcd& GS_, bool square) {
        int Nsite_ = variables_.NumberofSites;
        int Hsize = Basis_.basis.size();

        Eigen::VectorXcd out(Hsize);
        for(int i1=0; i1<Hsize; i1++) {
            int ket = Basis_.basis[i1];
            int keto1;
            dcomplex Szt=0;

            for (int s=0; s<Nsite_;s++) {
                dcomplex Szi;
                Basis_.ApplySz(s,ket,keto1,Szi);
                assert(ket==keto1);
                assert(Szi.imag()==0);
                Szt+=Szi;
            }
            if(square) {
                out(ket) = Szt*Szt*GS_(ket);
            } else {
                out(ket) = Szt*GS_(ket);
            }

        }
        return out;
    }

    // -----------------------------------
    Eigen::VectorXcd ApplySxTotal(Eigen::VectorXcd& GS_, bool square) {
        int Nsite_ = variables_.NumberofSites;
        int Hsize = Basis_.basis.size();
        Eigen::VectorXcd out(Hsize);
        out.setZero();

        for (int s=0; s<Nsite_;s++) {
            out += ApplySx(GS_,s);
        }
        return out;
    }

    // -----------------------------------
    Eigen::VectorXcd ApplySyTotal(Eigen::VectorXcd& GS_, bool square) {
        int Nsite_ = variables_.NumberofSites;
        int Hsize = Basis_.basis.size();
        Eigen::VectorXcd out(Hsize);
        out.setZero();

        for (int s=0; s<Nsite_;s++) {
            out += ApplySy(GS_,s);
        }
        return out;
    }

private:
    ConstVariables& variables_;
    Lattice& Lat_;
    QBasis& Basis_;
    Hamiltonian& Hamil_;
    //VectorXf& GS_;
};
#endif // OBSERVABLES_H
