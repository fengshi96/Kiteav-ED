/*  
 * Author: Nirav D. Patel
 * file: ParametersEngine.h
 * Contains the constant variables
*/
#include <exception>

typedef std::mt19937_64 RDMGeneratortype;

struct ConstVariables {
	ConstVariables(string inputfile)
	    : NumberofSites(0)
	{
		Initialize(inputfile);
	}

	RDMGeneratortype generator_;
	int Lx, Ly, NumberofSites,SpecialSite;
	int LanczosSteps,RandomSeed, Threads;
	std::uniform_real_distribution<double> dis_;
	double LanczosEps;
	double Kxx,Kyy,Kzz,Bxx,Byy,Bzz,Temperature;
	bool IsPeriodicX,IsPeriodicY;
    bool IsAPeriodicX,IsAPeriodicY;
    string Model,Solver,LancType;

    bool LocalSkw,RIXSDynamics,SpinCurrDynamics;
	int RIXSIntermediate_Steps;

	void Initialize(string inputfile);
	double grnd(); // double random number generator
	double matchstring(string file,string match);
	string matchstring2(const string file,const string match);
};



/*
 * ***********
 *  Functions in Class Coordinates ------
 *  ***********
*/
void ConstVariables::Initialize(string inputfile) {

	cout << "____________________________________" << endl;
	cout << " - Reading the inputfile: " << inputfile << endl;
	cout << "____________________________________" << endl;
	Lx = int(matchstring(inputfile,"NumberofSitesLx"));
	Ly = int(matchstring(inputfile,"NumberofSitesLy"));
	NumberofSites = Lx*Ly;

	RandomSeed = matchstring(inputfile,"RandomSeed");
	RDMGeneratortype generator(RandomSeed);
	generator_=generator;

	Model = matchstring2(inputfile,"Model");
	Solver = matchstring2(inputfile,"Solver");
	Kxx = matchstring(inputfile,"Kxx");
	Kyy = matchstring(inputfile,"Kyy");
	Kzz = matchstring(inputfile,"Kzz");
	Bxx = matchstring(inputfile,"Bxx");
	Byy = matchstring(inputfile,"Byy");
	Bzz = matchstring(inputfile,"Bzz");
	Temperature = matchstring(inputfile,"Temperature");
	LanczosEps = matchstring(inputfile,"LanczosEps");
	IsPeriodicX = bool(matchstring(inputfile,"IsPeriodicX"));
    IsPeriodicY = bool(matchstring(inputfile,"IsPeriodicY"));
    IsAPeriodicX = bool(matchstring(inputfile,"IsAPeriodicX"));
    IsAPeriodicY = bool(matchstring(inputfile,"IsAPeriodicY"));
    LancType = matchstring2(inputfile,"LancType");
    if(LancType=="OnlyEgs") cout << "WARNING:: Reortho will be turned off for LancType=OnlyEgs (set by default) " << endl;
	Threads = int(matchstring(inputfile,"Threads"));



	// -- Lanczos Steps -- by default set to 200
	try {
		LanczosSteps = matchstring(inputfile,"LanczosSteps");
	} catch (exception& e) {
		LanczosSteps = 200;
		cout << "LanczosSteps = " << LanczosSteps <<  " (set by default) " << endl;
	}

	try {
		SpecialSite = matchstring(inputfile,"SpecialSite");
	} catch (exception& e) {
		SpecialSite = NumberofSites/2;
		cout << "SpecialSite = " << SpecialSite <<  " (set by default) " << endl;
    }


	// ========== Dynamical correlations =====================
    try { // --- RIXS ----
		RIXSDynamics = matchstring(inputfile,"RIXSDynamics");
		RIXSIntermediate_Steps = matchstring(inputfile,"RIXSIntermediate_Steps");
	} catch (exception& e) {
		RIXSDynamics = false;
		RIXSIntermediate_Steps = 0;
	}
    if(RIXSDynamics==true && LancType!="All") {
        string errorout ="You need to set LancType=All for Dynamical calculations --- you dummy!";
        throw std::invalid_argument(errorout);
    }

    try { // --- LocalSkw ----
		LocalSkw = matchstring(inputfile,"LocalSkw");
	} catch (exception& e) {
		LocalSkw = false;
	}
    if(LocalSkw==true && LancType!="All") {
        string errorout ="You need to set LancType=All for Dynamical calculations --- you dummy!";
        throw std::invalid_argument(errorout);
    }

    try { // --- Spin Current ----
        SpinCurrDynamics = matchstring(inputfile,"SpinCurrDynamics");
    } catch (exception& e) {
        SpinCurrDynamics = false;
    }
    if(SpinCurrDynamics==true && LancType!="All") {
        string errorout ="You need to set LancType=All for Dynamical calculations --- you dummy!";
        throw std::invalid_argument(errorout);
    }





	cout << "____________________________________" << endl << endl;
}


double ConstVariables::grnd(){
    return dis_(generator_);
} // ----------


double ConstVariables::matchstring(string file, const string match) {
	string test;
	string line;
	ifstream readFile(file);
	double amount;
	bool pass=false;
	while (std::getline(readFile, line)) {
		std::istringstream iss(line);
		if (std::getline(iss, test, '=') && pass==false) {
			// ---------------------------------
			if (iss >> amount && test==match) {
				// cout << amount << endl;
				pass=true;
			}
			else {
				pass=false;
			}
			// ---------------------------------
			if(pass) break;
		}
	}
	if (pass==false) {
		string errorout=match;
		errorout+="= argument is missing in the input file!";
		throw std::invalid_argument(errorout);
	}
	cout << match << " = " << amount << endl;
	return amount;
} // ----------

string ConstVariables::matchstring2(string file, const string match) {
	string test;
	string line;
	ifstream readFile(file);
	string amount;
	bool pass=false;

	while (std::getline(readFile, line)) {
		std::istringstream iss(line);
		if (std::getline(iss, test, '=') && pass==false) {
			// ---------------------------------
			if (iss >> amount && test==match) {
				// cout << amount << endl;
				pass=true;
			}
			else {
				pass=false;
			}
			// ---------------------------------
			if(pass) break;
		}
	}
	if (pass==false) {
		string errorout=match;
		errorout+="= argument is missing in the input file!";
		throw std::invalid_argument(errorout);
	}
	cout << match << " = " << amount << endl;
	return amount;
} // ----------

