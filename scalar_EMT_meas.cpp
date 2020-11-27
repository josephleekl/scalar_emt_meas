#include <Grid/Grid.h>
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

#include "Input_Params.h"
#include "Utils.cpp"


using namespace std;
using namespace Grid;
using namespace Hadrons;

const double pi = 3.14159265358979323846;


int main(int argc, char** argv) {
	if (argc <2) {
		std::cout << "usage: " << argv[0] << " <parameter file>" << std::endl;
		return EXIT_FAILURE;
	}

	Scalar_EMT_Par par;
	XmlReader reader(argv[1]);

	read(reader,  "trajPar",    par.trajPar);
	read(reader,  "IOPar",      par.ioPar);
	read(reader,  "ConfigPar",      par.configPar);

	// Inititalise Grid
	Grid_init(&argc, &argv);
	HadronsLogError.Active(GridLogError.isActive());
	HadronsLogWarning.Active(GridLogWarning.isActive());
	HadronsLogMessage.Active(GridLogMessage.isActive());
	HadronsLogIterative.Active(GridLogIterative.isActive());
	HadronsLogDebug.Active(GridLogDebug.isActive());
	LOG(Message) << "Grid initialized" << std::endl;

	LOG(Message) << "======= Parameters from \"" << argv[1] << "\" =======" << endl;
 	LOG(Message) << par.trajPar;
	LOG(Message) << par.ioPar;
	LOG(Message) << par.configPar;

	//Run parameter
	string runName       = par.ioPar.runName;
	string outPath       = par.ioPar.outputStem;
	string scheduleFile  = par.ioPar.scheduleFile;
	//string scheduleFile  = "";
	unsigned int configStart  = par.trajPar.start;
	unsigned int configEnd    = par.trajPar.end;
	unsigned int configStep   = par.trajPar.step;


	// Global Parameters
	Application application;
	Application::GlobalPar globalPar;
	globalPar.runId              = runName;
	globalPar.trajCounter.start  = configStart;
	globalPar.trajCounter.end    = configEnd;
	globalPar.trajCounter.step   = configStep;
	application.setPar(globalPar);


	//PIflowtimes = 1/(g*sqrt(t))
	//flowtime = t
	std::vector<std::string> PIflowtimes;
	if (par.configPar.g == 0.1)
	{
		PIflowtimes ={"1.0", "2.0", "3.0", "4.0", "5.0", "6.0", "7.0", "8.0", "9.0"};
	} else if (par.configPar.g == 0.2)
	{
		PIflowtimes ={"1.0" ,"1.4" ,"1.8", "2.2", "2.6", "3.0", "3.4", "3.8", "4.2"};
	} else if (par.configPar.g == 0.3)
	{
		PIflowtimes ={"1.0", "1.3" ,"1.6" ,"1.9", "2.2", "2.5", "2.8", "3.1", "3.4"};
	}

	unsigned int nt;

	//Load SU(N) field
	std::string config_filename = par.ioPar.configFile;
	MIO::LoadCosmHolSU2::Par loadPar;
	loadPar.file = config_filename;
	application.createModule<MIO::LoadCosmHolSU2>("phi", loadPar); 

	//Measurement: trphi
	MScalarSUN::TrPhiSU2::Par trphipar;
	trphipar.field = "phi";
	trphipar.maxPow = 6;
	trphipar.output = "onept/cosmhol-trphi";
	application.createModule<MScalarSUN::TrPhiSU2>("trphi", trphipar);
	
	//Measurement: trkinc
	MScalarSUN::TrKineticSU2::Par trkincpar;
	trkincpar.field = "phi";
	trkincpar.output = "onept/cosmhol-trkinc";
	trkincpar.type = 3; // Difftype 3 = central
	application.createModule<MScalarSUN::TrKineticSU2>("trkinc", trkincpar);

	//Measurement: impc
	MScalarSUN::TransProjSU2::Par impcpar;
	impcpar.op = "trphi_2";
	impcpar.output = "onept/cosmhol-impc";
	impcpar.type = 3; // Difftype 3 = central
	application.createModule<MScalarSUN::TransProjSU2>("impc", impcpar);

	//Measurement: emtc
	MScalarSUN::EMTSU2::Par emtcpar;
	emtcpar.kinetic = "trkinc";
	emtcpar.phiPow = "trphi";
	emtcpar.improvement = "impc";
	emtcpar.m2 = par.configPar.m2;
	emtcpar.g = par.configPar.g;
	emtcpar.xi = 0.;
	emtcpar.lambda = 1.0;
	emtcpar.output = "onept/cosmhol-emtc";
	application.createModule<MScalarSUN::EMTSU2>("emtc", emtcpar);


	//Measurement: wilson flowtime 
	vector<MScalarSUN::WilsonFlowSU2::Par>   wilson_field_par(PIflowtimes.size());
	vector<MScalarSUN::TrPhiSU2::Par>        wilson_trphi_par(PIflowtimes.size());	
 	vector<MScalarSUN::TwoPointSU2::Par> wilson_twopt_par_dir0(PIflowtimes.size());
 	vector<MScalarSUN::TwoPointSU2::Par> wilson_twopt_par_dir1(PIflowtimes.size());
 	vector<MScalarSUN::TwoPointSU2::Par> wilson_twopt_par_dir2(PIflowtimes.size());


	 for (int i=0; i<PIflowtimes.size(); i++)
	{
		wilson_field_par[i].field = "phi"; 
		wilson_field_par[i].flowtime = pow(stod(PIflowtimes[i]) * par.configPar.g , -2);
		application.createModule<MScalarSUN::WilsonFlowSU2>("wilson_"+PIflowtimes[i]+"_phi", wilson_field_par[i]);

		wilson_trphi_par[i].field = "wilson_"+PIflowtimes[i]+"_phi";
		wilson_trphi_par[i].maxPow = 6;
		wilson_trphi_par[i].output = "wilson_onept/cosmhol-trphi";
		application.createModule<MScalarSUN::TrPhiSU2>("wilson_"+PIflowtimes[i]+"_trphi", wilson_trphi_par[i]);

		wilson_twopt_par_dir0[i].op.push_back(std::make_pair("emtc_0_0", "wilson_"+PIflowtimes[i]+"_trphi_2"));
		wilson_twopt_par_dir1[i].op.push_back(std::make_pair("emtc_1_1", "wilson_"+PIflowtimes[i]+"_trphi_2"));
		wilson_twopt_par_dir2[i].op.push_back(std::make_pair("emtc_2_2", "wilson_"+PIflowtimes[i]+"_trphi_2"));
        wilson_twopt_par_dir0[i].op.push_back(std::make_pair("trphi_2", "wilson_"+PIflowtimes[i]+"_trphi_2"));
		wilson_twopt_par_dir1[i].op.push_back(std::make_pair("trphi_2", "wilson_"+PIflowtimes[i]+"_trphi_2"));
		wilson_twopt_par_dir2[i].op.push_back(std::make_pair("trphi_2", "wilson_"+PIflowtimes[i]+"_trphi_2"));

		if (par.configPar.L == 64)
		{
			wilson_twopt_par_dir0[i].mom.push_back("1 0");
			wilson_twopt_par_dir1[i].mom.push_back("0 1");
		} else if (par.configPar.L == 128) 
		{
			wilson_twopt_par_dir0[i].mom.push_back("2 0");
			wilson_twopt_par_dir1[i].mom.push_back("0 2");
		} else if (par.configPar.L == 256) 
		{
			wilson_twopt_par_dir0[i].mom.push_back("4 0");
			wilson_twopt_par_dir1[i].mom.push_back("0 4");
		}
		wilson_twopt_par_dir2[i].mom.push_back("0 0");


		wilson_twopt_par_dir0[i].output = "twopt/flowtime_"+PIflowtimes[i]+"/dir0/cosmhol-twopt"+PIflowtimes[i];
		wilson_twopt_par_dir1[i].output = "twopt/flowtime_"+PIflowtimes[i]+"/dir1/cosmhol-twopt"+PIflowtimes[i];
		wilson_twopt_par_dir2[i].output = "twopt/flowtime_"+PIflowtimes[i]+"/dir2/cosmhol-twopt"+PIflowtimes[i];
		application.createModule<MScalarSUN::TwoPointSU2>("wilson_twopt"+PIflowtimes[i]+"dir0", wilson_twopt_par_dir0[i]);
		application.createModule<MScalarSUN::TwoPointSU2>("wilson_twopt"+PIflowtimes[i]+"dir1", wilson_twopt_par_dir1[i]);
		application.createModule<MScalarSUN::TwoPointSU2>("wilson_twopt"+PIflowtimes[i]+"dir2", wilson_twopt_par_dir2[i]);
	} 


	// Save Xml
	mkdir(outPath+"/Xml");
	application.saveParameterFile(outPath+"/Xml/"+runName+"_"+to_string(configStart)+"-"+to_string(configEnd)+"-"+to_string(configStep)+".xml");
	
	// Scheduler
	if (scheduleFile == "") {
		application.schedule();
		application.saveSchedule(runName+".sched");
	} else {
		application.loadSchedule(scheduleFile);
	}
	
	// Execute
	application.run();

	// Epilogue
	LOG(Message) << "Grid is finalizing now" << std::endl;
	Grid_finalize();


	return EXIT_SUCCESS;
}
