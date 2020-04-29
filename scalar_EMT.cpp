#include <Grid/Grid.h>
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

#include "Input_Params.h"
#include "Utils.cpp"


using namespace std;
using namespace Grid;
using namespace Hadrons;

const double pi = 3.14159265358979323846;
const double alpha = 7.2973525693e-3; // PDG value
const double e = sqrt(4*pi*alpha);


int main(int argc, char** argv) {
	if (argc <2) {
		std::cout << "usage: " << argv[0] << " <parameter file>" << std::endl;
		return EXIT_FAILURE;
	}

	Omega_QED_Par par;
	XmlReader reader(argv[1]);

	read(reader,  "trajPar",    par.trajPar);
	read(reader,  "timePar",    par.timePar);
	read(reader,  "IOPar",      par.ioPar);
	read(reader,  "sourcesPar", par.sourcesPar);
	read(reader,  "quarksPar",  par.quarksPar);
	read(reader,  "contractionPar",  par.contractionPar);


	unsigned int configStart  = par.trajPar.start;
	unsigned int configEnd    = par.trajPar.end;
	unsigned int configStep   = par.trajPar.step;

	unsigned int tiStart  = par.timePar.tiStart;
	unsigned int tiEnd    = par.timePar.tiEnd;
	unsigned int tiStep   = par.timePar.tiStep;

	string runName       = par.ioPar.runName;
	string outPath       = par.ioPar.outputStem;
	string scheduleFile  = par.ioPar.scheduleFile;
	
	bool doStochQED         = par.contractionPar.doStochQED;
	bool doPertQED          = par.contractionPar.doPertQED;
	bool doScalarInsertion  = par.contractionPar.doScalarInsertion;
	bool doBaryons          = par.contractionPar.doBaryons;
	bool doMesons           = par.contractionPar.doMesons;
	bool doNParity          = par.contractionPar.doNParity;

	if ( (!doBaryons && !doMesons) ) {
		cout << "Nothing to be done" << endl;
		return 0;
	}

	vector<string> flavour; 
	vector<string> epack;
	vector<double> mass;
	vector<int>    Ls;
	vector<double> M5;
	vector<double> scale;
	vector<double> residual;
	for (int i=0; i<par.quarksPar.quarks.size(); i++) {
		flavour.push_back(  par.quarksPar.quarks[i].type     );
		epack.push_back(         ""                          );
		mass.push_back(     par.quarksPar.quarks[i].mass     );
		Ls.push_back(       par.quarksPar.quarks[i].Ls       );
		M5.push_back(       par.quarksPar.quarks[i].M5       );
		scale.push_back(    par.quarksPar.quarks[i].scale    );
		residual.push_back( par.quarksPar.quarks[i].residual );
	}

	vector<string> source;
	vector<double> width;
	vector<string> source_name;
	for (int i=0; i<par.sourcesPar.sources.size(); i++) {
		source.push_back(  par.sourcesPar.sources[i].type  );
		width.push_back(   par.sourcesPar.sources[i].width );
		source_name.push_back(   source[i]+"_w"+to_string(width[i]) );
	}

	string boundary  = "1 1 1 -1";
	string twist     = "0. 0. 0. 0.";

	string gaugeFile = par.ioPar.gaugeFile;


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
	LOG(Message) << par.timePar;
	LOG(Message) << par.ioPar;
	LOG(Message) << par.sourcesPar;
	LOG(Message) << par.quarksPar;


	// Global Parameters
	Application application;
	Application::GlobalPar globalPar;
	globalPar.runId              = runName;
	globalPar.trajCounter.start  = configStart;
	globalPar.trajCounter.end    = configEnd;
	globalPar.trajCounter.step   = configStep;
	application.setPar(globalPar);


	unsigned int nt;
	// Setup Gauge field
	if (gaugeFile == "") {
		application.createModule<MGauge::Unit>("gauge");
		nt = GridDefaultLatt()[Tp];
	} else {
		MIO::LoadNersc::Par loadPar;
		loadPar.file       = gaugeFile;
		application.createModule<MIO::LoadNersc>("gauge", loadPar);
		nt = GridDefaultLatt()[Tp];
	}

	// Fix gauge (coulomb)
	MGauge::GaugeFix::Par gaugeFixPar;
	gaugeFixPar.gauge      = "gauge";
	gaugeFixPar.alpha      = 0.05;
	gaugeFixPar.maxiter    = 1000000;
	gaugeFixPar.Omega_tol  = 1e-8;
	gaugeFixPar.Phi_tol    = 1e-8;
	gaugeFixPar.gaugeFix   = MGauge::Fix::coulomb;
	gaugeFixPar.Fourier    = true;
	application.createModule<MGauge::GaugeFix>("gaugeFix", gaugeFixPar);

	// Cast gauge fields to single precision
	MUtilities::GaugeSinglePrecisionCast::Par gaugefPar;

	gaugefPar.field     = "gauge";
	application.createModule<MUtilities::GaugeSinglePrecisionCast>("gaugef", gaugefPar);

	gaugefPar.field     = "gaugeFix";
	application.createModule<MUtilities::GaugeSinglePrecisionCast>("gaugefFix", gaugefPar);

if(doStochQED || doPertQED) {
	// Stochastic photon field
	MGauge::StochEm::Par photonPar;
	photonPar.gauge =     PhotonR::Gauge::feynman;
	photonPar.zmScheme =  PhotonR::ZmScheme::qedL;
	application.createModule<MGauge::StochEm>("photon_field", photonPar);
}
if (doStochQED) {
	// Electrify gauge field with +/-e (q=-1/3 for strange quark)
	MGauge::Electrify::Par electrifyPar;
	electrifyPar.gauge    = "gaugeFix";
	electrifyPar.emField  = "photon_field";
	electrifyPar.e        = e;
	electrifyPar.charge   = -1.0/3;
	application.createModule<MGauge::Electrify>("gaugeFix_stoch+",electrifyPar);

	electrifyPar.gauge    = "gaugeFix";
	electrifyPar.emField  = "photon_field";
	electrifyPar.e        = -e;
	electrifyPar.charge   = -1.0/3;
	application.createModule<MGauge::Electrify>("gaugeFix_stoch-",electrifyPar);

	gaugefPar.field     = "gaugeFix_stoch+";
	application.createModule<MUtilities::GaugeSinglePrecisionCast>("gaugefFix_stoch+", gaugefPar);

	gaugefPar.field     = "gaugeFix_stoch-";
	application.createModule<MUtilities::GaugeSinglePrecisionCast>("gaugefFix_stoch-", gaugefPar);
}

	// Setup action (DWF) and solver (mixed precision)
	for (int i=0; i<flavour.size(); i++) {
		// Actions double precision
		MAction::ScaledDWF::Par actionPar;

		actionPar.gauge     = "gaugeFix";
		actionPar.Ls        = Ls[i];
		actionPar.M5        = M5[i];
		actionPar.mass      = mass[i];
		actionPar.scale     = scale[i];
		actionPar.boundary  = boundary;
		actionPar.twist     = twist;
		application.createModule<MAction::ScaledDWF>("dwf_stoch0_"+flavour[i], actionPar);

		// Actions single precision
		MAction::ScaledDWF::Par actionfPar;

		actionfPar.gauge     = "gaugefFix";
		actionfPar.Ls        = Ls[i];
		actionfPar.M5        = M5[i];
		actionfPar.mass      = mass[i];
		actionfPar.scale      = scale[i];
		actionfPar.boundary  = boundary;
		actionfPar.twist     = twist;
		application.createModule<MAction::ScaledDWFF>("dwff_stoch0_"+flavour[i], actionfPar);
if (doStochQED) {
		actionPar.gauge     = "gaugeFix_stoch+";
		actionPar.Ls        = Ls[i];
		actionPar.M5        = M5[i];
		actionPar.mass      = mass[i];
		actionPar.scale     = scale[i];
		actionPar.boundary  = boundary;
		actionPar.twist     = twist;
		application.createModule<MAction::ScaledDWF>("dwf_stoch+_"+flavour[i], actionPar);

		actionPar.gauge     = "gaugeFix_stoch-";
		actionPar.Ls        = Ls[i];
		actionPar.M5        = M5[i];
		actionPar.mass      = mass[i];
		actionPar.scale     = scale[i];
		actionPar.boundary  = boundary;
		actionPar.twist     = twist;
		application.createModule<MAction::ScaledDWF>("dwf_stoch-_"+flavour[i], actionPar);


		actionfPar.gauge     = "gaugefFix_stoch+";
		actionfPar.Ls        = Ls[i];
		actionfPar.M5        = M5[i];
		actionfPar.mass      = mass[i];
		actionfPar.scale     = scale[i];
		actionfPar.boundary  = boundary;
		actionfPar.twist     = twist;
		application.createModule<MAction::ScaledDWFF>("dwff_stoch+_"+flavour[i], actionfPar);

		actionfPar.gauge     = "gaugefFix_stoch-";
		actionfPar.Ls        = Ls[i];
		actionfPar.M5        = M5[i];
		actionfPar.mass      = mass[i];
		actionfPar.scale     = scale[i];
		actionfPar.boundary  = boundary;
		actionfPar.twist     = twist;
		application.createModule<MAction::ScaledDWFF>("dwff_stoch-_"+flavour[i], actionfPar);
}

		// Solvers
		MSolver::MixedPrecisionRBPrecCG::Par solverPar;

		solverPar.innerAction        = "dwff_stoch0_"+flavour[i];
		solverPar.outerAction        = "dwf_stoch0_"+flavour[i];
		solverPar.maxInnerIteration  = 30000;
		solverPar.maxOuterIteration  = 100;
		solverPar.residual           = residual[i];
		solverPar.eigenPack          = epack[i];
		application.createModule<MSolver::MixedPrecisionRBPrecCG>("mcg_stoch0_"+flavour[i], solverPar);

if (doStochQED) {
		solverPar.innerAction        = "dwff_stoch+_"+flavour[i];
		solverPar.outerAction        = "dwf_stoch+_"+flavour[i];
		solverPar.maxInnerIteration  = 30000;
		solverPar.maxOuterIteration  = 100;
		solverPar.residual           = residual[i];
		solverPar.eigenPack          = epack[i];
		application.createModule<MSolver::MixedPrecisionRBPrecCG>("mcg_stoch+_"+flavour[i], solverPar);

		solverPar.innerAction        = "dwff_stoch-_"+flavour[i];
		solverPar.outerAction        = "dwf_stoch-_"+flavour[i];
		solverPar.maxInnerIteration  = 30000;
		solverPar.maxOuterIteration  = 100;
		solverPar.residual           = residual[i];
		solverPar.eigenPack          = epack[i];
		application.createModule<MSolver::MixedPrecisionRBPrecCG>("mcg_stoch-_"+flavour[i], solverPar);
}
	}

	// Sources
	for (int ti=tiStart; ti<tiEnd; ti+=tiStep) {
		for (int k=0; k < source.size(); k++){
			if (source[k] == "point") {
				MSource::Point::Par sourcePointPar;
				sourcePointPar.position       = "0 0 0 "+to_string(ti);
				application.createModule<MSource::Point>(
					"source_"+source_name[k]+"_ti"+to_string(ti), sourcePointPar);
			} else if (source[k] == "gauss") {
				MSource::Gauss::Par sourceGaussPar;
				sourceGaussPar.position  = "0 0 0";
				sourceGaussPar.mom       = "0 0 0 0";
				sourceGaussPar.tA        = ti;
				sourceGaussPar.tB        = ti;
				sourceGaussPar.width     = width[k];
				application.createModule<MSource::Gauss>(
					"source_"+source_name[k]+"_ti"+to_string(ti), sourceGaussPar);
			} else if (source[k] == "wall") {
				MSource::Wall::Par sourceWallPar;
				sourceWallPar.tW         = ti;
				sourceWallPar.mom        = "0. 0. 0. 0.";
				application.createModule<MSource::Wall>(
					"source_"+source_name[k]+"_ti"+to_string(ti), sourceWallPar);
			} else {
				LOG(Message) << "Source type \"" << source[k] << "\" not valid";
				return EXIT_FAILURE;
			}
		}
	}

	// sink
	MSink::Point::Par sinkPointPar;
	sinkPointPar.mom  = "0 0 0";
	application.createModule<MSink::ScalarPoint>("sink_scalar_pnt", sinkPointPar);

	// Propagators

	/*
	*	(1) prop_e0 --------                                      (10) prop_e+      (12) prop_e-
	*	      |             \                                           |                |
	*         |              \                                          |                |
	*	(2) prop_S_e0     (5) prop_A_e0 -----------               (11) prop_S_e+    (13) prop_S_e-
	*	      |                \                   \
	*         |                 \                   \
	*	(3) prop_SA_e0       (6) prop_AS_e0      (8) prop_AA_e0
	*	      |                   \                   \
	*         |                    \                   \
	*	(4) prop_SAA_e0         (7) prop_ASA_e0      (9) prop_AAS_e0
	*
	*/

	vector<vector<map<int,vector<string>>>> prop_names(flavour.size());

	for (int i=0; i<flavour.size(); i++) {
		prop_names[i].resize(source.size());

		for (int k=0; k<source.size(); k++) {
			for (int ti=tiStart; ti<tiEnd; ti+=tiStep) {
				vector<string> tmp_prop_names;
				// (1) prop_stoch0_N
				createProp(application,
					flavour[i]+"_"+source_name[k]+"_stoch0_N_ti"+to_string(ti),
					tmp_prop_names,
					nt);
if (doScalarInsertion) {
				// (2) prop_stoch0_S
				createProp(application,
					flavour[i]+"_"+source_name[k]+"_stoch0_S_ti"+to_string(ti),
					tmp_prop_names,
					nt);
}
if (doPertQED && doScalarInsertion) {
				// (3) prop_stoch0_SA
				createProp(application,
					flavour[i]+"_"+source_name[k]+"_stoch0_SA_ti"+to_string(ti),
					tmp_prop_names,
					nt);
}
if (doPertQED && doScalarInsertion) {
				// (4) prop_stoch0_SAA
				createProp(application,
					flavour[i]+"_"+source_name[k]+"_stoch0_SAA_ti"+to_string(ti),
					tmp_prop_names,
					nt);
}
if (doPertQED) {
				// (5) prop_stoch0_A
				createProp(application,
					flavour[i]+"_"+source_name[k]+"_stoch0_A_ti"+to_string(ti),
					tmp_prop_names,
					nt);
}
if (doPertQED && doScalarInsertion) {
				// (6) prop_stoch0_AS
				createProp(application,
					flavour[i]+"_"+source_name[k]+"_stoch0_AS_ti"+to_string(ti),
					tmp_prop_names,
					nt);
}
if (doPertQED && doScalarInsertion) {
				// (7) prop_stoch0_ASA
				createProp(application,
					flavour[i]+"_"+source_name[k]+"_stoch0_ASA_ti"+to_string(ti),
					tmp_prop_names,
					nt);
}
if (doPertQED) {
				// (8) prop_stoch0_AA
				createProp(application,
					flavour[i]+"_"+source_name[k]+"_stoch0_AA_ti"+to_string(ti),
					tmp_prop_names,
					nt);
}
if (doPertQED && doScalarInsertion) {
				// (9) prop_stoch0_AAS
				createProp(application,
					flavour[i]+"_"+source_name[k]+"_stoch0_AAS_ti"+to_string(ti),
					tmp_prop_names,
					nt);
}
if (doStochQED) {
				// (10) prop_stoch+_N
				createProp(application,
					flavour[i]+"_"+source_name[k]+"_stoch+_N_ti"+to_string(ti),
					tmp_prop_names,
					nt);
}
if (doStochQED && doScalarInsertion) {
				// (11) prop_stoch+_S
				createProp(application,
					flavour[i]+"_"+source_name[k]+"_stoch+_S_ti"+to_string(ti),
					tmp_prop_names,
					nt);
}
if (doStochQED) {
				// (12) prop_stoch-_N
				createProp(application,
					flavour[i]+"_"+source_name[k]+"_stoch-_N_ti"+to_string(ti),
					tmp_prop_names,
					nt);
}
if (doStochQED && doScalarInsertion) {
				// (13) prop_stoch-_S
				createProp(application,
					flavour[i]+"_"+source_name[k]+"_stoch-_S_ti"+to_string(ti),
					tmp_prop_names,
					nt);
}
				prop_names[i][k][ti] = tmp_prop_names;
			}
		}
	}

	for (int i=0; i<flavour.size(); i++) {
		for (int k=0; k<source.size(); k++) {
			for (int ti=tiStart; ti<tiEnd; ti+=tiStep) {

				// Contract Baryons
				int cB=0;
if (doBaryons) {
				for (int l=0; l<prop_names[i][k][ti].size(); l++) {
				for (int m=0; m<prop_names[i][k][ti].size(); m++) {
				for (int n=0; n<prop_names[i][k][ti].size(); n++) {
					std::regex rex("([a-zA-Z])_([a-z0-9._]+)_stoch([0+-])_([A-Z]+)_ti([0-9]+)");
					std::smatch smRex1;
					std::smatch smRex2;
					std::smatch smRex3;
					std::regex_search(prop_names[i][k][ti][l], smRex1, rex);
					std::regex_search(prop_names[i][k][ti][m], smRex2, rex);
					std::regex_search(prop_names[i][k][ti][n], smRex3, rex);

					if (smRex1[3] != smRex2[3] || smRex2[3] != smRex3[3]) // Check all have the same stocastic charge
						continue;

					string stoch  = smRex1[3].str();
					string p1     = smRex1[4].str();
					string p2     = smRex2[4].str();
					string p3     = smRex3[4].str();

					// Check if 0 or 1 scalar insertions
					int nScalar = std::count(p1.begin(), p1.end(), 'S')
								+std::count(p2.begin(), p2.end(), 'S')
								+std::count(p3.begin(), p3.end(), 'S');

					if (nScalar != 0 && nScalar != 1)
						continue;

					// Check if 0 or 2 Aslash insertions
					int nAslash = std::count(p1.begin(), p1.end(), 'A')
								+std::count(p2.begin(), p2.end(), 'A')
								+std::count(p3.begin(), p3.end(), 'A');

					if (nAslash != 0 && nAslash != 2)
						continue;

					createBaryon(application,
						flavour[i]+flavour[i]+flavour[i]+"_"+source_name[k]+"_stoch"+stoch+"_"+p1+"_"+p2+"_"+p3+"_ti"+to_string(ti)+"_P+",
						"sink_scalar_pnt",
						outPath);

					if (doNParity) {
						createBaryon(application,
							flavour[i]+flavour[i]+flavour[i]+"_"+source_name[k]+"_stoch"+stoch+"_"+p1+"_"+p2+"_"+p3+"_ti"+to_string(ti)+"_P-",
							"sink_scalar_pnt",
							outPath);
						cB++;
					}					
					cB++;
				}}}
}

				// Contraction Mesons
				int cM=0;
if (doMesons) {
				for (int l=0; l<prop_names[i][k][ti].size(); l++) {
				for (int m=0; m<prop_names[i][k][ti].size(); m++) {
					std::regex rex("([a-zA-Z])_([a-z0-9._]+)_stoch([0+-])_([A-Z]+)_ti([0-9]+)");
					std::smatch smRex1;
					std::smatch smRex2;
					std::regex_search(prop_names[i][k][ti][l], smRex1, rex);
					std::regex_search(prop_names[i][k][ti][m], smRex2, rex);

					if (smRex1[3] != smRex2[3]) // Check both have the same stocastic charge
						continue;

					string stoch  = smRex1[3].str();
					string p1     = smRex1[4].str();
					string p2     = smRex2[4].str();

					// Check if 0 or 1 scalar insertions
					int nScalar = std::count(p1.begin(), p1.end(), 'S')
								+std::count(p2.begin(), p2.end(), 'S');

					if (nScalar != 0 && nScalar != 1)
						continue;

					// Check if 0 or 2 Aslash insertions
					int nAslash = std::count(p1.begin(), p1.end(), 'A')
								 +std::count(p2.begin(), p2.end(), 'A');

					if (nAslash != 0 && nAslash != 2)
						continue;

					createMeson(application,
						flavour[i]+flavour[i]+"_"+source_name[k]+"_stoch"+stoch+"_"+p1+"_"+p2+"_ti"+to_string(ti),
						"sink_scalar_pnt",
						outPath);
					cM++;
				}}
}
			}
		}
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
