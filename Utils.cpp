using namespace std;
using namespace Grid;
using namespace Hadrons;

void createProp(	Application& app,
					string propName,
					std::vector<string>& prop_names,
					int nt) {

	std::regex rex("([a-zA-Z])_([a-z0-9._]+)_stoch([0+-])_([A-Z]+)_ti([0-9]+)");
	std::smatch smRex;
	std::regex_search(propName, smRex, rex);

	string q      = smRex[1].str();
	string source = smRex[2].str();
	string stoch  = smRex[3].str();
	string p      = smRex[4].str();
	string ti     = smRex[5].str();
	
	MFermion::GaugeProp::Par quarkPar;

	quarkPar.solver    = "mcg_stoch"+stoch+"_"+q;

	if (p == "N") {
		quarkPar.source    = "source_"+source+"_ti"+ti;
	} else {
		string preType;
		string postType;
		if (p.size() == 1) {
			preType = "N";
			postType = p;
		} else {
			preType = p.substr(0,p.size()-1);
			postType = p.substr(p.size()-1,1);
		}
		string preProp = q+"_"+source+"_stoch"+stoch+"_"+preType+"_ti"+ti;

		if (postType == "S") {
			quarkPar.source    = "prop_"+preProp;
		} else if (postType == "A") {
			MSource::SeqAslash::Par  ASlashSourcePar;
			ASlashSourcePar.tA = 0;
			ASlashSourcePar.tB = nt-1;
			ASlashSourcePar.emField = "photon_field";
			ASlashSourcePar.mom = "0. 0. 0. 0.";

			ASlashSourcePar.q  = "prop_"+preProp;
			app.createModule<MSource::SeqAslash>(
				"source_"+preProp, ASlashSourcePar);

			quarkPar.source    = "source_"+preProp;
		}
	}

	app.createModule<MFermion::GaugeProp>("prop_"+propName, quarkPar);
	
	prop_names.push_back("prop_"+propName);
}

void createBaryon(	Application& app,
					string conName,
					string sink,
					string outPath) {

	std::regex rex("([a-zA-Z]+)_([a-z0-9._]+)_stoch([0+-])_([A-Z]+)_([A-Z]+)_([A-Z]+)_ti([0-9]+)_P([+-])");
	std::smatch smRex;
	std::regex_search(conName, smRex, rex);


	string quarks = smRex[1].str();
	string q1(      1,quarks[0]);
	string q2(      1,quarks[1]);
	string q3(      1,quarks[2]);
	string source = smRex[2].str();
	string stoch  = smRex[3].str();
	string p1     = smRex[4].str();
	string p2     = smRex[5].str();
	string p3     = smRex[6].str();
	string ti     = smRex[7].str();
	string parity = smRex[8].str();
	
	MContraction::Baryon::Par barPar;
	barPar.gammas      = "(j12 j12) (j32X j32X) (j32Y j32Y) (j32Z j32Z)";
	barPar.quarks      = quarks;
	barPar.prefactors  = "1.0";
	barPar.sink        = sink;
	barPar.parity      = (parity == "-")? "-1" : "+1";

	barPar.output = outPath+"/Baryon/"+quarks+"_"+source+"_stoch"+stoch+"_P"+parity+"/"+p1+"_"+p2+"_"+p3+"/"+conName;
	barPar.q1          = "prop_"+q1+"_"+source+"_stoch"+stoch+"_"+p1+"_ti"+ti;
	barPar.q2          = "prop_"+q2+"_"+source+"_stoch"+stoch+"_"+p2+"_ti"+ti;
	barPar.q3          = "prop_"+q3+"_"+source+"_stoch"+stoch+"_"+p3+"_ti"+ti;
	app.createModule<MContraction::Baryon>("baryon_"+conName, barPar);

}

void createMeson(	Application& app,
					string conName,
					string sink,
					string outPath) {

	std::regex rex("([a-zA-Z]+)_([a-z0-9._]+)_stoch([0+-])_([A-Z]+)_([A-Z]+)_ti([0-9]+)");
	std::smatch smRex;
	std::regex_search(conName, smRex, rex);

	string quarks = smRex[1].str();
	string q1(      1,quarks[0]);
	string q2(      1,quarks[1]);
	string source = smRex[2].str();
	string stoch  = smRex[3].str();
	string p1     = smRex[4].str();
	string p2     = smRex[5].str();
	string ti     = smRex[6].str();
	
	MContraction::Meson::Par mesonPar;
	mesonPar.gammas    = "all";
	mesonPar.sink        = sink;

	mesonPar.q1          = "prop_"+q1+"_"+source+"_stoch"+stoch+"_"+p1+"_ti"+ti;
	mesonPar.q2          = "prop_"+q2+"_"+source+"_stoch"+stoch+"_"+p2+"_ti"+ti;
	
	string p2r = p2;
	std::reverse(p2.begin(), p2.end());

	mesonPar.output = outPath+"/Meson/"+quarks+"_"+source+"_stoch"+stoch+"/"+p1+"_"+p2r+"/"+conName;
	
	app.createModule<MContraction::Meson>("meson_"+conName, mesonPar);

}