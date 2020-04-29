using namespace std;
using namespace Grid;
using namespace Hadrons;

namespace Omega_QED_IO {
	class TrajPar : Serializable {
	public:
		GRID_SERIALIZABLE_CLASS_MEMBERS(TrajPar,
			unsigned int, start,
			unsigned int, end,
			unsigned int, step
		);
	};

	class TimePar : Serializable {
	public:
		GRID_SERIALIZABLE_CLASS_MEMBERS(TimePar,
			unsigned int, tiStart,
			unsigned int, tiEnd,
			unsigned int, tiStep
		);
	};

	class IOPar : Serializable {
	public:
		GRID_SERIALIZABLE_CLASS_MEMBERS(IOPar,
			std::string,  runName,
			std::string,  gaugeFile,
			std::string,  outputStem,
			std::string,  scheduleFile
		);
	};

	class SourcePar :Serializable {
	public:
		GRID_SERIALIZABLE_CLASS_MEMBERS(SourcePar,
			std::string, type,
			double,      width
		);
	};

	class SourcesPar :Serializable {
	public:
		GRID_SERIALIZABLE_CLASS_MEMBERS(SourcesPar,
			std::vector<SourcePar>, sources
		);
	};

	class QuarkPar : Serializable {
	public:
		GRID_SERIALIZABLE_CLASS_MEMBERS(QuarkPar,
			std::string, type,
			double,      mass,
			int,         Ls,
			double,      M5,
			double,      scale,
			double,      residual
		);
	};

	class QuarksPar : Serializable {
	public:
		GRID_SERIALIZABLE_CLASS_MEMBERS(QuarksPar,
			std::vector<QuarkPar>, quarks
		);
	};

	class ContractionPar : Serializable {
	public:
		GRID_SERIALIZABLE_CLASS_MEMBERS(ContractionPar,
			bool, doStochQED,
			bool, doPertQED,
			bool, doScalarInsertion,
			bool, doBaryons,
			bool, doMesons,
			bool, doNParity
		);
	};
}



struct Omega_QED_Par{
	Omega_QED_IO::TrajPar         trajPar;
	Omega_QED_IO::TimePar         timePar;
	Omega_QED_IO::IOPar           ioPar;
	Omega_QED_IO::SourcesPar      sourcesPar;
	Omega_QED_IO::QuarksPar       quarksPar;
	Omega_QED_IO::ContractionPar  contractionPar;
};
