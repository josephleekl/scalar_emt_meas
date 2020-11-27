using namespace std;
using namespace Grid;
using namespace Hadrons;

namespace Scalar_EMT_IO {
	class TrajPar : Serializable {
	public:
		GRID_SERIALIZABLE_CLASS_MEMBERS(TrajPar,
			unsigned int, start,
			unsigned int, end,
			unsigned int, step
		);
	};


	class IOPar : Serializable {
	public:
		GRID_SERIALIZABLE_CLASS_MEMBERS(IOPar,
			std::string,  runName,
			std::string,  outputStem,
			std::string,  scheduleFile,
			std::string,  configFile
		);
	};

	class ConfigPar : Serializable {
	public:
		GRID_SERIALIZABLE_CLASS_MEMBERS(ConfigPar,
			double, g,
			double, m2,
			unsigned int, L,
			unsigned int, N
		);
	};

}



struct Scalar_EMT_Par{
	Scalar_EMT_IO::TrajPar         trajPar;
	Scalar_EMT_IO::IOPar           ioPar;
	Scalar_EMT_IO::ConfigPar       configPar;

};
