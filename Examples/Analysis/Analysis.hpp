#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// Problem specific includes
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "LoadBalance.H"
#include <H5Cpp.h>
using namespace H5;

class sim_parameters
{
    public:
        int max_level;
        int num_levels;
        Real time;
        int num_components;
        Vector<Real> dx;
        Vector<int> ref_ratio;
};

void read_HDF5_key_attributes(const std::string &filename, sim_parameters &a_params);

void read_HDF5_LevelProblemDomain(const std::string &filename, ProblemDomain &probDomain, int &level_int);
void read_HDF5_ProblemDomain(const std::string &filename, Vector<ProblemDomain> &probDomain, sim_parameters &a_params);

void read_HDF5_LevelData(const std::string &filename, ProblemDomain &probDomain, LevelData<FArrayBox> &a_data, int &level_int);
void read_HDF5_Data(const std::string &filename, Vector<ProblemDomain> &probDomain, Vector<LevelData<FArrayBox> *> &a_data, sim_parameters &a_params);

void scale_factor(HDF5Handle & a_handle, BHAMR bh_amr);