#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// Problem specific includes
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "LoadBalance.H"
#include <vector>
#include <algorithm>
#include <H5Cpp.h>
using namespace H5;

#ifndef SIM_PARAMS_CLASS
#define SIM_PARAMS_CLASS
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
#endif

void read_HDF5_key_attributes(const std::string &filename, sim_parameters &a_params);

void read_HDF5_LevelProblemDomain(const std::string &filename, ProblemDomain &probDomain, int &level_int);
void read_HDF5_ProblemDomain(const std::string &filename, Vector<ProblemDomain> &probDomain, sim_parameters &a_params);
void read_grid(const std::string &filename, Vector<DisjointBoxLayout> &a_grid);

void read_HDF5_LevelData(const std::string &filename, ProblemDomain &probDomain, LevelData<FArrayBox> &a_data, int &level_int);
void read_HDF5_Data(const std::string &filename, Vector<LevelData<FArrayBox> *> &a_data, Vector<DisjointBoxLayout> &a_grid, sim_parameters &a_params);
void add_ghosts(Vector<LevelData<FArrayBox> *> &a_data, const IntVect &ghost_vect);
void fill_ghosts(Vector<LevelData<FArrayBox> *> &a_data, Vector<DisjointBoxLayout> &a_grid, Vector<ProblemDomain> &probDomain, sim_parameters &a_params);
void comp_field(Vector<LevelData<FArrayBox> *> &a_data, Vector<LevelData<FArrayBox> *> &a_output, int component);
void writing_to_notepad(const std::string filename, const Vector<Vector<Real>> &a_data, const std::vector<std::string> &headers);