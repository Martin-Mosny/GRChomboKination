#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// Problem specific includes
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "LoadBalance.H"
#include <H5Cpp.h>
using namespace H5;

void read_ProblemDomain(const std::string filename, ProblemDomain probDomain);

void scale_factor(HDF5Handle & a_handle, BHAMR bh_amr);

DisjointBoxLayout loadBalance(const Vector<Box> &a_grids);
void readHDF5(HDF5Handle &a_handle, int &num_levels, Vector<DisjointBoxLayout> &all_grids, Vector<LevelData<FArrayBox>> &all_states);