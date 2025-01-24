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
#include "InputOutput.hpp"
using namespace H5;

void create_layout_data(Vector<LevelData<FArrayBox> *> &a_output, Vector<LevelData<FArrayBox> *> &a_input, Real &val);
Real average_integral(Vector<LevelData<FArrayBox> *> &a_data, sim_parameters &a_params, int &comp);
Real find_max_value(Vector<LevelData<FArrayBox> *> &a_data, int &comp);
Real find_diff_value(Vector<LevelData<FArrayBox> *> &a_data, int &comp);
void rho_calculator(Vector<LevelData<FArrayBox> *> &rho, Vector<LevelData<FArrayBox> *> &a_data, Vector<LevelData<FArrayBox> *> &a_phi, sim_parameters &a_param);
Vector<LevelData<FArrayBox> *> Ham_calc(const Vector<DisjointBoxLayout> &grids, Vector<LevelData<FArrayBox> *> a_data, Vector<LevelData<FArrayBox> *> &a_phi, const sim_parameters &a_params);
Vector<LevelData<FArrayBox> *> normalized_Ham_calc(const Vector<DisjointBoxLayout> &grids, Vector<LevelData<FArrayBox> *> a_data, Vector<LevelData<FArrayBox> *> &a_phi, const sim_parameters &a_params);