/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Chombo includes
#include "parstream.H" //Gives us pout()

// System includes
#include <iostream>

// Our general includes
#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"
#include "Analysis.hpp"

// Problem specific includes:
#include "KinationLevel.hpp"

// Chombo namespace
#include "UsingNamespace.H"

int runGRChomboAnalysis(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    char *in_file = argv[1];
    GRParmParse pp(argc - 2, argv + 2, NULL, in_file);
    SimulationParameters sim_params(pp);

    if (sim_params.just_check_params)
        return 0;

    // Here we now run various analysis on the results of the simulation
    // First we loop over all points
    const std::string filename = "/home/mosny/Desktop/scratch/Kination_000002.3d.hdf5";
    
    ProblemDomain probDomain;
    read_ProblemDomain(filename, probDomain);

    
    HDF5Handle a_handle(filename, HDF5Handle::OPEN_RDONLY);
    
    Vector<Box> boxes;
    std::string a_var1 = "/level_0/boxes";
    std::string &a_name1 = a_var1;
    read(a_handle, boxes, a_name1);
    pout() << "Hello" << endl;
    Vector<int> procAssign(boxes.size());
    LoadBalance(procAssign, boxes);
    DisjointBoxLayout a_layout(boxes, procAssign, probDomain);
    LevelData<FArrayBox> a_data;
    std::string a_variable = "/level_0/data";
    std::string &a_name = a_variable;
    read(a_handle, a_data, a_name, a_layout);

    a_handle.close();

    pout() << a_data.nComp() << endl;

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runGRChomboAnalysis(argc, argv);

    if (status == 0)
        pout() << "GRChombo finished." << std::endl;
    else
        pout() << "GRChombo failed with return code " << status << std::endl;

    mainFinalize();
    return status;
}
