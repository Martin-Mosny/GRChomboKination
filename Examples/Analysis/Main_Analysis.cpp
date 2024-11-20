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
    // Initialize necessary objects
    const std::string filename = "/home/mosny/Desktop/scratch/Kination_000002.3d.hdf5";
    sim_parameters params;
    read_HDF5_key_attributes(filename, params);

    // Initialize necessary objects and assign memory for pointer objects
    Vector<ProblemDomain> probDomain(params.num_levels);
    Vector<LevelData<FArrayBox> *> data(params.num_levels, NULL);
    for (int i = 0; i < 2; i++)
    {
        data[i] = new LevelData<FArrayBox>();
    }

    // Read HDF5 Data
    read_HDF5_Data(filename, probDomain, data, params);

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
