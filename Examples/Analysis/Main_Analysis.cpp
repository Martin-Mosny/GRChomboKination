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

    // Relevant initial values
    int num_checkpoints = 30;
    int start_checkpoint_num = 1;
    std::string filename;

    // Objects we wish to calculate
    Vector<Real> time(num_checkpoints);
    Vector<Real> dx(num_checkpoints);
    Vector<Real> lapse(num_checkpoints);
    Vector<Real> E_fold(num_checkpoints);
    Vector<Real> Hubble(num_checkpoints);
    Vector<Real> Phi(num_checkpoints);
    Vector<Real> Pi(num_checkpoints);

    std::string filename_prefix = "/home/mosny/Desktop/scratch/Kination";

    for (int i = 0; i < num_checkpoints; i++)
    {
        int element = i + start_checkpoint_num;

        if (element < 10)
        {
            filename = filename_prefix + "_00000" + std::to_string(element) + ".3d.hdf5";
        }
        else if ( element < 100)
        {
            filename = filename_prefix + "_0000" + std::to_string(element) + ".3d.hdf5";
        }
        else if (element < 1000)
        {
            filename = filename_prefix + "_000" + std::to_string(element) + ".3d.hdf5";   
        }

        // First we must read the different information from our HDF5 file
        // ---------------------------------------------------------------

        sim_parameters params;
        read_HDF5_key_attributes(filename, params);
        time[i] = params.time;
        dx[i] = params.dx[0];

        // Initialize necessary objects and assign memory for pointer objects
        Vector<ProblemDomain> probDomain(params.num_levels);
        Vector<LevelData<FArrayBox> *> data(params.num_levels, NULL);
        for (int j = 0; j < 2; j++)
        {
            data[j] = new LevelData<FArrayBox>();
        }

        // Read HDF5 Data
        read_HDF5_Data(filename, probDomain, data, params);

        // Now we can start performing the analysis on the data
        // ----------------------------------------------------

        int c_chi = 0;
        int c_K = 7;
        int c_lapse = 18;
        int c_phi = 25;
        int c_pi = 26;
        E_fold[i] = -0.5*std::log(average_integral(data, params, c_chi));
        Hubble[i] = -average_integral(data, params, c_K)/3.0;
        lapse[i] = average_integral(data, params, c_lapse);
        Phi[i] = average_integral(data, params, c_phi);
        Pi[i] = average_integral(data, params, c_pi);
    }


    // Here we combine the data into a form that can be printed
    // --------------------------------------------------------

    // Output file
    std::string output_filename = "/home/mosny/Desktop/GRChombo/Examples/Analysis/output_data.txt";

    // Select which data we want to print out
    Vector<Vector<Real>> output_data;
    std::vector<std::string> output_header = {"time", "dx", "lapse", "e-fold", "Hubble", "Phi", "Pi"};
    output_data.push_back(time);
    output_data.push_back(dx);
    output_data.push_back(lapse);
    output_data.push_back(E_fold);
    output_data.push_back(Hubble);
    output_data.push_back(Phi);
    output_data.push_back(Pi);

    // Write the output in a txt fle
    writing_to_notepad(output_filename, output_data, output_header);

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runGRChomboAnalysis(argc, argv);

    if (status == 0)
        pout() << "GRChombo analysis finished." << std::endl;
    else
        pout() << "GRChombo analysis failed with return code " << status << std::endl;

    mainFinalize();
    return status;
}
