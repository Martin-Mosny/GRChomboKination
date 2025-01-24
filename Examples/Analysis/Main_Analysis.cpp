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
#include "InputOutput.hpp"
#include "computeSum.H"


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
    int max_checkpoint = 20;
    int start_checkpoint_num = 1;
    int checkpoint_interval = 1;
    int num_checkpoints = (max_checkpoint - start_checkpoint_num)/checkpoint_interval;
    std::string filename;

    // Objects we wish to calculate
    Vector<Real> time(num_checkpoints);
    Vector<Real> dx(num_checkpoints);
    Vector<Real> lapse(num_checkpoints);
    Vector<Real> E_fold(num_checkpoints);
    Vector<Real> Hubble(num_checkpoints);
    Vector<Real> Average_Phi(num_checkpoints);
    Vector<Real> Diff_Phi(num_checkpoints);
    Vector<Real> Average_Pi(num_checkpoints);
    Vector<Real> Diff_Pi(num_checkpoints);
    
    Vector<Real> Average_rho(num_checkpoints);
    Vector<Real> Diff_rho(num_checkpoints);
    Vector<Real> Rho(num_checkpoints);
    Vector<Real> Delta(num_checkpoints);

    Vector<Real> Ham(num_checkpoints);

    int c_chi = 0;
    int c_K = 7;
    int c_lapse = 18;
    int c_phi = 25;
    int c_pi = 26;
    int c_rho = 0;

    std::string filename_prefix = "/home/mosny/Desktop/scratch_Experiments/Kination4/Kination";
    for (int i = 0; i < num_checkpoints; i++)
    {
        int element = checkpoint_interval * i + start_checkpoint_num;

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
        Vector<DisjointBoxLayout> grid(params.num_levels);
        Vector<LevelData<FArrayBox> *> data(params.num_levels, NULL);
        for (int j = 0; j < params.num_levels; j++)
        {
            data[j] = new LevelData<FArrayBox>();
        }

        // Read HDF5 grid and data
        Vector<LevelData<FArrayBox>*> phi;
        Vector<LevelData<FArrayBox>*> chi;
        IntVect ghost_vect = IntVect::Unit * 3;
        read_grid(filename, grid);
        read_HDF5_ProblemDomain(filename, probDomain, params);
        read_HDF5_Data(filename, data, grid, params);

        comp_field(data, phi, c_phi);
        add_ghosts(phi, ghost_vect);
        fill_ghosts(phi, grid, probDomain, params);

        comp_field(data, chi, c_chi);
        add_ghosts(chi, ghost_vect);
        fill_ghosts(chi, grid, probDomain, params);

        LevelData<FArrayBox> *levelData1 = phi[0];
        IntVect ghost = levelData1 -> ghostVect();

        // Now we can start performing the analysis on the data
        // ----------------------------------------------------

        // Create the energy denisty object with correct layout
        Vector<LevelData<FArrayBox> *> rho(params.num_levels, NULL);
        for (int j = 0; j < params.num_levels; j++)
        {
            rho[j] = new LevelData<FArrayBox>(grid[j], 1, IntVect::Zero);
        }

        rho_calculator(rho, data, phi, params);

        // Perform all necessary single number calculations on this time step
        E_fold[i] = -0.5*std::log(average_integral(data, params, c_chi));
        Real constant_K = average_integral(data, params, c_K);
        Hubble[i] = -constant_K/3.0;
        lapse[i] = average_integral(data, params, c_lapse);
        Average_Phi[i] = average_integral(data, params, c_phi);
        Diff_Phi[i] = find_diff_value(data, c_phi);
        Average_Pi[i] = average_integral(data, params, c_pi);
        Diff_Pi[i] = find_diff_value(data, c_pi);
        
        Diff_rho[i] = find_diff_value(rho, c_rho);
        Average_rho[i] = average_integral(rho, params, c_rho);
        Rho[i] = Average_rho[i];
        Delta[i] = Diff_rho[i]/Average_rho[i];

        // Diagnostics 

        // Vector<LevelData<FArrayBox> *> Ham_array = Ham_calc(grid, data, phi, params);
        Vector<LevelData<FArrayBox> *> normalized_Ham_array = normalized_Ham_calc(grid, data, chi, params);

        int comp_val = 0;
        Real integral = average_integral(normalized_Ham_array, params, comp_val);
        Ham[i] = integral;

        pout() << "Normalized Hamiltonian constraint is: " << integral << endl;


    }


    // Here we combine the data into a form that can be printed
    // --------------------------------------------------------

    // Output file
    std::string output_filename = "/home/mosny/Desktop/GRChombo/Examples/Analysis/output_data.txt";

    // Select which data we want to print out
    Vector<Vector<Real>> output_data;
    std::vector<std::string> output_header = {"time", "dx", "lapse", "e-fold", "Hubble", "Average_Phi", "Diff_Phi", "Average_Pi", "Diff_Pi", "Rho", "Delta", "Ham"};
    output_data.push_back(time);
    output_data.push_back(dx);
    output_data.push_back(lapse);
    output_data.push_back(E_fold);
    output_data.push_back(Hubble);
    output_data.push_back(Average_Phi);
    output_data.push_back(Diff_Phi);
    output_data.push_back(Average_Pi);
    output_data.push_back(Diff_Pi);
    output_data.push_back(Rho);
    output_data.push_back(Delta);
    output_data.push_back(Ham);

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
