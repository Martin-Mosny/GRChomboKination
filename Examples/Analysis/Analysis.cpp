// Here we define the various functions used to analyse the data

#include "Analysis.hpp"
#include <iostream>
#include "LoadBalance.H"
#include "GRAMRLevel.hpp"

struct ProbDomVector
{
    int small[3];
    int big[3];
};

void read_ProblemDomain(const std::string filename, ProblemDomain probDomain)
{
     
    // Load in the problem domain box size attributes
    const std::string groupName = "/level_0";
    const std::string attributeName = "prob_domain";
    H5File file(filename, H5F_ACC_RDONLY);
    Group group = file.openGroup(groupName);
    Attribute attribute = group.openAttribute(attributeName);
 
    // Create a composite type to read the hdf5 data
    CompType HighLowType(sizeof(ProbDomVector));
    hsize_t dims[1] = {3};
    ArrayType smallType(PredType::NATIVE_INT, 1, dims);
    HighLowType.insertMember("small", HOFFSET(ProbDomVector, small), smallType);
    HighLowType.insertMember("big", HOFFSET(ProbDomVector, big), smallType);
   
    // Read the data into a PDV oject of composite data type
    // Then split it up into two IntVectors
    ProbDomVector PDV;
    attribute.read(HighLowType, &PDV);
    IntVect lower(PDV.small[0], PDV.small[1], PDV.small[2]);
    IntVect upper(PDV.big[0], PDV.big[1], PDV.big[2]);


    bool periodicity_array[3] = {0,0,0};
    // Set the problem domain periodicity
    for (int i = 0; i < 3; i++)
    {
        const std::string periodicity_att_name = "is_periodic_" + std::to_string(i);
        Attribute periodicity_att = group.openAttribute(periodicity_att_name);
        DataType type = periodicity_att.getDataType();
        int periodicity;
        periodicity_att.read(type, &periodicity);
        periodicity_array[i] = periodicity;
    }

    // Define the problem domain
    probDomain.define(lower, upper, periodicity_array);
}













void scale_factor(HDF5Handle & a_handle, BHAMR bh_amr)
{
    // Calculate the scale factor for a particular slice

    //Calculate the volume integrand


    //Real integral = computeSum();
}

/*

void readHDF5(HDF5Handle &a_handle, int &num_levels, Vector<DisjointBoxLayout> &all_grids, Vector<LevelData<FArrayBox>> &all_states)
{
    // Read the number of levels from the HDF5 header
    HDF5HeaderData header;
    header.readFromFile(a_handle);

    if (header.m_int.find("num_levels") == header.m_int.end())
    {
        MayDay::Error("Checkpoint file missing num_levels information.");
    }

    num_levels = header.m_int["num_levels"];
    std::cout << "Number of levels: " << num_levels << std::endl;

    all_grids.resize(num_levels);
    all_states.resize(num_levels);

    for (int level = 0; level < num_levels; ++level)
    {
        // Create the level group name
        char level_str[20];
        sprintf(level_str, "level_%d", level);
        a_handle.setGroup(level_str);

        // Read the grids for this level
        Vector<Box> grids;
        if (read(a_handle, grids) != 0)
        {
            MayDay::Error("Error: Could not read grids from HDF5 file.");
        }

        // Load grids into DisjointBoxLayout
        all_grids[level] = loadBalance(grids);

        // Prepare storage for state data
        all_states[level].define(all_grids[level], NUM_VARS, IntVect::Zero);

        // Read the state data
        Interval comps(0, NUM_VARS - 1);
        bool redefine_data = false;
        if (read<FArrayBox>(a_handle, all_states[level], "data", all_grids[level], comps, redefine_data) != 0)
        {
            MayDay::Error("Error: Could not read state data from HDF5 file.");
        }

        std::cout << "Level " << level << ": " << grids.size() << " grids read." << std::endl;
    }
}
