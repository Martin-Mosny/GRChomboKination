// Here we define the various functions used to analyse the data

#include "InputOutput.hpp"
#include <iostream>
#include "LoadBalance.H"
#include "GRAMRLevel.hpp"
#include "computeSum.H"
#include "QuadCFInterp.H"


struct ProbDomVector
{
    // Datatype for reading HDF5 ProblemDomain vectors
    int lo_i;
    int lo_j;
    int lo_k;
    int hi_i;
    int hi_j;
    int hi_k;
};

void read_HDF5_key_attributes(const std::string &filename, sim_parameters &a_params)
{
    // Load in the HDF5 file
    H5File file(filename, H5F_ACC_RDONLY);

    // Read the max level
    Group group = file.openGroup("/");
    Attribute attribute = group.openAttribute("max_level");
    DataType type = attribute.getDataType();
    attribute.read(type, &a_params.max_level);

    // Read the num_levels
    attribute = group.openAttribute("num_levels");
    type = attribute.getDataType();
    attribute.read(type, &a_params.num_levels);

    // Read the time
    attribute = group.openAttribute("time");
    type = attribute.getDataType();
    attribute.read(type, &a_params.time);

    // Read the num_components
    attribute = group.openAttribute("num_components");
    type = attribute.getDataType();
    attribute.read(type, &a_params.num_components);

    // Resize the dx and ref_ratio vectors
    Real dx_val;
    int ref_ratio_val;
    a_params.dx.resize(a_params.num_levels, 0);
    a_params.ref_ratio.resize(a_params.num_levels, 0);

    // Read the dx and ref_ratios
    for (int i = 0; i < a_params.num_levels; i++)
    {
        std::string groupName = "/level_" + std::to_string(i);
        group = file.openGroup(groupName);
        
        // First the dx
        attribute = group.openAttribute("dx");
        type = attribute.getDataType();
        attribute.read(type, &dx_val);
        a_params.dx[i] = dx_val;
        
        // Now for the ref_ratios
        attribute = group.openAttribute("ref_ratio");
        type = attribute.getDataType();
        attribute.read(type, &ref_ratio_val);
        a_params.ref_ratio[i] = ref_ratio_val;
    }

}

void read_HDF5_LevelProblemDomain(const std::string &filename, ProblemDomain &probDomain, int &level_int)
{
     
    // Load in the problem domain box size attributes
    const std::string groupName = "/level_" + std::to_string(level_int);
    H5File file(filename, H5F_ACC_RDONLY);
    Group group = file.openGroup(groupName);
    Attribute attribute = group.openAttribute("prob_domain");

    // Create the necessary data type to read in the high and low vectors of PorblemDomain
    CompType HighLowType(sizeof(ProbDomVector));
    HighLowType.insertMember("lo_i", HOFFSET(ProbDomVector, lo_i), PredType::NATIVE_INT);
    HighLowType.insertMember("lo_j", HOFFSET(ProbDomVector, lo_j), PredType::NATIVE_INT);
    HighLowType.insertMember("lo_k", HOFFSET(ProbDomVector, lo_k), PredType::NATIVE_INT);
    HighLowType.insertMember("hi_i", HOFFSET(ProbDomVector, hi_i), PredType::NATIVE_INT);
    HighLowType.insertMember("hi_j", HOFFSET(ProbDomVector, hi_j), PredType::NATIVE_INT);
    HighLowType.insertMember("hi_k", HOFFSET(ProbDomVector, hi_k), PredType::NATIVE_INT);


    // Read the data into a PDV oject of composite data type
    // Then split it up into two IntVectors
    ProbDomVector PDV;
    attribute.read(HighLowType, &PDV);
    IntVect lower(PDV.lo_i, PDV.lo_j, PDV.lo_k);
    IntVect upper(PDV.hi_i, PDV.hi_j, PDV.hi_k);
    
    // Set the problem domain periodicity
    bool periodicity_array[3] = {0,0,0};
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

void read_HDF5_LevelData(const std::string &filename, ProblemDomain &probDomain, LevelData<FArrayBox> &a_data, int &level_int)
{
    // Read in the problem domain
    read_HDF5_LevelProblemDomain(filename, probDomain, level_int);

    // Open the HDF5 file handle
    HDF5Handle a_handle(filename, HDF5Handle::OPEN_RDONLY);
    
    // Read in the boxes
    Vector<Box> boxes;
    std::string level = "/level_0";// + std::to_string(level_int);
    std::string a_var1 = level + "/boxes";
    std::string &a_name1 = a_var1;
    read(a_handle, boxes, a_name1);

    // Assign processors using LoadBalance and create DisjointBoxLayout
    Vector<int> procAssign(boxes.size());
    LoadBalance(procAssign, boxes);
    DisjointBoxLayout a_layout(boxes, procAssign, probDomain);

    // Read in the data
    std::string a_variable = level + "/data";
    std::string &a_name = a_variable;
    read(a_handle, a_data, a_name, a_layout);

    a_handle.close();
    
}

void read_grid(const std::string &filename, Vector<DisjointBoxLayout> &a_grid)
{
    // This function reads in the grid format from the HDF5 file

    // First open the file
    HDF5Handle a_handle(filename, HDF5Handle::OPEN_RDONLY);

    for (int i = 0; i < a_grid.size(); i++)
    {
    
    // Read the problem domain
    ProblemDomain probDomain;
    read_HDF5_LevelProblemDomain(filename, probDomain, i);

    // Read in the boxes
    Vector<Box> boxes;
    std::string level = "/level_" + std::to_string(i);
    std::string a_var1 = level + "/boxes";
    std::string &a_name1 = a_var1;
    read(a_handle, boxes, a_name1);

    // Assign processors and define DisjointBoxLayout
    Vector<int> procAssign(boxes.size());
    LoadBalance(procAssign, boxes);
    DisjointBoxLayout a_layout(boxes, procAssign, probDomain);
    a_grid[i] = a_layout;
    }

    a_handle.close();
}


void read_HDF5_ProblemDomain(const std::string &filename, Vector<ProblemDomain> &probDomain, sim_parameters &a_params)
{
    // Run over the various levels
    ProblemDomain dom;
    for (int i = 0; i < a_params.num_levels; i++)
    {
        read_HDF5_LevelProblemDomain(filename, dom, i);
        probDomain[i] = dom;
    }
}


void read_HDF5_Data(const std::string &filename, Vector<LevelData<FArrayBox> *> &a_data, Vector<DisjointBoxLayout> &a_grid, sim_parameters &a_params)
{
    // This simply reads the HDF5 data over all levels

    // Open the HDF5 file handle and loop over levels
    HDF5Handle a_handle(filename, HDF5Handle::OPEN_RDONLY);
    for (int i = 0; i < a_params.num_levels; i++)
    {
        // Read in the data
        std::string level = "/level_" + std::to_string(i);
        std::string a_variable = level + "/data";
        std::string &a_name = a_variable;
        read(a_handle, *a_data[i], a_name, a_grid[i]);
    }

    a_handle.close();
}

void add_ghosts(Vector<LevelData<FArrayBox> *> &a_data, const IntVect &ghost_vect)
{
    // This function adds ghosts to a ghostless data set
    for (int level = 0; level < a_data.size(); ++level)
    {
        // Create a copy of the data and layout
        LevelData<FArrayBox> *old_data = a_data[level];
        const DisjointBoxLayout &old_layout = old_data -> disjointBoxLayout();
        int num_components = old_data -> nComp();

        // Create a new LevelData but with ghosts
        LevelData<FArrayBox> *new_data = new LevelData<FArrayBox>(old_layout, num_components, ghost_vect);

        // Copy the data into the new LevelData
        old_data -> copyTo(new_data -> interval(), *new_data, new_data -> interval());

        // Replace the old LevelData with the new LevelData in the same object
        delete old_data;
        a_data[level] = new_data;
    }
}


void fill_ghosts(Vector<LevelData<FArrayBox> *> &a_data, Vector<DisjointBoxLayout> &a_grid, Vector<ProblemDomain> &probDomain, sim_parameters &a_params)
{
    // This fills in the ghosts for a 1-component object; needs to be modified in QuadCFInterp 1 -> nComp for more
    for (int level = 1; level < a_data.size(); level++)
    {
        QuadCFInterp quadCFI(a_grid[level], &a_grid[level - 1], a_params.dx[level], a_params.ref_ratio[level], 1, probDomain[level]);
        quadCFI.coarseFineInterp(*a_data[level], *a_data[level - 1]);
    }
}


void comp_field(Vector<LevelData<FArrayBox> *> &a_data, Vector<LevelData<FArrayBox> *> &a_output, int component)
{
    // This is a function that creates an a_output that only consists of the component of a_data
    a_output.clear();

    for (int level = 0; level < a_data.size(); ++level)
    {
        // Create a copy of the data and layout
        const LevelData<FArrayBox> *inputLevel = a_data[level];

        // Create a new LevelData
        LevelData<FArrayBox> *outputLevel = new LevelData<FArrayBox>(inputLevel -> disjointBoxLayout(), 1, inputLevel -> ghostVect());

        for (DataIterator dit = inputLevel -> dataIterator(); dit.ok(); ++dit)
        {
            const FArrayBox &inputFab = (*inputLevel)[dit];
            FArrayBox &outputFab = (*outputLevel)[dit];

            outputFab.copy(inputFab, component, 0, 1);
        }

        a_output.push_back(outputLevel);
    }
}


void writing_to_notepad(const std::string filename, const Vector<Vector<Real>> &a_data, const std::vector<std::string> &headers)
{
    // This function writes out our data in the format row 1: header[0],header[1],... 
    // and then each following row i is the next element data[0][i], data[1][i], ...

    // Write our the header names
    std::ofstream outputFile(filename);
    for (int i = 0; i < headers.size(); i++)
    {
        outputFile << headers[i];
        if (i < headers.size() - 1)
        {
            outputFile << ",";
        }

    }
    outputFile << "\n";

    // Write out the data
    for (int j = 0; j < a_data[0].size(); ++j)
    {
        for (int i = 0; i < a_data.size(); ++i)
        {
            outputFile << a_data[i][j];
            if (i <  a_data.size() - 1)
            {
                outputFile << ",";
            }
        }

        outputFile << "\n";
    }

}


