// Here we define the various functions used to analyse the data

#include "Analysis.hpp"
#include <iostream>
#include "LoadBalance.H"
#include "GRAMRLevel.hpp"

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

    // Read the dx and ref_ratios
    
    Real dx_val;
    int ref_ratio_val;
    for (int i = 0; i < a_params.num_levels; i++)
    {
        std::string groupName = "/level_" + std::to_string(i);
        group = file.openGroup(groupName);
        
        // First the dx
        attribute = group.openAttribute("dx");
        type = attribute.getDataType();
        attribute.read(type, &dx_val);
        
        a_params.dx.append(dx_val);
        
        // Now for the ref_ratios
        attribute = group.openAttribute("ref_ratio");
        type = attribute.getDataType();
        attribute.read(type, &ref_ratio_val);
        a_params.ref_ratio.append(ref_ratio_val);
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
    pout() << "Problem domain low vector is: (" << PDV.lo_i << "," << PDV.lo_j << "," << PDV.lo_k << ")" << endl;
    pout() << "Problem domain high vector is: (" << PDV.hi_i << "," << PDV.hi_j << "," << PDV.hi_k << ")" << endl;
    
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
    for (int i = 0; i < 8; i++)
    {
        pout() << boxes[i].smallEnd()[0] << " " << boxes[i].smallEnd()[1] << " " << boxes[i].smallEnd()[2] << endl;
    }
    
    for (int i = 0; i < 8; i++)
    {
        pout() << boxes[i].bigEnd()[0] << " " << boxes[i].bigEnd()[1] << " " << boxes[i].bigEnd()[2] << endl;
    }
    //pout() << "Box " << i << " intersects: " << probDomain.intersects(boxes[i], boxes[i+1]) << endl;;

    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            pout() << "Boxes " << i << " and " << j << " intersect: " <<probDomain.intersects(boxes[i], boxes[j]) << endl;
        }
    }

    // Assign processors using LoadBalance and create DisjointBoxLayout
    Vector<int> procAssign(boxes.size());
    LoadBalance(procAssign, boxes);
    DisjointBoxLayout a_layout(boxes, procAssign, probDomain);

    // Read in the data
    std::string a_variable = level + "/data";
    std::string &a_name = a_variable;
    pout() << "HELLO" << endl;
    read(a_handle, a_data, a_name, a_layout);
    pout() << "DID THIS WORK" <<endl;
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


void read_HDF5_Data(const std::string &filename, Vector<ProblemDomain> &probDomain, Vector<LevelData<FArrayBox> *> &a_data, sim_parameters &a_params)
{
    for (int i = 0; i < 2; i++)
    {
        read_HDF5_LevelData(filename, probDomain[i], *a_data[i], i);
    }
    
}






void scale_factor(HDF5Handle & a_handle, BHAMR bh_amr)
{
    // Calculate the scale factor for a particular slice

    //Calculate the volume integrand


    //Real integral = computeSum();
}


    /*
    // Example code of how you iterate over our data = object
    const DisjointBoxLayout &layout = object.disjointBoxLayout();
    DataIterator dit = layout.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox& fab = object[dit()];
        for (BoxIterator bit(fab.box()); bit.ok(); ++bit)
        {
            const IntVect& iv = bit();
            pout() << fab(iv, 1) << endl;
        }
    }
    */