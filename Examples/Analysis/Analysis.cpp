// Here we define the various functions used to analyse the data

#include "Analysis.hpp"
#include <iostream>
#include "LoadBalance.H"
#include "GRAMRLevel.hpp"
#include "computeSum.H"

void create_layout_data(Vector<LevelData<FArrayBox> *> &a_output, Vector<LevelData<FArrayBox> *> &a_input, Real &val)
{
    // This function creates an a_output object of the same layout dimension as a_input 
    // except only with one component, where every value is val.

    // First we create the correct layout for a_output
    for (int i = 0; i < a_input.size(); i++)
    {
        const DisjointBoxLayout &layout = a_input[i] -> disjointBoxLayout();
        a_output[i] = new LevelData<FArrayBox>(layout, 1, IntVect::Zero);
    }

    // Now we loop through all its elements and set them to 1
    for (int i = 0; i < a_output.size(); i++)
    {
        const DisjointBoxLayout &layout = a_output[i] -> disjointBoxLayout();
        DataIterator dit = layout.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
        {
            FArrayBox& fab = (*a_output[i])[dit()];
            for (BoxIterator bit(fab.box()); bit.ok(); ++bit)
            {
                const IntVect& iv = bit();
                fab(iv, 0) = val;
            }
        }
    }
    
}

Real average_integral(Vector<LevelData<FArrayBox> *> &a_data, sim_parameters &a_params, int &comp)
{
    // This function returns the average value of a_data over our ProblemDomain
    Real volume = pow(a_params.dx[0]*32, 3);
    Real integral = computeSum(a_data, a_params.ref_ratio, a_params.dx[0], Interval(comp, comp));
    return integral/volume;
}

Real find_max_value(Vector<LevelData<FArrayBox> *> &a_data, int &comp)
{
    // This function finds the maximal value of a parameter
    std::vector<Real> max_elements;

    // Now we loop through all its FArrayBoxes
    for (int i = 0; i < a_data.size(); i++)
    {
        const DisjointBoxLayout &layout = a_data[i] -> disjointBoxLayout();
        DataIterator dit = layout.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
        {
            Real max_FArray_comp = (*a_data[i])[dit()].max(comp);
            max_elements.push_back(max_FArray_comp);
        }
    }

    Real max_val = *std::max_element(max_elements.begin(), max_elements.end());
    return max_val;
    
}

Real find_diff_value(Vector<LevelData<FArrayBox> *> &a_data, int &comp)
{
    // This function finds the maximal value of a parameter
    std::vector<Real> max_elements;
    std::vector<Real> min_elements;

    // Now we loop through all its FArrayBoxes
    for (int i = 0; i < a_data.size(); i++)
    {
        const DisjointBoxLayout &layout = a_data[i] -> disjointBoxLayout();
        DataIterator dit = layout.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
        {
            Real max_FArray_comp = (*a_data[i])[dit()].max(comp);
            Real min_FArray_comp = (*a_data[i])[dit()].min(comp);
            max_elements.push_back(max_FArray_comp);
            min_elements.push_back(min_FArray_comp);

        }
    }

    Real max_val = *std::max_element(max_elements.begin(), max_elements.end());
    Real min_val = *std::min_element(min_elements.begin(), min_elements.end());
    Real diff = max_val - min_val;
    return diff;
    
}

inline Real get_grad_phi_sq(const IntVect &a_iv, const FArrayBox &a_phi_fab,
                            const Real &a_dx, const FArrayBox &a_data)
{
    Real grad_phi_sq = 0.0;
    int c_index;
    for (int idir_1 = 0; idir_1 < SpaceDim; ++idir_1)
    {
        IntVect iv_offset11 = a_iv;
        IntVect iv_offset12 = a_iv;
        iv_offset11[idir_1] -= 1;
        iv_offset12[idir_1] += 1;
    
        // 2nd order stencils for now
        Real dphi_dx_1 =
            0.5 * (a_phi_fab(iv_offset12) - a_phi_fab(iv_offset12)) / a_dx;

        for (int idir_2 = 0; idir_2 > SpaceDim; ++idir_2)
        {
            IntVect iv_offset21 = a_iv;
            IntVect iv_offset22 = a_iv;
            iv_offset21[idir_2] -= 1;
            iv_offset22[idir_2] += 1;

            // 2nd order stencils for now
            Real dphi_dx_2 =
                0.5 * (a_phi_fab(iv_offset21) - a_phi_fab(iv_offset22)) / a_dx;

            if (idir_1 != 0 && idir_2 != 0)
            {
                c_index = idir_1 + idir_2 + 2;
            }
            else if (idir_1 = 0)
            {
                c_index = idir_2 + 1;
            }
            else
            {
                c_index = idir_1 + 1;
            }

            grad_phi_sq += a_data(a_iv, c_index) * dphi_dx_1 * dphi_dx_2;
        }
    }
    return grad_phi_sq;
}

void rho_calculator(Vector<LevelData<FArrayBox> *> &rho, Vector<LevelData<FArrayBox> *> &a_data, Vector<LevelData<FArrayBox> *> &a_phi, sim_parameters &a_param)
{
    // Calculates the energy density
    for (int i = 0; i < rho.size(); i++)
    {
        const DisjointBoxLayout &layout = rho[i] -> disjointBoxLayout();
        DataIterator dit = layout.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
        {
            FArrayBox& f_rho = (*rho[i])[dit()];
            FArrayBox& f_data = (*a_data[i])[dit()];
            FArrayBox& f_phi = (*a_phi[i])[dit()];
            for (BoxIterator bit(f_rho.box()); bit.ok(); ++bit)
            {
                const IntVect& iv = bit();
                f_rho(iv, 0) = pow(f_data(iv, 26), 2) / 2.0 + f_data(iv, 0) * get_grad_phi_sq(iv, f_phi, a_param.dx[i], f_data) / 2.0;
            }
        }
    }
}

// computes the Laplacian of psi at a point in a box
inline Real get_laplacian_psi_from_chi(const IntVect &a_iv, const FArrayBox &a_psi_fab,
                              const Real &a_dx)
{
    Real laplacian_of_psi = 0.0;
    for (int idir = 0; idir < SpaceDim; ++idir)
    {
        IntVect iv_offset1 = a_iv;
        IntVect iv_offset2 = a_iv;
        iv_offset1[idir] -= 1;
        iv_offset2[idir] += 1;

        // 2nd order stencil for now
        Real d2psi_dxdx = 1.0 / (a_dx * a_dx) *
                          (+1.0 * pow(a_psi_fab(iv_offset2), -0.25) -
                           2.0 * pow(a_psi_fab(a_iv), -0.25) + 1.0 * pow(a_psi_fab(iv_offset1), -0.25));
        laplacian_of_psi += d2psi_dxdx;
    }
    return laplacian_of_psi;
} // end get_laplacian_psi

Vector<LevelData<FArrayBox> *> Ham_calc(const Vector<DisjointBoxLayout> &grids, Vector<LevelData<FArrayBox> *> a_data, Vector<LevelData<FArrayBox> *> &a_phi, const sim_parameters &a_params)
{   
    Vector<LevelData<FArrayBox> *> Ham(a_data.size(), NULL);

    // Initialize the Hamiltonian constraint girds
    for (int ilev = 0; ilev < Ham.size(); ilev++)
    {
        Ham[ilev] = new LevelData<FArrayBox>(grids[ilev], 1, IntVect::Zero);
    }

    // Fill in the Hamiltonian constraint
    for (int ilev = 0; ilev < Ham.size(); ilev++)
        {
            LevelData<FArrayBox> *Fields = a_data[ilev];
            LevelData<FArrayBox> *Ham_fields = Ham[ilev];
            LevelData<FArrayBox> *phi_fields = a_phi[ilev];

            DataIterator dit = grids[ilev].dataIterator();
            for (dit.begin(); dit.ok(); ++dit)
            {
                FArrayBox& Fields_box = (*Fields)[dit()];
                FArrayBox& Ham_box = (*Ham_fields)[dit()];
                FArrayBox& psi_fab = (*phi_fields)[dit()];
                const Box& b = grids[ilev][dit];
                BoxIterator bit(b);
                for (bit.begin(); bit.ok(); ++bit)
                {
                    // Must take the modular coordinate to assign ghost values correctly
                    IntVect iv = bit();
                    Real laplacian_value = get_laplacian_psi_from_chi(iv, psi_fab, a_params.dx[ilev]);
                    Ham_box(iv) = laplacian_value 
                       - (pow(psi_fab(iv), -1.25) * pow(Fields_box(iv, 7), 2) / 12.0 
                       - M_PI * pow(Fields_box(iv, 26), 2) * pow(psi_fab(iv), -1.25));
                }
            }
        }
    return Ham;
}

Vector<LevelData<FArrayBox> *> normalized_Ham_calc(const Vector<DisjointBoxLayout> &grids, Vector<LevelData<FArrayBox> *> a_data, Vector<LevelData<FArrayBox> *> &a_phi, const sim_parameters &a_params)
{   
    Vector<LevelData<FArrayBox> *> normalized_Ham(a_data.size(), NULL);

    // Initialize the Hamiltonian constraint girds
    for (int ilev = 0; ilev < normalized_Ham.size(); ilev++)
    {
        normalized_Ham[ilev] = new LevelData<FArrayBox>(grids[ilev], 1, IntVect::Zero);
    }

    // Fill in the Hamiltonian constraint
    for (int ilev = 0; ilev < normalized_Ham.size(); ilev++)
        {
            LevelData<FArrayBox> *Fields = a_data[ilev];
            LevelData<FArrayBox> *normalized_Ham_fields = normalized_Ham[ilev];
            LevelData<FArrayBox> *phi_fields = a_phi[ilev];

            DataIterator dit = grids[ilev].dataIterator();
            for (dit.begin(); dit.ok(); ++dit)
            {
                FArrayBox& Fields_box = (*Fields)[dit()];
                FArrayBox& normalized_Ham_box = (*normalized_Ham_fields)[dit()];
                FArrayBox& psi_fab = (*phi_fields)[dit()];
                const Box& b = grids[ilev][dit];
                BoxIterator bit(b);
                for (bit.begin(); bit.ok(); ++bit)
                {
                    // Must take the modular coordinate to assign ghost values correctly
                    IntVect iv = bit();
                    Real laplacian_value = get_laplacian_psi_from_chi(iv, psi_fab, a_params.dx[ilev]);
                    Real Ham_value = laplacian_value 
                       - (pow(psi_fab(iv), -1.25) * pow(Fields_box(iv, 7), 2) / 12.0 
                       - M_PI * pow(Fields_box(iv, 26), 2) * pow(psi_fab(iv), -1.25));
                    Real Ham_ref = sqrt(pow(laplacian_value, 2) 
                       + pow((pow(psi_fab(iv), -1.25) * pow(Fields_box(iv, 7), 2) / 12.0), 2)
                       + pow(M_PI * pow(Fields_box(iv, 26), 2) * pow(psi_fab(iv), -1.25), 2));
                    
                    // pout() << "Lapacian: " << laplacian_value << endl;
                    // pout() << "Term 1:   " << pow(psi_fab(iv), -1.25) * pow(Fields_box(iv, 7), 2) / 12.0  << endl;
                    // pout() << "Psi Term: " << psi_fab(iv) << endl;
                    // pout() << "Term 2:   " << M_PI * pow(Fields_box(iv, 26), 2) * pow(psi_fab(iv), -1.25) << endl;
                    // pout() << "Ref value: " << Ham_ref << endl;
                    // pout() << "Act value: " << Ham_value << endl;

                    normalized_Ham_box(iv) = Ham_value/Ham_ref; 

                }
            }
        }
    return normalized_Ham;
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

/*

bool contains_vector(const DisjointBoxLayout& finer_grids, const IntVector& iv)
{
    // This function checks whether a vector iv is contained in a DisjointBoxLayout
    for (LayoutIterator it = finerGrids.layoutIterator(); it.ok(); ++it)
    {
        const Box& b = finerGrids[it];
        if (box.contains(iv))
        {
            return true;
        }
    }
    return false;
}

Real refined_average(const Vector<DisjointBoxLayout> &grids, Vector<LevelData<FArrayBox> *> &data, int comp, sim_params &a_params)
{
    Real sum = 0.0;
    // Iterate over levels starting at the most refined
    for (int ilev = data.size(); ilev >= 0; --ilev)
    {
        LevelData<FArrayBox> *Fields = data[ilev];
        const DisjointBoxLayout& levelGrids = grids[ilev];
        DataIterator dit = grids[ilev].dataIterator();
        // Iterate over the boxes on a particular level
        for (dit.begin(); dit.ok(); ++dit())
        {
            FArrayBox& Fields_box = (*Fields)[dit()];
            FArrayBox Field_box(Interval(comp, comp), Fields_box);
            const Box& b = grids[ilev][dit];
            // Iterate over all vectors in a paricular box
            for (bit.begin(); bit.ok(); ++bit)
            {
                IntVect iv = bit();
                bool covered = false;

                // No need to check the most refined level
                if (ilev != data.size())
                {
                    // Check that iv is not contained in the more refined level
                    DisjointBoxLayout &finer_grids = grids[ilev + 1];
                    covered = contains_vector(finer_grids, iv);
                }
                
                // If not covered, add the weighted value to the sum
                if (!covered) 
                {
                    Real value = Field_box(iv);
                    sum += value / pow(8.0, ilev);
                }

                
            }

        }

    }

    return sum;
}

*/