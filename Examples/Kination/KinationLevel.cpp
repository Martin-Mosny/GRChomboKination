/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "KinationLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "MatterCCZ4RHS.hpp"

// For constraints calculation
#include "NewMatterConstraints.hpp"

// For tag cells
#include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "GammaCalculator.hpp"
#include "InitialScalarData.hpp"
#include "KerrBH.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "SetValue.hpp"
#include "computeSum.H"

// Things to do at each advance step, after the RK4 is calculated
void KinationLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(
            NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance"),
            m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void KinationLevel::initialData()
{
    CH_TIME("KinationLevel::initialData");
    if (m_verbosity)
        pout() << "KinationLevel::initialData " << m_level << endl;

    // First set everything to zero then initial conditions for scalar field -
    // here a Kerr BH and a scalar field profile
    BoxLoops::loop(
        make_compute_pack(SetValue(0.), KerrBH(m_p.kerr_params, m_dx),
                          InitialScalarData(m_p.initial_params, m_dx)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);
}

#ifdef CH_USE_HDF5
// Things to do before outputting a checkpoint file
void KinationLevel::prePlotLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(
        MatterConstraints<ScalarFieldWithPotential>(
            scalar_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom, c_Mom)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    pout() << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    Vector<int> ref_rat(2, m_ref_ratio);
    DisjointBoxLayout layout = m_state_diagnostics.disjointBoxLayout();
    const DisjointBoxLayout* layout_ptr = &layout;
    Real integral = computeSum(m_state_diagnostics, &layout, m_ref_ratio, m_dx, Interval(c_Ham,c_Ham));
    Real volume = pow(1.5625*32, 3);
    pout() << "Hamiltonian constraint is " << integral/volume << endl;
}
#endif

// Things to do in RHS update, at each RK4 step
void KinationLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    if (m_p.max_spatial_derivative_order == 4)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      FourthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      SixthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// Things to do at ODE update, after soln + rhs
void KinationLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void KinationLevel::preTagCells()
{
    // we don't need any ghosts filled for the fixed grids tagging criterion
    // used here so don't fill any
}

void KinationLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)
{
    BoxLoops::loop(
        FixedGridsTaggingCriterion(m_dx, m_level, 2.0 * m_p.L, m_p.center),
        current_state, tagging_criterion);
}
void KinationLevel::specificPostTimeStep()
{
#ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
#endif
}
/*
// computes the gradient of the scalar field squared at a point in a box
// i.e. \delta^{ij} d_i phi d_j phi
inline Real get_grad_phi_sq(const IntVect &a_iv, const FArrayBox &a_phi_fab,
                            const RealVect &a_dx)
{
    Real grad_phi_sq = 0.0;
    for (int idir = 0; idir < SpaceDim; ++idir)
    {
        IntVect iv_offset1 = a_iv;
        IntVect iv_offset2 = a_iv;
        iv_offset1[idir] -= 1;
        iv_offset2[idir] += 1;

        // 2nd order stencils for now
        Real dphi_dx =
            0.5 * (a_phi_fab(iv_offset2) - a_phi_fab(iv_offset1)) / a_dx[idir];

        grad_phi_sq += dphi_dx * dphi_dx;
    }
    return grad_phi_sq;
} // end get_grad_phi_sq

void Hamiltonian_constraint_level(DisjointBoxLayout &grids, LevelData<FArrayBox> multigrid_vars, const RealVect &a_dx,
             const PoissonParameters &a_params, const Real constant_K)
{   
    Real Ham_constraint = 0;
    int Num = 0;
    Real volume = pow(1.5625*32, 3);
    LevelData<FArrayBox> *Fields = multigrid_vars[ilev];
    //LevelData<FArrayBox> *rhs_level_data = rhs[ilev];
    //DataIterator dit = Fields -> dataIterator();
    DataIterator dit = grids[ilev].dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
        {
        FArrayBox& Fields_box = (*Fields)[dit()];
        // FArrayBox& rhs_box = (*rhs_level_data)[dit()];
        FArrayBox psi_fab(Interval(c_psi_reg, c_psi_reg), Fields_box);
        const Box& b = grids[ilev][dit];
        BoxIterator bit(b);
        for (bit.begin(); bit.ok(); ++bit)
        {
            // Must take the modular coordinate to assign ghost values correctly
            IntVect iv = bit();
            Real laplacian_value = get_laplacian_psi(iv, psi_fab, a_dx[ilev]);
            Num += 1;
            Ham_constraint += laplacian_value 
               - (pow(psi_fab(iv), 5) * pow(constant_K, 2) / 12.0 
                - M_PI * pow(Fields_box(iv, c_Pi_0), 2) * pow(psi_fab(iv), 5));
        }
    }

    pout() << "The Hamiltonian constraint at level is give by " << Ham_constraint / Num << endl;
} // end set_rhs
*/
