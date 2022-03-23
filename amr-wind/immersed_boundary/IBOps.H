#ifndef IBOPS_H
#define IBOPS_H

#include "amr-wind/immersed_boundary/IBTypes.H"
#include "amr-wind/immersed_boundary/complex_terrain/ComplexTerrain.H"
#include "AMReX_Vector.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/actuator/actuator_utils.H"
#include "amr-wind/utilities/tensor_ops.H"

namespace amr_wind {
namespace ib {
namespace ops {

/** Read user inputs from an input file for a given immersed boundary object.
 *
 *  \ingroup immersed boundary
 *
 *  This operator is called within Physics::pre_init_actions() and, therefore,
 *  does not have access to mesh data. Actions that require mesh information for
 *  initialization should be deferred until ops::init_data_structures call.
 */
template <typename GeomTrait, typename = void>
struct ReadInputsOp;

/** Perform one-time initialization of the immersed boundary object.
 *
 *  \ingroup immersed boundary
 *
 *  This operator is responsible for allocating memory ensuring that the
 *  immersed boundary instance has been initialized and ready to be used by the
 * immersed boundary physics instance. This function is called once during
 * IB::post_init_actions
 *
 */
template <typename GeomTrait, typename = void>
struct InitDataOp;

/** Update the positions of the immersed boundary nodes.
 *
 *  \ingroup immersed boundary
 *
 *  After this call, the interface expects updated position vectors in the
 *  vel_pos array of the immersed boundary grid instance. These positions will
 * be used the sample the velocity field during the simulation.
 */
template <typename GeomTrait, typename = void>
struct UpdatePosOp;

/** Process updated velocities from CFD flow field.
 *
 *  \ingroup immersed boundary
 *
 *  The immersed boundary physics interface will populate the `vel` array of the
 *  immersed boundary grid instance before calling this function. The updated
 * velocities can be further processed by specializing this function.
 */
template <typename GeomTrait, typename = void>
struct UpdateVelOp;

/** Compute aerodynamic forces at the immersed boundary grid points during a
 * simulation.
 *
 *  \ingroup immersed boundary
 */
template <typename GeomTrait, typename = void>
struct ComputeForceOp;

template <typename GeomTrait, typename = void>
struct ComputeSourceOp;

template <typename IbType>
struct ComputeSourceOp<
    IbType,
    std::enable_if_t<std::is_base_of<ComplexTerrainType, IbType>::value>>
{
    void operator()(
        typename IbType::DataType& data,
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Geometry& geom)
    {
        Field& source = data.sim().repo().get_field("ib_src_term");
        const Field& level = data.sim().repo().get_field("ib_levelset");
        const Field& normal = data.sim().repo().get_field("ib_normal");
        const Field& velocity = data.sim().repo().get_field("velocity");
        const Field& density = data.sim().repo().get_field("density");

        const auto& bx = mfi.tilebox();
        const auto& varr = velocity(lev).array(mfi);
        const auto& darr = density(lev).array(mfi);
        const auto& phi_arr = level(lev).array(mfi);
        const auto& norm_arr = normal(lev).array(mfi);
        auto sarr = source(lev).array(mfi);

        const amrex::Real cd = data.meta().drag_coefficient;
        const amrex::Real eps = data.meta().eps_spread;

        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // A lot of this is unnecessary for a flat plate, but is
                // general enough that it could be broken off into its own
                // function at some point.
                const amrex::Vector<amrex::Real> vel = {
                    {varr(i, j, k, 0), varr(i, j, k, 1), varr(i, j, k, 2)}};
                const amrex::Real vmag = utils::vec_mag(vel.data());
                const amrex::Vector<amrex::Real> norm = {
                    {norm_arr(i, j, k, 0), norm_arr(i, j, k, 1),
                     norm_arr(i, j, k, 2)}};
                // Outside the band if the normal vector is assigned to be zero
                const amrex::Real nmag = utils::vec_mag(norm.data());
                if (nmag > 0.5) {
                    const amrex::Real vdotn =
                        utils::dot_prod(norm.data(), vel.data());
                    const auto gauss_fac = actuator::utils::gaussian1d(
                        amrex::Math::abs(phi_arr(i, j, k)), eps);
                    // v_tan = v - (v dot n) n
                    amrex::Vector<amrex::Real> vtan = {
                        {vel[0] - vdotn * norm[0], vel[1] - vdotn * norm[1],
                         vel[2] - vdotn * norm[2]}};
                    // Get tangent velocity magnitude
                    const amrex::Real vtmag = utils::vec_mag(vtan.data());
                    // Calculate drag force C_D * 0.5 * rho * vtmag^2
                    // vtan already has one factor of vtmag in it, so there's no
                    // need to normalize, just leave off a factor in the
                    // coefficient.
                    const amrex::Real force =
                        -0.5 * darr(i, j, k) * gauss_fac * vtmag * cd;
                    sarr(i, j, k, 0) += force * vtan[0];
                    sarr(i, j, k, 1) += force * vtan[1];
                    sarr(i, j, k, 2) += force * vtan[2];
                }
            });
    }
};

template <typename GeomTrait, typename = void>
struct ProcessOutputsOp;

} // namespace ops
} // namespace ib
} // namespace amr_wind

#endif /* IBOPS_H */
