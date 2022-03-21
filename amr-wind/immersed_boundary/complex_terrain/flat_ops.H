#ifndef FLAT_OPS_H
#define FLAT_OPS_H

#include "amr-wind/immersed_boundary/complex_terrain/Flat.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/immersed_boundary/IB.H"
#include "amr-wind/immersed_boundary/complex_terrain/complex_terrain_ops.H"
// nmatula Probably can use a local version of the gaussian function in the long
// run
#include "amr-wind/wind_energy/actuator/actuator_utils.H"

namespace amr_wind {
namespace ib {

namespace ops {

template <>
struct ReadInputsOp<Flat>
{
    void
    operator()(Flat::DataType& data, const ::amr_wind::utils::MultiParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();

        complex_terrain::read_inputs(wdata, info, pp);
        pp.get("height", wdata.height);
        pp.get("drag_coefficient", wdata.drag_coefficient);
        pp.get("eps_spread", wdata.eps_spread);
    }
};

template <>
struct InitDataOp<Flat>
{
    //! Immersed boundary force term for wall modeled terrain

    void operator()(Flat::DataType& data)
    {
        const auto& wdata = data.meta();
        auto& sim = data.sim();
        // cppcheck-suppress constVariable
        auto& mask_node = sim.repo().get_int_field("mask_node");
        // cppcheck-suppress constVariable
        auto& levelset = sim.repo().get_field("ib_levelset");

        auto nlevels = sim.repo().num_active_levels();
        auto geom = sim.mesh().Geom();

        amrex::Real h = wdata.height;

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& problo = geom[lev].ProbLoArray();
            const auto& dx = geom[lev].CellSizeArray();

            for (amrex::MFIter mfi(levelset(lev)); mfi.isValid(); ++mfi) {
                const auto& bx = mfi.growntilebox();
                auto phi = levelset(lev).array(mfi);
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                        const amrex::Real phi_glob = phi(i, j, k);
                        amrex::Real phi_loc = z - h;
                        phi(i, j, k) = std::min(phi_loc, phi_glob);
                    });

                const auto& nbx = mfi.nodaltilebox();
                auto epsilon_node = mask_node(lev).array(mfi);
                amrex::ParallelFor(
                    nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real z = problo[2] + k * dx[2];
                        if (z <= h) {
                            epsilon_node(i, j, k) = 0;
                        }
                    });
            }
        }
    }
};

template <>
struct ComputeSourceOp<Flat>
{
    void operator()(
        Flat::DataType& data,
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Geometry& geom)
    {
        Field& source = data.sim().repo().get_field("ib_src_term");
        const Field& level = data.sim().repo().get_field("ib_levelset");
        const Field& normal = data.sim().repo().get_field("ib_normal");
        const Field& velocity = data.sim().repo().get_field("velocity");

        const auto& bx = mfi.tilebox();
        const auto& varr = velocity(lev).array(mfi);
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
                const auto gauss_fac = actuator::utils::gaussian1d(
                    amrex::Math::abs(phi_arr(i, j, k)), eps);
                const amrex::Real u = varr(i, j, k, 0);
                const amrex::Real v = varr(i, j, k, 1);
                const amrex::Real w = varr(i, j, k, 2);
                const amrex::Real vmag = std::sqrt(u * u + v * v + w * w);
                if (vmag > 0.0) {
                    // Get v dot n
                    const amrex::Real nx = norm_arr(i, j, k, 0);
                    const amrex::Real ny = norm_arr(i, j, k, 1);
                    const amrex::Real nz = norm_arr(i, j, k, 2);
                    const amrex::Real vdotn = u * nx + v * ny + w * nz;
                    // Get v_tan = v - (v dot n) n
                    const amrex::Real ut = u - vdotn * nx;
                    const amrex::Real vt = v - vdotn * ny;
                    const amrex::Real wt = w - vdotn * nz;
                    // Get tangent unit vector
                    const amrex::Real vtmag =
                        std::sqrt(ut * ut + vt * vt + wt * wt);
                    const amrex::Real tx = ut / vtmag;
                    const amrex::Real ty = vt / vtmag;
                    const amrex::Real tz = wt / vtmag;
                    // Calculate drag force
                    const amrex::Real force = -gauss_fac * vtmag * vtmag * cd;
                    sarr(i, j, k, 0) += force * tx;
                    sarr(i, j, k, 1) += force * ty;
                    sarr(i, j, k, 2) += force * tz;
                }
            });
    }
};

} // namespace ops
} // namespace ib
} // namespace amr_wind

#endif /* FLAT_OPS_H */
