#ifndef COMPLEX_TERRAIN_H
#define COMPLEX_TERRAIN_H

#include "amr-wind/immersed_boundary/IBTypes.H"
#include "amr-wind/core/vs/vector_space.H"

#include <string>

namespace amr_wind {
namespace ib {

/** Base data representation of an IB Complex Terrain geometry.
 *
 *  \ingroup ib
 */
struct ComplexTerrainBaseData
{

    bool is_mms{false};

    bool has_wall_model{false};

    //! IB velocity in the interior (solid cells)
    amrex::Vector<amrex::Real> vel_bc{{0.0, 0.0, 0.0}};

    //! Total integrated drag forces over the complex terrain
    amrex::Vector<amrex::Real> frc{{0.0, 0.0, 0.0}};
    amrex::Real band_width;
    amrex::Real sample_height;
    amrex::Real drag_coefficient;
    amrex::Real kappa{0.41};
    amrex::Real roughness_height;
};

struct ComplexTerrainType : public IBType
{};

} // namespace ib
} // namespace amr_wind

#endif /* COMPLEX_TERRAIN_H */
