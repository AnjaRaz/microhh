/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef THERMO_MOIST_FUNCTIONS
#define THERMO_MOIST_FUNCTIONS

// In case the code is compiled with NVCC, add the macros for CUDA
#ifdef __CUDACC__
#  define CUDA_MACRO __host__ __device__
#else
#  define CUDA_MACRO
#endif

#include "constants.h"

namespace Thermo_moist_functions
{
    using namespace Constants;

    // INLINE FUNCTIONS
    CUDA_MACRO inline double buoyancy(const double exn, const double thl, const double qt, const double ql, const double thvref)
    {
        return grav * ((thl + Lv*ql/(cp*exn)) * (1. - (1. - Rv/Rd)*qt - Rv/Rd*ql) - thvref) / thvref;
    }

    CUDA_MACRO inline double buoyancy_no_ql(const double thl, const double qt, const double thvref)
    {
        return grav * (thl * (1. - (1. - Rv/Rd)*qt) - thvref) / thvref;
    }

    CUDA_MACRO inline double buoyancy_flux_no_ql(const double thl, const double thlflux, const double qt, const double qtflux, const double thvref)
    {
        return grav/thvref * (thlflux * (1. - (1.-Rv/Rd)*qt) - (1.-Rv/Rd)*thl*qtflux);
    }

    CUDA_MACRO inline double esat(const double T)
    {
        #ifdef __CUDACC__
        const double x=fmax(-80.,T-T0);
        #else
        const double x=std::max(-80.,T-T0);
        #endif

        return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))));
    }

    CUDA_MACRO inline double qsat(const double p, const double T)
    {
        return ep*esat(T)/(p-(1-ep)*esat(T));
    }

    CUDA_MACRO inline double exner(const double p)
    {
        return pow((p/p0),(Rd/cp));
    }
}
#endif
