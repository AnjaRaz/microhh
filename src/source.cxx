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

#include <cstdio>
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "source.h"
#include "defines.h"

namespace
{
    // Function to add external source to 3d fields. Add additional required variables...
    // void add_source(double* restrict at, const double* restrict a)
    // {
    // }
    /*
    isize_block = 4;
    jsize_block = 5
    ksize_block = 6

    block(isize_block*jsize_block*ksize_block);

    ijk = i + j*isize_block + k*isize_block*jsize_block;

    block[ijk] = a;
    */

    struct Shape
    {
        std::vector<double> range_x;
        std::vector<double> range_y;
        std::vector<double> range_z;
    };


    std::vector<double> calc_shape(const double* restrict xh, const double x0, const double sigma_x, int istart, int iend)
    {
        std::vector<double> range(2);

        int i = istart;
        for (; i<iend; ++i)
        {
            if ( std::abs(xh[i]-x0) < 4*sigma_x )
            {
                range[0] = i;
                break;
            }
        }

        for (; i<iend; ++i)
            if ( std::abs(xh[i]-x0) > 4*sigma_x )
            {
                range[1] = i;
                break;
            }

        return range;
    }
}

// Constructor: read values from ini file that do not need info from other classes
Source::Source(Master& master, Grid& grid, Fields& fields, Input& input) :
    master(master), grid(grid), fields(fields)
{
    int nerror = 0;
    nerror += input.get_item(&swsource, "source", "swsource", "", "0");

    if (swsource == "1")
    {
        nerror += input.get_list(&source_x0  , "source", "source_x0"  , "");
        nerror += input.get_list(&source_y0  , "source", "source_y0"  , "");
        nerror += input.get_list(&source_z0  , "source", "source_z0"  , "");
        nerror += input.get_list(&sigma_x    , "source", "sigma_x"    , "");
        nerror += input.get_list(&sigma_y    , "source", "sigma_y"    , "");
        nerror += input.get_list(&sigma_z    , "source", "sigma_z"    , "");
        nerror += input.get_list(&strength   , "source", "strength"   , "");
    }

    if (nerror)
        throw 1;
}

Source::~Source()
{
}

// Init function: allocate memory that you need
void Source::init()
{
    if (swsource == "1")
    {

    }
}

// Create function: read information from ini file that does need info from other class.
void Source::create(Input& input)
{
    double norm;

    std::vector<Shape> shape(source_x0.size());
    // Next functions are called in exec()?

    for (int n=0; source_x0.size(); ++n)
    {
        // Shape of the source in each direction
        shape[n].range_x = calc_shape(grid.xh, source_x0[n], sigma_x[n], grid.istart, grid.iend);
        shape[n].range_y = calc_shape(grid.yh, source_y0[n], sigma_y[n], grid.jstart, grid.jend);
        shape[n].range_z = calc_shape(grid.zh, source_z0[n], sigma_z[n], grid.kstart, grid.kend);

        // The norm
        norm = calc_norm(grid.xh, source_x0[n], sigma_x[n], grid.yh, source_y0[n], sigma_y[n], grid.zh, source_z0[n], sigma_z[n], shape[n].range_x, shape[n].range_y, shape[n].range_z);

    }
}

// Add the source to the fields. This function is called in the main time loop.
void Source::exec()
{
}


double Source::calc_norm(const double* const restrict xh, const double x0, const double sigma_x,
                         const double* const restrict yh, const double y0, const double sigma_y,
                         const double* const restrict zh, const double z0, const double sigma_z,
                         std::vector<double> range_x, std::vector<double>range_y, std::vector<double> range_z)
{
    double sum = 0;
    double blob_norm;

    for (int k = grid.kstart; k = grid.kend; ++k)
        for (int j = grid.jstart; j = grid.jend; ++j)
            for (int i = grid.istart; i = grid.iend; ++i)
            {
                // const int ijk = i + j*jj + k*kk;

                double blob_norm = 0;
                if (i>range_x[0] and i<range_x[1] and j>range_y[0] and j<range_y[1] and k>range_z[0] and k<range_z[1])
                    blob_norm = exp(-pow(xh[i]-x0,2)/pow(sigma_x,2) - pow(yh[i]-y0,2)/pow(sigma_y,2) - pow(zh[k]-z0,2)/pow(sigma_z,2));

                sum += blob_norm*grid.dx*grid.dy*grid.dz[k];
        }

    return sum;
}


std::vector<double> Source::calc_source(const double* const restrict xh, const double x0, const double sigma_x,
                            const double* const restrict yh, const double y0, const double sigma_y,
                            const double* const restrict zh, const double z0, const double sigma_z,
                            double* range_x, double* range_y, double* range_z,
                            const int strength, const double norm)
{

    double sum = 0;
    const int ii = 1;
    const int jj = grid.icells;
    const int kk = grid.ijcells;
    std::vector<double> blob;


    for (int k = grid.kstart; k = grid.kend; ++k)
        for (int j = grid.jstart; j = grid.jend; ++j)
            for (int i = grid.istart; i = grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                if (i>range_x[0] and i<range_x[1] and j>range_y[0] and j<range_y[1] and k>range_z[0] and k<range_z[1])
                    blob[ijk] = strength/norm*exp(-pow(xh[i]-x0,2)/pow(sigma_x,2) - pow(yh[i]-y0,2)/pow(sigma_y,2) - pow(zh[k]-z0,2)/pow(sigma_z,2));

                else
                    blob[ijk] = 0;

                sum+=sum+blob[ijk];
            }

    return blob;
}
