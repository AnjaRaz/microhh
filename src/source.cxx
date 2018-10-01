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
#include <fstream>
#include <cmath>
#include <vector>
#include <array>
#include <list>
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "source.h"
#include "defines.h"

namespace
{

    struct Shape
    {
        std::vector<double> range_x;
        std::vector<double> range_y;
        std::vector<double> range_z;
    };

    std::vector<Shape> shape;
    std::vector<double> blob;
    std::vector<double> norm;

    std::vector<double> calc_shape(const double* restrict x, const double x0, const double sigma_x, int istart, int iend)
    {
        std::vector<double> range(2);

        int i = istart;
        for (; i<iend; ++i)
        {
            if ( std::abs(x[i]-x0) < 4*sigma_x )
            {
                range[0] = i;
                break;
            }
        }

        for (; i<iend; ++i)
            if ( std::abs(x[i]-x0) > 4*sigma_x )
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
        nerror += input.get_list(&sourcelist  , "source", "sourcelist"   , "");
        nerror += input.get_list(&source_x0   , "source", "source_x0"    , "");
        nerror += input.get_list(&source_y0   , "source", "source_y0"    , "");
        nerror += input.get_list(&source_z0   , "source", "source_z0"    , "");
        nerror += input.get_list(&sigma_x     , "source", "sigma_x"      , "");
        nerror += input.get_list(&sigma_y     , "source", "sigma_y"      , "");
        nerror += input.get_list(&sigma_z     , "source", "sigma_z"      , "");
        nerror += input.get_list(&strength    , "source", "strength"     , "");
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
    shape.resize(source_x0.size());
    norm.resize(source_x0.size());
    
    for (int n=0; n<source_x0.size(); ++n)
    {
        // Shape of the source in each direction
        shape[n].range_x = calc_shape(grid.x, source_x0[n], sigma_x[n], grid.istart, grid.iend);
        shape[n].range_y = calc_shape(grid.y, source_y0[n], sigma_y[n], grid.jstart, grid.jend);
        shape[n].range_z = calc_shape(grid.z, source_z0[n], sigma_z[n], grid.kstart, grid.kend);

        norm[n] = calc_norm(grid.x, source_x0[n], sigma_x[n], grid.y, source_y0[n], sigma_y[n], grid.z, source_z0[n], sigma_z[n], shape[n].range_x, shape[n].range_y, shape[n].range_z);
        
    }
    blob.resize(grid.ncells);
}

// Add the source to the fields. This function is called in the main time loop.
void Source::exec()
{
    
    for (int n=0; n<sourcelist.size(); ++n)
    {
        calc_source(blob.data(), grid.x, source_x0[n], sigma_x[n], grid.y, source_y0[n], sigma_y[n], grid.z, source_z0[n], sigma_z[n], shape[n].range_x, shape[n].range_y, shape[n].range_z,
                   strength[n], norm[n]);
            

        add_source(fields.st[sourcelist[n]]->data, blob.data());
    }
}


double Source::calc_norm(const double* const restrict x, const double x0, const double sigma_x,
                         const double* const restrict y, const double y0, const double sigma_y,
                         const double* const restrict z, const double z0, const double sigma_z,
                         std::vector<double> range_x, std::vector<double>range_y, std::vector<double> range_z)
{
    double sum = 0.;
    double blob_norm = 0.;

    for (int k = grid.kstart; k < grid.kend; ++k)
        for (int j = grid.jstart; j < grid.jend; ++j)
            for (int i = grid.istart; i < grid.iend; ++i)
            {
                if (i>=range_x[0] && i<=range_x[1] && j>=range_y[0] && j<=range_y[1] && k>=range_z[0] && k<=range_z[1])
                {
                    blob_norm = exp(-pow(x[i]-x0,2.0)/pow(sigma_x,2.0) - pow(y[j]-y0,2.0)/pow(sigma_y,2.0) - pow(z[k]-z0,2.0)/pow(sigma_z,2.0));
                }
                sum += blob_norm*grid.dx*grid.dy*grid.dz[k];
            }

    return sum;
}


void Source::calc_source(double* const restrict blob,
                         const double* const restrict x, const double x0, const double sigma_x,
                         const double* const restrict y, const double y0, const double sigma_y,
                         const double* const restrict z, const double z0, const double sigma_z,
                         std::vector<double> range_x, std::vector<double>range_y, std::vector<double> range_z,
                         const double strength, double norm)
{

    double sum = 0.;
    const int ii = 1;
    const int jj = grid.icells;
    const int kk = grid.ijcells;



    for (int k = grid.kstart; k < grid.kend; ++k)
        for (int j = grid.jstart; j < grid.jend; ++j)
            for (int i = grid.istart; i < grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                if (i>=range_x[0] && i<=range_x[1] && j>=range_y[0] && j<=range_y[1] && k>=range_z[0] && k<=range_z[1])
                {
                    blob[ijk] = strength/norm*exp(-pow(x[i]-x0,2.0)/pow(sigma_x,2.0) - pow(y[j]-y0,2.0)/pow(sigma_y,2.0) - pow(z[k]-z0,2.0)/pow(sigma_z,2.0));
                }
                else
                    blob[ijk] = 0.;


                sum += blob[ijk]*grid.dx*grid.dy*grid.dz[k];
                
            }

}

void Source:: add_source(double* const restrict st, const double* const restrict)
{
    const int jj = grid.icells;
    const int kk = grid.ijcells;

    for (int k=grid.kstart; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                st[ijk] += blob[ijk];


            }

}
