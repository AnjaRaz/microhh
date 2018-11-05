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
        std::vector<int> range_x;
        std::vector<int> range_y;
        std::vector<int> range_z;
    };

    std::vector<Shape> shape;
    std::vector<double> blob;
    std::vector<double> norm;

    bool is_blob(const double x, const double x0, const double line_x, const double sigma_x)
    {
        bool in = false;

        if ( std::abs(x-x0) < 4*sigma_x && x-x0-line_x < 4*sigma_x )
            in = true;

        return in;
    }

    std::vector<int> calc_shape(const double* restrict x, const double x0, const double sigma_x, const double line_x, int istart, int iend, const int mpiid)
    {
        std::vector<int> range(2);
        /*
        range[0] = istart;
        range[1] = istart;

        int i = istart;

        for (; i<iend; ++i)
        {
            if (is_blob(x[i], x0, line_x, sigma_x) == 0 && x[i] <= x0)
            {
                ++range[0];
                ++range[1];
            }
            else if ()

            while (is_blob(x[i], x0, line_x, sigma_x) == 0 )
                range[0] = range[0] + 1;
            while (is_blob(x[i], x0, line_x, sigma_x) == 1 )
                range[1] = range[0] + 1;

        }*/


        //std::cout<<"istart i ("<<mpiid<<"):"<<i<<std::endl;
        int i = istart;
        range[0] = iend;
        for (; i<iend; ++i)
        {
            std::cout<<"bla(" << mpiid << "): " <<"First loop, i"<<" "<<i<<" "<< "x[i]:"<<x[i]<<std::endl;
            //std::cout<<"i in the first loop:"<<i <<std::endl;
            if ( x[i]-x0 + 4*sigma_x > 0 )
            {
                range[0] = i;
                //std::cout<<"range[0] for ("<<mpiid<<"):"<< i << std::endl;
                break;
            }

        }

        //std::cout<<"After first loop istart, iend and i for bla (" << mpiid << "): "<<istart<<" " <<iend<<" "<<i<<std::endl;
        i = istart;
        //std::cout<<"New i:"<<i<<std::endl;
        for (; i<iend; ++i)
        {
            range[1] = iend;
            //std::cout<<"iend in the second loop, bla (" << mpiid << "): " <<istart<<" "<<iend<<std::endl;
            //std::cout<<"i in the second loop, before if, bla (" << mpiid << "): "<< i <<std::endl;
            std::cout<<"bla(" << mpiid << "): " <<"Second loop, i"<<" "<<i<<" "<< "x[i]:"<<x[i]<<std::endl;

            if (x[i]-x0-line_x - 4*sigma_x > 0)
            {
                range[1] = i;
                break;
            }

            //std::cout<<"i in the second loop, after if, bla (" << mpiid << "): "<< i <<std::endl;

         }
        /*if (range[0] == 0 && range[1] == 0 )
            istart = iend;*/

        std::cout<<"bla (" << mpiid << "): " << range[0]<<" "<<range[1]<<" "<<istart<<" "<<iend<<std::endl;
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
        nerror += input.get_list(&line_x      , "source", "line_x"       , "");
        nerror += input.get_list(&line_y      , "source", "line_y"       , "");
        nerror += input.get_list(&line_z      , "source", "line_z"       , "");
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
        shape[n].range_x = calc_shape(grid.x, source_x0[n], sigma_x[n], line_x[n], grid.istart, grid.iend, master.mpiid);
        shape[n].range_y = calc_shape(grid.y, source_y0[n], sigma_y[n], line_y[n], grid.jstart, grid.jend, master.mpiid);
        shape[n].range_z = calc_shape(grid.z, source_z0[n], sigma_z[n], line_z[n], grid.kstart, grid.kend, master.mpiid);

        norm[n] = calc_norm(grid.x, source_x0[n], sigma_x[n], line_x[n], 
                            grid.y, source_y0[n], sigma_y[n], line_y[n], 
                            grid.z, source_z0[n], sigma_z[n], line_z[n],
                            shape[n].range_x, shape[n].range_y, shape[n].range_z);
        
    }
    blob.resize(grid.ncells);
}

// Add the source to the fields. This function is called in the main time loop.
void Source::exec()
{
    
    for (int n=0; n<sourcelist.size(); ++n)
    {
        calc_source(blob.data(), grid.x, source_x0[n], sigma_x[n], line_x[n], 
                    grid.y, source_y0[n], sigma_y[n], line_y[n], 
                    grid.z, source_z0[n], sigma_z[n], line_z[n], 
                    shape[n].range_x, shape[n].range_y, shape[n].range_z,
                    strength[n], norm[n]);
            

        add_source(fields.st[sourcelist[n]]->data, blob.data(), shape[n].range_x, shape[n].range_y, shape[n].range_z);
    }
}


double Source::calc_norm(const double* const restrict x, const double x0, const double sigma_x, const double line_x,
                         const double* const restrict y, const double y0, const double sigma_y, const double line_y,
                         const double* const restrict z, const double z0, const double sigma_z, const double line_z,
                         std::vector<int> range_x, std::vector<int>range_y, std::vector<int> range_z)
{
    double sum = 0.;
    double blob_norm = 0.;
    const int jj = grid.icells;
    const int kk = grid.ijcells;

    for (int k = grid.kstart; k < grid.kend; ++k)
        for (int j = grid.jstart; j < grid.jend; ++j)
            for (int i = grid.istart; i < grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                if (i>=range_x[0] && i<=range_x[1] && j>=range_y[0] && j<=range_y[1] && k>=range_z[0] && k<=range_z[1])
                {
                    if (line_x != 0)
                    {
                        if (x[i] >= x0+line_x)
                            blob_norm = exp(-pow(x[i]-x0-line_x,2.0)/pow(sigma_x,2.0) - pow(y[j]-y0-line_y,2.0)/pow(sigma_y,2.0) - pow(z[k]-z0-line_z,2.0)/pow(sigma_z,2.0));
                        else if (x[i]<=x0)
                            blob_norm = exp(-pow(x[i]-x0,2.0)/pow(sigma_x,2.0) - pow(y[j]-y0-line_y,2.0)/pow(sigma_y,2.0) - pow(z[k]-z0-line_z,2.0)/pow(sigma_z,2.0));
                        else
                            blob_norm = exp(- pow(y[j]-y0-line_y,2.0)/pow(sigma_y,2.0) - pow(z[k]-z0-line_z,2.0)/pow(sigma_z,2.0));
                    }
                    else if (line_y != 0)
                    {
                        if (y[j] >= y0+line_y)
                            blob_norm = exp(-pow(x[i]-x0-line_x,2.0)/pow(sigma_x,2.0) - pow(y[j]-y0-line_y,2.0)/pow(sigma_y,2.0) - pow(z[k]-z0-line_z,2.0)/pow(sigma_z,2.0));
                        else if (y[j]<=y0)
                            blob_norm = exp(-pow(x[i]-x0-line_x,2.0)/pow(sigma_x,2.0) - pow(y[j]-y0,2.0)/pow(sigma_y,2.0) - pow(z[k]-z0-line_z,2.0)/pow(sigma_z,2.0));
                        else
                            blob_norm = exp(-pow(x[i]-x0-line_x,2.0)/pow(sigma_x,2.0) - pow(z[k]-z0-line_z,2.0)/pow(sigma_z,2.0));
                        
                    }
                    else if (line_z != 0)
                    {
                        if (z[k] >= z0+line_z)
                            blob_norm = exp(-pow(x[i]-x0-line_x,2.0)/pow(sigma_x,2.0) - pow(y[j]-y0-line_y,2.0)/pow(sigma_y,2.0) - pow(z[k]-z0-line_z,2.0)/pow(sigma_z,2.0));
                        else if (z[k]<=z0)
                            blob_norm = exp(-pow(x[i]-x0-line_x,2.0)/pow(sigma_x,2.0) - pow(y[j]-y0-line_y,2.0)/pow(sigma_y,2.0) - pow(z[k]-z0,2.0)/pow(sigma_z,2.0));
                        else
                            blob_norm = exp(-pow(x[i]-x0-line_x,2.0)/pow(sigma_x,2.0) - pow(y[j]-y0-line_y,2.0)/pow(sigma_y,2.0));
                    }
                    else
                        blob_norm = exp(-pow(x[i]-x0-line_x,2.0)/pow(sigma_x,2.0) - pow(y[j]-y0-line_y,2.0)/pow(sigma_y,2.0) - pow(z[k]-z0-line_z,2.0)/pow(sigma_z,2.0));
                }
                sum += blob_norm*grid.dx*grid.dy*grid.dz[k];
            }

    master.sum(&sum, 1);

    return sum;
}


void Source::calc_source(double* const restrict blob,
                         const double* const restrict x, const double x0, const double sigma_x, const double line_x,
                         const double* const restrict y, const double y0, const double sigma_y, const double line_y,
                         const double* const restrict z, const double z0, const double sigma_z, const double line_z,
                         std::vector<int> range_x, std::vector<int> range_y, std::vector<int> range_z,
                         const double strength, double norm)
{

    double sum = 0.;
    const int ii = 1;
    const int jj = grid.icells;
    const int kk = grid.ijcells;

    std::ofstream file;
    file.open("tuga.txt");

    const int range[3] = {range_x[1]-range_x[0], range_y[1]-range_y[0], range_z[1]-range_z[0]};

/*
    for (int k = grid.kstart; k < grid.kend; ++k)
        for (int j = grid.jstart; j < grid.jend; ++j)
            for (int i = grid.istart; i < grid.iend; ++i)*/
 /*   for (int k = range_z[0]; k < range_z[1]; ++k)
        for (int j = range_y[0]; j < range_z[1]; ++j)
            for (int i = range_x[0]; i < range_x[1]; ++i)*/
    for(int k = 0; k < range[2]; ++k)
        for(int j = 0; j < range[1]; ++j)
            for(int i = 0; i < range[0]; ++i)
            {
                //const int ijk = i + j*jj + k*kk;  // i + j*sizei + k*sizei*sizej
                const int ijk = i + j*range[0] + k*range[0]*range[1];

                if (range_x[0] == range_x[1] || range_y[0] == range_y[1] || range_z[0] == range_z[1])
                    break;

                //std::cout<<i<<"  "<<j<<"  "<<k<<std::endl;

               // if (i>=range_x[0] && i<=range_x[1] && j>=range_y[0] && j<=range_y[1] && k>=range_z[0] && k<=range_z[1])
               // {
                 if (line_x != 0)
                 {
                     if (x[i + range_x[0]] >= x0+line_x)
                        blob[ijk] = strength/norm*exp(-pow(x[i + range_x[0]]-x0-line_x,2.0)/pow(sigma_x,2.0) - pow(y[j + range_y[0]]-y0-line_y,2.0)/pow(sigma_y,2.0) - pow(z[k + range_z[0]]-z0-line_z,2.0)/pow(sigma_z,2.0));
                     else if (x[i + range_x[0]]<=x0)
                        blob[ijk] = strength/norm*exp(-pow(x[i + range_x[0]]-x0,2.0)/pow(sigma_x,2.0) - pow(y[j + range_y[0]]-y0-line_y,2.0)/pow(sigma_y,2.0) - pow(z[k + range_z[0]]-z0-line_z,2.0)/pow(sigma_z,2.0));
                     else
                        blob[ijk] = strength/norm*exp(-pow(y[j + range_y[0]]-y0-line_y,2.0)/pow(sigma_y,2.0) - pow(z[k + range_z[0]]-z0-line_z,2.0)/pow(sigma_z,2.0));
                 }
                 else if (line_y != 0)
                 {
                    if (y[j + range_y[0]] >= y0+line_y)
                        blob[ijk] = strength/norm*exp(-pow(x[i + range_x[0]]-x0-line_x,2.0)/pow(sigma_x,2.0) - pow(y[j + range_y[0]]-y0-line_y,2.0)/pow(sigma_y,2.0) - pow(z[k + range_z[0]]-z0-line_z,2.0)/pow(sigma_z,2.0));
                    else if (y[j + range_y[0]]<=y0)
                        blob[ijk] = strength/norm*exp(-pow(x[i + range_x[0]]-x0-line_x,2.0)/pow(sigma_x,2.0) - pow(y[j + range_y[0]]-y0,2.0)/pow(sigma_y,2.0) - pow(z[k + range_z[0]]-z0-line_z,2.0)/pow(sigma_z,2.0));
                    else
                        blob[ijk] = strength/norm*exp(-pow(x[i + range_x[0]]-x0-line_x,2.0)/pow(sigma_x,2.0) - pow(z[k + range_z[0]]-z0-line_z,2.0)/pow(sigma_z,2.0));
                        
                 }
                 else if (line_z != 0)
                 {

                    if (z[k + range_z[0]] >= z0+line_z)
                        blob[ijk] = strength/norm*exp(-pow(x[i + range_x[0]]-x0-line_x,2.0)/pow(sigma_x,2.0) - pow(y[j + range_y[0]]-y0-line_y,2.0)/pow(sigma_y,2.0) - pow(z[k + range_z[0]]-z0-line_z,2.0)/pow(sigma_z,2.0));
                    else if (z[k + range_z[0]]<=z0)
                        blob[ijk] = strength/norm*exp(-pow(x[i + range_x[0]]-x0-line_x,2.0)/pow(sigma_x,2.0) - pow(y[j + range_y[0]]-y0-line_y,2.0)/pow(sigma_y,2.0) - pow(z[k + range_z[0]]-z0,2.0)/pow(sigma_z,2.0));
                    else
                        blob[ijk] = strength/norm*exp(-pow(x[i + range_x[0]]-x0-line_x,2.0)/pow(sigma_x,2.0) - pow(y[j + range_y[0]]-y0-line_y,2.0)/pow(sigma_y,2.0));
                 }
                 else
                    blob[ijk] = strength/norm*exp(-pow(x[i + range_x[0]]-x0-line_x,2.0)/pow(sigma_x,2.0) - pow(y[j + range_y[0]]-y0-line_y,2.0)/pow(sigma_y,2.0) - pow(z[k + range_z[0]]-z0-line_z,2.0)/pow(sigma_z,2.0));
                //}
                //else
                //    blob[ijk] = 0.;

                sum += blob[ijk]*grid.dx*grid.dy*grid.dz[k + range_z[0]];
                file<<blob[ijk]<<std::endl;

                

            }
    master.sum(&sum,1);
    std::cout <<sum<<std::endl;


    file.close();

}

void Source:: add_source(double* const restrict st, const double* const restrict,
                         std::vector<int> range_x, std::vector<int>range_y, std::vector<int> range_z)
{
    const int jj = grid.icells;
    const int kk = grid.ijcells;
    const int range[3] =  {range_x[1]-range_x[0], range_y[1]-range_y[0], range_z[1]-range_z[0]};


    for (int k = range_z[0]; k < range_z[1]; ++k)
        for (int j = range_y[0]; j < range_y[1]; ++j)
            for (int i = range_x[0]; i < range_x[1]; ++i)
            {
                if (range_x[0] == range_x[1] || range_y[0] == range_y[1] || range_z[0] == range_z[1])
                    break;

                const int ijk = i + j*jj + k*kk;
                const int cdf = i - range_x[0] + (j-range_y[0])*range[0] + (k - range_z[0])*range[0]*range[1];

                //std::cout<<ijk<<"  "<<cdf<<std::endl;

                st[ijk] += blob[cdf];

            }

}
