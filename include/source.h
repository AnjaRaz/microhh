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

#ifndef SOURCE_H
#define SOURCE_H

class Master;
class Grid;
class Fields;
class Input;

class Source
{
    public:
        Source(Master&, Grid&, Fields&, Input&); // Constructor of the buffer class.
        ~Source(); // Destructor of the buffer class.

        void init();
        void create(Input&);

        void exec(); // Add the tendencies created by the damping.


    private:
        Master& master; // Reference to master class.
        Grid&   grid;   // Reference to grid class.
        Fields& fields; // Reference to fields class.

        std::string swsource;


        double x0;
        double y0;
        double z0;
        std::vector<std::string> sourcelist;
        std::vector<double>      source_x0;
        std::vector<double>      source_y0;
        std::vector<double>      source_z0;
        std::vector<double>      sigma_x;
        std::vector<double>      sigma_y;
        std::vector<double>      sigma_z;
        std::vector<double>      strength;
        std::vector<double>      line_x;
        std::vector<double>      line_y;
        std::vector<double>      line_z;
        std::vector<double>      blob;

        double calc_norm(const double* const, const double, const double, const double,
                        const double* const, const double, const double, const double,
                        const double* const, const double, const double, const double,
                        std::vector<int>, std::vector<int>, std::vector<int>);

        void calc_source(double* const, const double* const, const double, const double, const double,
                        const double* const, const double, const double, const double,
                        const double* const, const double, const double, const double,
                        std::vector<int>, std::vector<int>, std::vector<int>,
                        const double, double);
                        
        void add_source(double* const, const double* const, std::vector<int>, std::vector<int>, std::vector<int>);
};
#endif
