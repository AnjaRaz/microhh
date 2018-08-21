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

#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "source.h"
#include "defines.h"

namespace
{
    // Function to add external source to 3d fields. Add additional required variables...
    void add_source(double* restrict at, const double* restrict a)
    {
    }
}

// Constructor: read values from ini file that do not need info from other classes
Source::Source(Master& master, Grid& grid, Fields& fields, Input& input) :
    master(master), grid(grid), fields(fields)
{
    int nerror = 0;
    // nerror += inputin.get_item(&swbuffer, "buffer", "swbuffer", "", "0");

    if (nerror)
        throw 1;
}

Source::~Source()
{
}

// Init function: allocate memory that you need
void Source::init()
{
}

// Create function: read information from ini file that does need info from other class.
void Source::create(Input& input)
{
}

// Add the source to the fields. This function is called in the main time loop.
void Source::exec()
{
}
