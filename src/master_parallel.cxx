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

#ifdef USEMPI
#include <mpi.h>
#include <stdexcept>
#include "grid.h"
#include "defines.h"
#include "master.h"

Master::Master()
{
    initialized = false;
    allocated   = false;

    // set the mpiid, to ensure that errors can be written if MPI init fails
    mpiid = 0;
}

Master::~Master()
{
    if (allocated)
    {
        delete[] reqs;
        MPI_Comm_free(&commxy);
        MPI_Comm_free(&commx);
        MPI_Comm_free(&commy);
    }

    print_message("Finished run on %d processes\n", nprocs);

    if (initialized)
        MPI_Finalize();
}

void Master::start(int argc, char *argv[])
{
    // initialize the MPI
    int n = MPI_Init(NULL, NULL);
    if (check_error(n))
        throw 1;

    wall_clock_start = get_wall_clock_time();

    initialized = true;

    // get the rank of the current process
    n = MPI_Comm_rank(MPI_COMM_WORLD, &mpiid);
    if (check_error(n))
        throw 1;

    // get the total number of processors
    n = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (check_error(n))
        throw 1;

    // store a temporary copy of COMM_WORLD in commxy
    n = MPI_Comm_dup(MPI_COMM_WORLD, &commxy);
    if (check_error(n))
        throw 1;

    print_message("Starting run on %d processes\n", nprocs);

    // process the command line options
    if (argc <= 1)
    {
        print_error("Specify init, run or post mode\n");
        throw 1;
    }
    else
    {
        // check the execution mode
        mode = argv[1];
        if (mode != "init" && mode != "run" && mode != "post")
        {
            print_error("Specify init, run or post mode\n");
            throw 1;
        }
        // set the name of the simulation
        if (argc > 2)
            simname = argv[2];
        else
            simname = "microhh";
    }
}

void Master::init(Input *inputin)
{
    int nerror = 0;
    nerror += inputin->get_item(&npx, "master", "npx", "", 1);
    nerror += inputin->get_item(&npy, "master", "npy", "", 1);

    // Get the wall clock limit with a default value of 1E8 hours, which will be never hit
    double wall_clock_limit;
    nerror += inputin->get_item(&wall_clock_limit, "master", "wallclocklimit", "", 1E8);

    if (nerror)
        throw 1;

    wall_clock_end = wall_clock_start + 3600.*wall_clock_limit;

    if (nprocs != npx*npy)
    {
        print_error("nprocs = %d does not equal npx*npy = %d*%d\n", nprocs, npx, npy);
        throw 1;
    }

    int n;
    int dims    [2] = {npy, npx};
    int periodic[2] = {true, true};

    // define the dimensions of the 2-D grid layout
    n = MPI_Dims_create(nprocs, 2, dims);
    if (check_error(n))
        throw 1;

    // create a 2-D grid communicator that is optimized for grid to grid transfer
    // first, free our temporary copy of COMM_WORLD
    n = MPI_Comm_free(&commxy);
    if (check_error(n))
        throw 1;

    // for now, do not reorder processes, blizzard gives large performance loss
    n = MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, false, &commxy);
    if (check_error(n))
        throw 1;

    n = MPI_Comm_rank(commxy, &mpiid);
    if (check_error(n))
        throw 1;

    // retrieve the x- and y-coordinates in the 2-D grid for each process
    int mpicoords[2];
    n = MPI_Cart_coords(commxy, mpiid, 2, mpicoords);
    if (check_error(n))
        throw 1;

    mpicoordx = mpicoords[1];
    mpicoordy = mpicoords[0];

    int dimx[2] = {false, true };
    int dimy[2] = {true , false};

    n = MPI_Cart_sub(commxy, dimx, &commx);
    if (check_error(n))
        throw 1;

    n = MPI_Cart_sub(commxy, dimy, &commy);
    if (check_error(n))
        throw 1;

    // find out who are the neighbors of this process to facilitate the communication routines
    n = MPI_Cart_shift(commxy, 1, 1, &nwest , &neast );
    if (check_error(n))
        throw 1;

    n = MPI_Cart_shift(commxy, 0, 1, &nsouth, &nnorth);
    if (check_error(n))
        throw 1;

    // create the requests arrays for the nonblocking sends
    int npmax;
    npmax = std::max(npx, npy);

    // have at least as many communicators as prognostic variables
    npmax = std::max(npmax, 8*4);
    reqs  = new MPI_Request[npmax*2];
    reqsn = 0;

    allocated = true;
}

double Master::get_wall_clock_time()
{
    return MPI_Wtime();
}

int Master::check_error(int n)
{
    char errbuffer[MPI_MAX_ERROR_STRING];
    int errlen;

    if (n != MPI_SUCCESS)
    {
        MPI_Error_string(n, errbuffer, &errlen);
        print_error("MPI: %s\n", errbuffer);
        return 1;
    }

    return 0;
}

void Master::wait_all()
{
    // wait for MPI processes and reset the number of pending requests
    MPI_Waitall(reqsn, reqs, MPI_STATUSES_IGNORE);
    reqsn = 0;
}

// do all broadcasts over the MPI_COMM_WORLD, to avoid complications in the input file reading
void Master::broadcast(char *data, int datasize)
{
    MPI_Bcast(data, datasize, MPI_CHAR, 0, commxy);
}

// overloaded broadcast functions
void Master::broadcast(int *data, int datasize)
{
    MPI_Bcast(data, datasize, MPI_INT, 0, commxy);
}

void Master::broadcast(unsigned long *data, int datasize)
{
    MPI_Bcast(data, datasize, MPI_UNSIGNED_LONG, 0, commxy);
}

void Master::broadcast(double *data, int datasize)
{
    MPI_Bcast(data, datasize, MPI_DOUBLE, 0, commxy);
}

void Master::sum(int *var, int datasize)
{
    MPI_Allreduce(MPI_IN_PLACE, var, datasize, MPI_INT, MPI_SUM, commxy);
}

void Master::sum(double *var, int datasize)
{
    MPI_Allreduce(MPI_IN_PLACE, var, datasize, MPI_DOUBLE, MPI_SUM, commxy);
}

void Master::max(double *var, int datasize)
{
    MPI_Allreduce(MPI_IN_PLACE, var, datasize, MPI_DOUBLE, MPI_MAX, commxy);
}

void Master::min(double *var, int datasize)
{
    MPI_Allreduce(MPI_IN_PLACE, var, datasize, MPI_DOUBLE, MPI_MIN, commxy);
}
#endif
