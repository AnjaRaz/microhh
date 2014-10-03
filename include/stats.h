/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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

#ifndef STATS
#define STATS

#include <netcdfcpp.h>

// forward declarations to reduce compilation time
class Master;
class Model;
class Grid;
class Fields;

// struct for profiles
struct profvar
{
  NcVar *ncvar;
  double *data;
};

// struct for time series
struct tseriesvar
{
  NcVar *ncvar;
  double data;
};

// typedefs for containers of profiles and time series
typedef std::map<std::string, profvar> profmap;
typedef std::map<std::string, tseriesvar> tseriesmap;

// structure
struct mask
{
  std::string name;
  NcFile *dataFile;
  NcDim  *z_dim, *zh_dim, *t_dim;
  NcVar  *t_var, *iter_var;
  profmap profs;
  tseriesmap tseries;
};

typedef std::map<std::string, mask> maskmap;

class Stats
{
  public:
    Stats(Model *, Input *);
    ~Stats();

    void init(double);
    void create(int);

    unsigned long gettimelim(unsigned long);
    int getmask(Field3d *, Field3d *, mask *);
    int exec(int, double, unsigned long);
    int dostats();
    std::string getsw();

    // container for all stats, masks as uppermost in hierarchy
    maskmap masks;
    int *nmask;
    int *nmaskh;
    int nmaskbot;

    // interface functions
    // profmap profs;
    // tseriesmap tseries;

    void addmask(const std::string);
    int addprof(std::string, std::string, std::string, std::string);
    int addfixedprof(std::string, std::string, std::string, std::string, double *);
    int addtseries(std::string, std::string, std::string);

    int calcarea   (double *, const int[3], int *);

    void calcmean  (double * const, const double * const,
                    const double, const int[3],
                    const double * const, const int * const);
    void calcmean2d(double * const, const double * const,
                    const double,
                    const double * const, const int * const);

    int calcmoment  (double *, double *, double *, double, const int[3], double *, int *);

    int calDiff_2nd(double *, double *, double *, double, const int[3], double *, int *);
    int calDiff_2nd(double *, double *, double *, double *, double *, double *, double *, double, const int[3], double *, int *);
    int calDiff_4th(double *, double *, double *, double, const int[3], double *, int *);

    int calcgrad_2nd(double *, double *, double *, const int[3], double *, int *);
    int calcgrad_4th(double *, double *, double *, const int[3], double *, int *);

    int calcflux_2nd(double *, double *, double *, double *, double *, double *, const int[3], double *, int *);
    int calcflux_4th(double *, double *, double *, double *, const int[3], double *, int *);

    int addfluxes   (double *, double *, double *);
    int calccount   (double *, double *, double, double *, int *);
    int calcpath    (double *, double *, int *, double *);
    int calccover   (double *, double *, int *, double *, double);

    int calcsortprof(double *, double *, double *);

  private:
    // NcFile *dataFile;
    // NcDim  *z_dim, *zh_dim, *t_dim;
    // NcVar  *t_var, *iter_var;

    int nstats;

    // mask calculations
    int calcmask(double *, double *, double *, int *, int *, int *);

  protected:
    Model  *model;
    Grid   *grid;
    Fields *fields;
    Master *master;

    double sampletime;
    unsigned long isampletime;

    std::string swstats;
};
#endif
