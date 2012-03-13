#ifndef FIELDS
#define FIELDS

#include "grid.h"
#include "field3d.h"

class cfields
{
  public:
    // functions
    cfields(cgrid *);
    ~cfields();
    int initfields();
    int createfields();

    int resettend();
    int boundary();

    int check();
    inline double interp2(const double, const double);

    // variables
    double *flow;
    double *flowt;

    cfield3d *u;
    cfield3d *v;
    cfield3d *w;

    cfield3d *ut;
    cfield3d *vt;
    cfield3d *wt;

    double visc;
  private:
    // variables
    cgrid *grid;

    // functions
    double calcmom(double *, double *, double *, double *);
    double calctke(double *, double *, double *, double *);
};
#endif
