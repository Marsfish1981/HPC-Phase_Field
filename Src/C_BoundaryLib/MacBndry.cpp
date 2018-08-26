
#include <winstd.H>
#include <LO_BCTYPES.H>
#include <MacBndry.H>

MacBndry::MacBndry ()
    :
    InterpBndryData()
{
    BoxLib::Abort("*** Calling default constructor for MacBndry()");
}

MacBndry::MacBndry (const BoxArray& _grids,
                    int             _ncomp,
                    const Geometry& _geom,
		    ParallelDescriptor::Color color)
    :
    InterpBndryData(_grids,_ncomp,_geom,color)
{}

MacBndry::~MacBndry () {}

void
MacBndry::setBndryConds (const BCRec& phys_bc,
                         int          ratio)
{
    const IntVect& ratio_vect = ratio * IntVect::TheUnitVector();
    setBndryConds(phys_bc, ratio_vect);
}

void
MacBndry::setBndryConds (const BCRec&   phys_bc,
                         const IntVect& ratio,
			 int            comp)
{
    m_phys_bc = phys_bc;

    //
    // ALL BCLOC VALUES ARE NOW DEFINED AS A LENGTH IN PHYSICAL
    // DIMENSIONS *RELATIVE* TO THE FACE, NOT IN ABSOLUTE PHYSICAL SPACE
    //
    const BoxArray& grids  = boxes();
    const Real*     dx     = geom.CellSize();
    const Box&      domain = geom.Domain();
    //
    // We note that all orientations of the FabSets have the same distribution.
    // We'll use the low 0 side as the model.
    //
    for (FabSetIter fsi(bndry[Orientation(0,Orientation::low)]); fsi.isValid(); ++fsi)
    {
        const int                  i     = fsi.index();
        const Box&                 grd   = grids[i];
        RealTuple&                 bloc  = bcloc[i];
        Array< Array<BoundCond> >& bctag = bcond[i];

        for (OrientationIter fi; fi; ++fi)
        {
            const Orientation face  = fi();
            const int         dir   = face.coordDir();

            if (domain[face] == grd[face] && !geom.isPeriodic(dir))
            {
                //
                // All physical bc values are located on face.
                //
                const int p_bc  = (face.isLow() ? phys_bc.lo(dir) : phys_bc.hi(dir));

                bctag[face][comp] = (p_bc == Outflow) ? LO_DIRICHLET : LO_NEUMANN;
                bloc[face]        = 0;
            }
            else
            {
                //
                // Internal bndry.
                //
                const Real delta = dx[dir]*ratio[dir];

                bctag[face][comp] = LO_DIRICHLET;
		bloc[face]        = 0.5*delta;
            }
        }
    }
}

void
MacBndry::setHomogValues (const BCRec&   bc,
                          const IntVect& ratio)
{
    setBndryConds(bc, ratio);
 
    for (OrientationIter fi; fi; ++fi)
    {
        const Orientation face  = fi();
        
        for (FabSetIter fsi(bndry[face]); fsi.isValid(); ++fsi)
        {
            bndry[face][fsi].setVal(0);
        }
    }
}

