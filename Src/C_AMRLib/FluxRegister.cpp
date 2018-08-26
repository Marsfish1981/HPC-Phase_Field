
#include <winstd.H>
#include <BArena.H>
#include <FluxRegister.H>
#include <Geometry.H>
#include <FLUXREG_F.H>
#include <ParallelDescriptor.H>
#include <BLProfiler.H>
#include <ccse-mpi.H>

#include <vector>

FluxRegister::FluxRegister ()
{
    fine_level = ncomp = -1;
    ratio = IntVect::TheUnitVector();
    ratio.scale(-1);
}

FluxRegister::FluxRegister (const BoxArray& fine_boxes, 
                            const IntVect&  ref_ratio,
                            int             fine_lev,
                            int             nvar)
{
    define(fine_boxes,ref_ratio,fine_lev,nvar);
}

FluxRegister::FluxRegister (const BoxArray&            fine_boxes, 
                            const IntVect&             ref_ratio,
                            int                        fine_lev,
                            int                        nvar,
                            const DistributionMapping& dm)
{
    define(fine_boxes,ref_ratio,fine_lev,nvar,dm);
}

const IntVect&
FluxRegister::refRatio () const
{
    return ratio;
}

int
FluxRegister::fineLevel () const
{
    return fine_level;
}

int
FluxRegister::crseLevel () const
{
    return fine_level-1;
}

int
FluxRegister::nComp () const
{
    return ncomp;
}

const BoxArray&
FluxRegister::coarsenedBoxes () const
{
    return grids;
}

void
FluxRegister::define (const BoxArray& fine_boxes, 
                      const IntVect&  ref_ratio,
                      int             fine_lev,
                      int             nvar)
{
    BL_ASSERT(fine_boxes.isDisjoint());
    BL_ASSERT(grids.size() == 0);

    ratio      = ref_ratio;
    fine_level = fine_lev;
    ncomp      = nvar;

    grids.define(fine_boxes);
    grids.coarsen(ratio);

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        const Orientation lo_face(dir,Orientation::low);
        const Orientation hi_face(dir,Orientation::high);

        IndexType typ(IndexType::TheCellType());

        typ.setType(dir,IndexType::NODE);

        BndryRegister::define(lo_face,typ,0,1,0,nvar);
        BndryRegister::define(hi_face,typ,0,1,0,nvar);
    }
}

void
FluxRegister::define (const BoxArray&            fine_boxes, 
                      const IntVect&             ref_ratio,
                      int                        fine_lev,
                      int                        nvar,
                      const DistributionMapping& dm)
{
    BL_ASSERT(fine_boxes.isDisjoint());
    BL_ASSERT(grids.size() == 0);

    ratio      = ref_ratio;
    fine_level = fine_lev;
    ncomp      = nvar;

    grids.define(fine_boxes);
    grids.coarsen(ratio);

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        const Orientation lo_face(dir,Orientation::low);
        const Orientation hi_face(dir,Orientation::high);

        IndexType typ(IndexType::TheCellType());

        typ.setType(dir,IndexType::NODE);

        BndryRegister::define(lo_face,typ,0,1,0,nvar,dm);
        BndryRegister::define(hi_face,typ,0,1,0,nvar,dm);
    }
}

FluxRegister::~FluxRegister () {}

Real
FluxRegister::SumReg (int comp) const
{
    Real sum = 0.0;

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        const FabSet& lofabs = bndry[Orientation(dir,Orientation::low) ];
        const FabSet& hifabs = bndry[Orientation(dir,Orientation::high)];

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif
        for (FabSetIter fsi(lofabs); fsi.isValid(); ++fsi)
        {
            sum += (lofabs[fsi].sum(comp) - hifabs[fsi].sum(comp));
        }
    }

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

void
FluxRegister::CrseInit (const MultiFab& mflx,
                        const MultiFab& area,
                        int             dir,
                        int             srccomp,
                        int             destcomp,
                        int             numcomp,
                        Real            mult,
                        FrOp            op)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= mflx.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    const Orientation face_lo(dir,Orientation::low);
    const Orientation face_hi(dir,Orientation::high);
 
    MultiFab mf(mflx.boxArray(),numcomp,0);

#ifdef _OPENMP
#pragma omp parallel
#endif    
    for (MFIter mfi(mflx,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();
	
        mf[mfi].copy(mflx[mfi],bx,srccomp,bx,0,numcomp);

        mf[mfi].mult(mult,bx,0,numcomp);

        for (int i = 0; i < numcomp; i++)
            mf[mfi].mult(area[mfi],bx,bx,0,i,1);
    }

    for (int pass = 0; pass < 2; pass++)
    {
        const Orientation face = ((pass == 0) ? face_lo : face_hi);

        if (op == FluxRegister::COPY)
        {
            bndry[face].copyFrom(mf,0,0,destcomp,numcomp);
        }
        else
        {
            FabSet fs(bndry[face].boxArray(),numcomp);

            fs.setVal(0);

            fs.copyFrom(mf,0,0,0,numcomp);

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (FabSetIter mfi(fs); mfi.isValid(); ++mfi)
                bndry[face][mfi].plus(fs[mfi],0,destcomp,numcomp);
        }
    }
}

void
FluxRegister::CrseInit (const MultiFab& mflx,
                        int             dir,
                        int             srccomp,
                        int             destcomp,
                        int             numcomp,
                        Real            mult,
                        FrOp            op)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= mflx.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    MultiFab area(mflx.boxArray(), 1, mflx.nGrow());

    area.setVal(1, 0, 1, area.nGrow());

    CrseInit(mflx,area,dir,srccomp,destcomp,numcomp,mult,op);
}

void
FluxRegister::FineAdd (const MultiFab& mflx,
                       int             dir,
                       int             srccomp,
                       int             destcomp,
                       int             numcomp,
                       Real            mult)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mflx); mfi.isValid(); ++mfi)
    {
        const int k = mfi.index();
        FineAdd(mflx[mfi],dir,k,srccomp,destcomp,numcomp,mult);
    }
}

void
FluxRegister::FineAdd (const MultiFab& mflx,
                       const MultiFab& area,
                       int             dir,
                       int             srccomp,
                       int             destcomp,
                       int             numcomp,
                       Real            mult)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mflx); mfi.isValid(); ++mfi)
    {
        const int k = mfi.index();
        FineAdd(mflx[mfi],area[k],dir,k,srccomp,destcomp,numcomp,mult);
    }
}

void
FluxRegister::FineAdd (const FArrayBox& flux,
                       int              dir,
                       int              boxno,
                       int              srccomp,
                       int              destcomp,
                       int              numcomp,
                       Real             mult)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= flux.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    const Box&  flxbox = flux.box();
    const int*  flo    = flxbox.loVect();
    const int*  fhi    = flxbox.hiVect();
    const Real* flxdat = flux.dataPtr(srccomp);

    FArrayBox& loreg = bndry[Orientation(dir,Orientation::low)][boxno];

#ifndef NDEBUG
    Box cbox = BoxLib::coarsen(flux.box(),ratio);
    BL_ASSERT(cbox.contains(loreg.box()));
#endif
    const int* rlo = loreg.box().loVect();
    const int* rhi = loreg.box().hiVect();
    Real* lodat = loreg.dataPtr(destcomp);
    FORT_FRFINEADD(lodat,ARLIM(rlo),ARLIM(rhi),
                   flxdat,ARLIM(flo),ARLIM(fhi),
                   &numcomp,&dir,ratio.getVect(),&mult);

    FArrayBox& hireg = bndry[Orientation(dir,Orientation::high)][boxno];

#ifndef NDEBUG
    BL_ASSERT(cbox.contains(hireg.box()));
#endif
    rlo = hireg.box().loVect();
    rhi = hireg.box().hiVect();
    Real* hidat = hireg.dataPtr(destcomp);
    FORT_FRFINEADD(hidat,ARLIM(rlo),ARLIM(rhi),
                   flxdat,ARLIM(flo),ARLIM(fhi),
                   &numcomp,&dir,ratio.getVect(),&mult);
}

void
FluxRegister::FineAdd (const FArrayBox& flux,
                       const FArrayBox& area,
                       int              dir,
                       int              boxno,
                       int              srccomp,
                       int              destcomp,
                       int              numcomp,
                       Real             mult)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= flux.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    const Real* area_dat = area.dataPtr();
    const int*  alo      = area.loVect();
    const int*  ahi      = area.hiVect();
    const Box&  flxbox   = flux.box();
    const int*  flo      = flxbox.loVect();
    const int*  fhi      = flxbox.hiVect();
    const Real* flxdat   = flux.dataPtr(srccomp);

    FArrayBox& loreg = bndry[Orientation(dir,Orientation::low)][boxno];

#ifndef NDEBUG
    Box cbox = BoxLib::coarsen(flux.box(),ratio);
    BL_ASSERT(cbox.contains(loreg.box()));
#endif
    const int* rlo = loreg.box().loVect();
    const int* rhi = loreg.box().hiVect();
    Real* lodat = loreg.dataPtr(destcomp);
    FORT_FRFAADD(lodat,ARLIM(rlo),ARLIM(rhi),
                 flxdat,ARLIM(flo),ARLIM(fhi),
                 area_dat,ARLIM(alo),ARLIM(ahi),
                 &numcomp,&dir,ratio.getVect(),&mult);

    FArrayBox& hireg = bndry[Orientation(dir,Orientation::high)][boxno];

#ifndef NDEBUG
    BL_ASSERT(cbox.contains(hireg.box()));
#endif
    rlo = hireg.box().loVect();
    rhi = hireg.box().hiVect();
    Real* hidat = hireg.dataPtr(destcomp);
    FORT_FRFAADD(hidat,ARLIM(rlo),ARLIM(rhi),
                 flxdat,ARLIM(flo),ARLIM(fhi),
                 area_dat,ARLIM(alo),ARLIM(ahi),
                 &numcomp,&dir,ratio.getVect(),&mult);
}

void 
FluxRegister::Reflux (MultiFab&       mf,
		      const MultiFab& volume,
		      Real            scale,
		      int             scomp,
		      int             dcomp,
		      int             ncomp,
		      const Geometry& geom)
{
    BL_PROFILE("FluxRegister::Reflux()");

    for (OrientationIter fi; fi; ++fi)
    {
	const Orientation& face = fi();
	int idir = face.coordDir();
	int islo = face.isLow();

	MultiFab flux(mf.boxArray(), ncomp, 0, Fab_allocate, IntVect::TheDimensionVector(idir));
	flux.setVal(0.0);

	flux.copy(bndry[face], scomp, 0, ncomp, 0, 0);
	geom.PeriodicCopy(flux, bndry[face], 0, scomp, ncomp);

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();

	    FArrayBox& sfab = mf[mfi];
	    const Box& sbox = sfab.box();

	    const FArrayBox& ffab = flux[mfi];
	    const Box& fbox = ffab.box();
	    
	    const FArrayBox& vfab = volume[mfi];
	    const Box& vbox = vfab.box();

	    FORT_FRREFLUX(bx.loVect(), bx.hiVect(),
			  sfab.dataPtr(dcomp), sbox.loVect(), sbox.hiVect(),
			  ffab.dataPtr(     ), fbox.loVect(), fbox.hiVect(),
			  vfab.dataPtr(     ), vfab.loVect(), vbox.hiVect(),
			  &ncomp, &scale, &idir, &islo);
			  
	}
    }
}

void 
FluxRegister::Reflux (MultiFab&       mf,
		      Real            scale,
		      int             scomp,
		      int             dcomp,
		      int             ncomp,
		      const Geometry& geom)
{
    const Real* dx = geom.CellSize();
    
    MultiFab volume(mf.boxArray(), 1, mf.nGrow());
    
    volume.setVal(D_TERM(dx[0],*dx[1],*dx[2]), 0, 1, mf.nGrow());

    Reflux(mf,volume,scale,scomp,dcomp,ncomp,geom);
}


void
FluxRegister::write (const std::string& name, std::ostream& os) const
{
    if (ParallelDescriptor::IOProcessor())
    {
        os << ratio      << '\n';
        os << fine_level << '\n';
        os << ncomp      << '\n';
    }

    const BndryRegister* br = this;

    br->write(name,os);
}


void
FluxRegister::read (const std::string& name, std::istream& is)
{

    is >> ratio;
    is >> fine_level;
    is >> ncomp;

    BndryRegister* br = this;

    br->read(name,is);
}
