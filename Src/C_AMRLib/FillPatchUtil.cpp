
#include <FillPatchUtil.H>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace BoxLib
{
    bool ProperlyNested (const IntVect& ratio, int blocking_factor, int ngrow,
			 const IndexType& boxType, Interpolater* mapper)
    {
	int ratio_max = ratio[0];
#if (BL_SPACEDIM > 1)
	ratio_max = std::max(ratio_max, ratio[1]);
#endif
#if (BL_SPACEDIM == 3)
	ratio_max = std::max(ratio_max, ratio[2]);
#endif
	// There are at least this many coarse cells outside fine grids 
	// (except at physical boundaries).
	int nbuf = blocking_factor / ratio_max;
	
	Box crse_box(IntVect(D_DECL(0 ,0 ,0 )), IntVect(D_DECL(4*nbuf-1,4*nbuf-1,4*nbuf-1)));
	crse_box.convert(boxType);
	Box fine_box(IntVect(D_DECL(  nbuf  ,  nbuf  ,  nbuf)),
		     IntVect(D_DECL(3*nbuf-1,3*nbuf-1,3*nbuf-1)));
	fine_box.convert(boxType);
	fine_box.refine(ratio_max);
	fine_box.grow(ngrow);
	const Box& fine_box_coarsened = mapper->CoarseBox(fine_box, ratio_max);
	return crse_box.contains(fine_box_coarsened);
    }

    void FillPatchSingleLevel (MultiFab& mf, Real time, 
			       const PArray<MultiFab>& smf, const std::vector<Real>& stime,
			       int scomp, int dcomp, int ncomp,
			       const Geometry& geom, PhysBCFunct& physbcf)
    {
	BL_PROFILE("FillPatchSingleLevel");

	BL_ASSERT(scomp+ncomp <= smf[0].nComp());
	BL_ASSERT(dcomp+ncomp <= mf.nComp());
	BL_ASSERT(smf.size() == stime.size());
	BL_ASSERT(smf.size() != 0);

	if (smf.size() == 1) 
	{
	    mf.copy(smf[0], scomp, dcomp, ncomp, 0, mf.nGrow());
	    geom.PeriodicCopy(mf, smf[0], dcomp, scomp, ncomp, mf.nGrow());
	} 
	else if (smf.size() == 2) 
	{
	    BL_ASSERT(smf[0].boxArray() == smf[1].boxArray());
	    PArray<MultiFab> raii(PArrayManage);
	    MultiFab * dmf;
	    int destcomp;
	    bool sameba;
	    if (mf.boxArray() == smf[0].boxArray()) {
		dmf = &mf;
		destcomp = dcomp;
		sameba = true;
	    } else {
		dmf = raii.push_back(new MultiFab(smf[0].boxArray(), ncomp, 0));
		destcomp = 0;
		sameba = false;
	    }

#ifdef _OPENMP
#pragma omp parallel 
#endif
	    for (MFIter mfi(*dmf,true); mfi.isValid(); ++mfi)
	    {
		const Box& bx = mfi.tilebox();
		(*dmf)[mfi].linInterp(smf[0][mfi],
				      scomp,
				      smf[1][mfi],
				      scomp,
				      stime[0],
				      stime[1],
				      time,
				      bx,
				      destcomp,
				      ncomp);
	    }
	    
	    if (sameba)
	    {
		mf.FillBoundary_nowait(dcomp,ncomp);
		geom.FillPeriodicBoundary_nowait(mf,dcomp,ncomp);

		mf.FillBoundary_finish();
		geom.FillPeriodicBoundary_finish(mf);
	    }
	    else
	    {
		int src_ngrow = 0;
		int dst_ngrow = mf.nGrow();

		mf.copy(*dmf, 0, dcomp, ncomp, src_ngrow, dst_ngrow);
		geom.PeriodicCopy(mf, *dmf, dcomp, 0, ncomp, dst_ngrow);
	    }
	}
	else {
	    BoxLib::Abort("FillPatchSingleLevel: high-order interpolation in time not implemented yet");
	}

	physbcf.doit(mf, dcomp, ncomp, time);
    }


    void FillPatchTwoLevels (MultiFab& mf, Real time,
			     const PArray<MultiFab>& cmf, const std::vector<Real>& ct,
			     const PArray<MultiFab>& fmf, const std::vector<Real>& ft,
			     int scomp, int dcomp, int ncomp,
			     const Geometry& cgeom, const Geometry& fgeom, 
			     PhysBCFunct& cbc, PhysBCFunct& fbc,
			     const IntVect& ratio, 
			     Interpolater* mapper, const Array<BCRec>& bcs)
    {
	BL_PROFILE("FillPatchTwoLevels");

	int ngrow = mf.nGrow();
	    
	if (ngrow > 0 || mf.getBDKey() != fmf[0].getBDKey()) 
	{
	    BoxArray ba_crse_patch;
	    DistributionMapping dm_crse_patch;
	    InterpolaterBoxCoarsener coarsener = mapper->BoxCoarsener(ratio);
	    
	    Box fdomain = fgeom.Domain();
	    fdomain.convert(mf.boxArray().ixType());
	    Box fdomain_g(fdomain);
	    for (int i = 0; i < BL_SPACEDIM; ++i) {
		if (fgeom.isPeriodic(i)) {
		    fdomain_g.grow(i,ngrow);
		}
	    }
	    
	    const FabArrayBase::FPC& fpc = FabArrayBase::TheFPC(fmf[0], mf, fdomain_g, ngrow, coarsener);

	    if (!fpc.ba_crse_patch.empty())
	    {
		MultiFab mf_crse_patch(fpc.ba_crse_patch, ncomp, 0, fpc.dm_crse_patch);
		
		FillPatchSingleLevel(mf_crse_patch, time, cmf, ct, scomp, 0, ncomp, cgeom, cbc);
		
		int idummy1=0, idummy2=0;
		bool cc = fpc.ba_crse_patch.ixType().cellCentered();
#ifdef _OPENMP
#pragma omp parallel if (cc)
#endif
		for (MFIter mfi(mf_crse_patch); mfi.isValid(); ++mfi)
		{
		    int li = mfi.LocalIndex();
		    int gi = fpc.dst_idxs[li];		
		    const Box& dbx = fpc.dst_boxes[li];
		    
		    Array<BCRec> bcr(ncomp);
		    BoxLib::setBC(dbx,fdomain,scomp,0,ncomp,bcs,bcr);
		    
		    mapper->interp(mf_crse_patch[mfi],
				   0,
				   mf[gi],
				   dcomp,
				   ncomp,
				   dbx,
				   ratio,
				   cgeom,
				   fgeom,
				   bcr,
				   idummy1, idummy2);
		}
	    }
	}

	FillPatchSingleLevel(mf, time, fmf, ft, scomp, dcomp, ncomp, fgeom, fbc);
    }
}
