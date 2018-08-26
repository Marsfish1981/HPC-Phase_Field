
#include <winstd.H>

#include <iostream>

#include <BoxArray.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <MultiFab.H>
#include <FArrayBox.H>
#include <BLProfiler.H>
#include <Utility.H>
#include <SPACE.H>

#ifdef BL_LAZY
#include <Lazy.H>
#endif

#ifdef BL_MEM_PROFILING
#include <MemProfiler.H>
long Geometry::fpb_cache_total_bytes     = 0L;
long Geometry::fpb_cache_total_bytes_hwm = 0L;
#endif

//
// The definition of some static data members.
//
int     Geometry::spherical_origin_fix = 0;
RealBox Geometry::prob_domain;
bool    Geometry::is_periodic[BL_SPACEDIM] = {D_DECL(0,0,0)};
int     Geometry::fpb_cache_max_size = 25;
Geometry::FPBMMap Geometry::m_FPBCache;

std::ostream&
operator<< (std::ostream&   os,
            const Geometry& g)
{
    os << (CoordSys&) g << g.ProbDomain() << g.Domain();
    return os;
}

std::istream&
operator>> (std::istream& is,
            Geometry&     g)
{
    Box     bx;
    RealBox rb;

    is >> (CoordSys&) g >> rb >> bx;

    g.Domain(bx);
    Geometry::ProbDomain(rb);

    return is;
}

Geometry::FPB::FPB ()
    :
    m_ngrow(-1),
    m_do_corners(false),
    m_reused(false),
    m_threadsafe_loc(false),
    m_threadsafe_rcv(false),
    m_LocTags(0),
    m_SndTags(0),
    m_RcvTags(0),
    m_SndVols(0),
    m_RcvVols(0) {}

Geometry::FPB::FPB (const BoxArray&            ba,
                    const DistributionMapping& dm,
                    const Box&                 domain,
                    int                        ngrow,
                    bool                       do_corners)
    :
    m_ba(ba),
    m_dm(dm),
    m_domain(domain),
    m_ngrow(ngrow),
    m_do_corners(do_corners),
    m_reused(false),
    m_threadsafe_loc(false),
    m_threadsafe_rcv(false),
    m_LocTags(0),
    m_SndTags(0),
    m_RcvTags(0),
    m_SndVols(0),
    m_RcvVols(0)
{
    BL_ASSERT(ngrow >= 0);
    BL_ASSERT(domain.ok());
}

Geometry::FPB::~FPB ()
{
    delete m_LocTags;
    delete m_SndTags;
    delete m_RcvTags;
    delete m_SndVols;
    delete m_RcvVols;
}

bool
Geometry::FPB::operator== (const FPB& rhs) const
{
    return
        m_ngrow == rhs.m_ngrow && m_do_corners == rhs.m_do_corners && m_domain == rhs.m_domain && m_ba == rhs.m_ba && m_dm == rhs.m_dm;
}

long 
Geometry::FPB::bytesOfMapOfComTagContainers (const Geometry::FPB::MapOfFPBComTagContainers& m)
{
    long r = sizeof(MapOfFPBComTagContainers);
    for (FPB::MapOfFPBComTagContainers::const_iterator it = m.begin(); it != m.end(); ++it) {
	r += sizeof(it->first) + BoxLib::bytesOf(it->second)
	    + BoxLib::gcc_map_node_extra_bytes;
    }
    return r;
}

long
Geometry::FPB::bytes () const
{
    long cnt = sizeof(Geometry::FPB);

    if (m_LocTags)
        cnt += BoxLib::bytesOf(*m_LocTags);

    if (m_SndTags)
	cnt += bytesOfMapOfComTagContainers(*m_SndTags);

    if (m_RcvTags)
        cnt += bytesOfMapOfComTagContainers(*m_RcvTags);

    if (m_SndVols)
	cnt += BoxLib::bytesOf(*m_SndVols);

    if (m_RcvVols)
	cnt += BoxLib::bytesOf(*m_RcvVols);

    return cnt;
}

void
Geometry::FillPeriodicBoundary (MultiFab& mf,
                                bool      do_corners,
                                bool      local) const
{
    FillPeriodicBoundary(mf,0,mf.nComp(),do_corners,local);
}

void
Geometry::FillPeriodicBoundary_nowait (MultiFab& mf,
				       bool      do_corners,
				       bool      local) const
{
    FillPeriodicBoundary_nowait(mf,0,mf.nComp(),do_corners,local);
}

void
Geometry::SumPeriodicBoundary (MultiFab& mf) const
{
    SumPeriodicBoundary(mf,0,mf.nComp());
}

void
Geometry::FillPeriodicBoundary (MultiFab& mf,
                                int       scomp,
                                int       ncomp,
                                bool      corners,
                                bool      local) const
{
    BL_PROFILE("Geometry::FillPeriodicBoundary()");

    if ( local )
    {
	FillPeriodicBoundary_local(mf, scomp, ncomp, corners);
    }
    else
    {
	if (!isAnyPeriodic() || mf.nGrow() == 0 || mf.size() == 0) return;
        BoxLib::FillPeriodicBoundary_nowait(*this, mf, scomp, ncomp, corners);
	BoxLib::FillPeriodicBoundary_finish(*this, mf);
    }
}

void
Geometry::FillPeriodicBoundary_local (MultiFab& mf,
				      int       scomp,
				      int       ncomp,
				      bool      corners) const
{
    if (!isAnyPeriodic() || mf.nGrow() == 0 || mf.size() == 0) return;

    BL_PROFILE("Geometry::FillPeriodicBoundary_local()");

    //
    // Do what you can with the FABs you own.  No parallelism allowed.
    //
    
    Box TheDomain = Domain();
    TheDomain.convert(mf.boxArray().ixType());

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<IntVect> pshifts(26);

        for (MFIter mfidst(mf); mfidst.isValid(); ++mfidst)
        {
            const Box& dst = mf[mfidst].box();

            BL_ASSERT(dst == BoxLib::grow(mfidst.validbox(), mf.nGrow()));

            if (TheDomain.contains(dst)) continue;

            // Turn off sharing among threads because this MFIter is inside another MFIter
	    unsigned char flags = MFIter::AllBoxes || MFIter::NoTeamBarrier;
            for (MFIter mfisrc(mf,flags); mfisrc.isValid(); ++mfisrc)
            {
                Box src = mfisrc.validbox() & TheDomain;

                if (corners)
                {
                    for (int i = 0; i < BL_SPACEDIM; i++)
                    {
                        if (!isPeriodic(i))
                        {
                            if (src.smallEnd(i) == Domain().smallEnd(i))
                                src.growLo(i,mf.nGrow());
                            if (src.bigEnd(i) == Domain().bigEnd(i))
                                src.growHi(i,mf.nGrow());
                        }
                    }
                }

                if (TheDomain.contains(src)) continue;

                periodicShift(dst, src, pshifts);

                for (Array<IntVect>::const_iterator it = pshifts.begin(), End = pshifts.end();
                     it != End;
                     ++it)
                {
                    const IntVect& iv = *it;

                    const Box& shft = src + iv;
                    const Box& dbx  = dst & shft;
                    const Box& sbx  = dbx - iv;

                    mf[mfidst].copy(mf[mfisrc],sbx,scomp,dbx,scomp,ncomp);
                }
            }
        }
    }
}

void
Geometry::FillPeriodicBoundary_nowait (MultiFab& mf,
				       int       scomp,
				       int       ncomp,
				       bool      corners,
				       bool      local) const
{
    BL_PROFILE("Geometry::FillPeriodicBoundary_nowait()");

    if ( local )
    {
	FillPeriodicBoundary_local(mf, scomp, ncomp, corners);
    }
    else
    {
	if (!isAnyPeriodic() || mf.nGrow() == 0 || mf.size() == 0) return;
        BoxLib::FillPeriodicBoundary_nowait(*this, mf, scomp, ncomp, corners);
    }
}

void
Geometry::FillPeriodicBoundary_finish (MultiFab& mf) const
{
    if (!isAnyPeriodic() || mf.nGrow() == 0 || mf.size() == 0) return;
    BoxLib::FillPeriodicBoundary_finish(*this, mf);
}

void
Geometry::SumPeriodicBoundary (MultiFab& mf,
                               int       scomp,
                               int       ncomp) const
{
    BL_PROFILE("SumPeriodicBoundaryInnards()");

    if (!isAnyPeriodic() || mf.nGrow() == 0 || mf.size() == 0) return;

    BL_ASSERT( mf.boxArray().ixType().cellCentered() ); 

    //
    // We're going to attempt to reuse the information in the FPB
    // cache.  The intersection info should be that same.  It's what we do
    // with it that's different.  Basically we have to send the m_RcvTags and
    // receive the m_SndTags.
    //
    const Geometry& geom = *this;

    const Geometry::FPB fpb(mf.boxArray(),mf.DistributionMap(),geom.Domain(),mf.nGrow(),false);

    Geometry::FPBMMapIter cache_it = Geometry::GetFPB(geom,fpb,mf);

    BL_ASSERT(cache_it != Geometry::m_FPBCache.end());

    const Geometry::FPB& TheFPB_ = cache_it->second;

    const FPB::FPBComTagsContainer&      LocTags = *(TheFPB_.m_LocTags);
    const FPB::MapOfFPBComTagContainers& SndTags = *(TheFPB_.m_RcvTags);
    const FPB::MapOfFPBComTagContainers& RcvTags = *(TheFPB_.m_SndTags);
    const std::map<int,int>&             SndVols = *(TheFPB_.m_RcvVols);
    const std::map<int,int>&             RcvVols = *(TheFPB_.m_SndVols);

    if (ParallelDescriptor::NProcs() == 1)
    {
	//
	// There can only be local work to do.
	//
	int N_loc = LocTags.size();
	// undafe to do OMP
	for (int i=0; i<N_loc; ++i)
        {
	    const Geometry::FPBComTag& tag = LocTags[i];
	    const Box& srcbox = tag.dbox;   // Note the switch of dst and src here.
	    const Box& dstbox = tag.sbox;
            mf[tag.srcIndex].plus(mf[tag.dstIndex],srcbox,dstbox,scomp,scomp,ncomp);
	}

        return;
    }

#ifdef BL_USE_MPI
    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    int SeqNum;
    {
	ParallelDescriptor::Color mycolor = mf.color();
	if (mycolor == ParallelDescriptor::DefaultColor()) {
	    SeqNum = ParallelDescriptor::SeqNum();
	} else if (mycolor == ParallelDescriptor::SubCommColor()) {
	    SeqNum = ParallelDescriptor::SubSeqNum();
	}
	// else I don't have any data and my SubSeqNum() should not be called.
    }

    if (LocTags.empty() && RcvTags.empty() && SndTags.empty())
	//
	// No work to do.
	//
	return;

    Array<MPI_Status>  stats;
    Array<int>         recv_from;
    Array<Real*>       recv_data;
    Array<MPI_Request> recv_reqs;
    //
    // Post rcvs. Allocate one chunk of space to hold'm all.
    //
    Real* the_recv_data = 0;

    FabArrayBase::PostRcvs(RcvVols,the_recv_data,recv_data,recv_from,recv_reqs,ncomp,SeqNum);

    //
    // Post send's
    //
    const int N_snds = SndTags.size();

    Array<Real*>                           send_data;
    Array<int>                             send_N;
    Array<int>                             send_rank;
    Array<const FPB::FPBComTagsContainer*> send_cctc;

    send_data.reserve(N_snds);
    send_N   .reserve(N_snds);
    send_rank.reserve(N_snds);
    send_cctc.reserve(N_snds);

    for (FPB::MapOfFPBComTagContainers::const_iterator m_it = SndTags.begin(),
             m_End = SndTags.end();
         m_it != m_End;
         ++m_it)
    {
        std::map<int,int>::const_iterator vol_it = SndVols.find(m_it->first);

        BL_ASSERT(vol_it != SndVols.end());

        const int N = vol_it->second*ncomp;

        BL_ASSERT(N < std::numeric_limits<int>::max());

        Real* data = static_cast<Real*>(BoxLib::The_Arena()->alloc(N*sizeof(Real)));

	send_data.push_back(data);
	send_N   .push_back(N);
	send_rank.push_back(m_it->first);
	send_cctc.push_back(&(m_it->second));
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i=0; i<N_snds; ++i)
    {
	Real*dptr = send_data[i];
	BL_ASSERT(dptr != 0);

	const FPB::FPBComTagsContainer& cctc = *send_cctc[i];

        for (FPB::FPBComTagsContainer::const_iterator it = cctc.begin();
             it != cctc.end(); ++it)
        {
            BL_ASSERT(mf.DistributionMap()[it->dstIndex] == ParallelDescriptor::MyProc());
            const Box& dbx = it->dbox;
            mf[it->dstIndex].copyToMem(dbx,scomp,ncomp,dptr);
            const int Cnt = dbx.numPts()*ncomp;
            dptr += Cnt;
        }
    }

    Array<MPI_Request> send_reqs;

    if (FabArrayBase::do_async_sends)
    {
	send_reqs.reserve(N_snds);
	for (int i=0; i<N_snds; ++i) {
            send_reqs.push_back(ParallelDescriptor::Asend
				(send_data[i],send_N[i],send_rank[i],SeqNum).req());
        }
    } else {
	for (int i=0; i<N_snds; ++i) {
            ParallelDescriptor::Send(send_data[i],send_N[i],send_rank[i],SeqNum);
            BoxLib::The_Arena()->free(send_data[i]);
        }
    }

    //
    // Do the local work.  Hope for a bit of communication/computation overlap.
    //
    int N_loc = LocTags.size();
    // undafe to do OMP
    for (int i=0; i<N_loc; ++i)
    {
        const FPBComTag& tag = LocTags[i];

        BL_ASSERT(mf.DistributionMap()[tag.dstIndex] == ParallelDescriptor::MyProc());
        BL_ASSERT(mf.DistributionMap()[tag.srcIndex] == ParallelDescriptor::MyProc());

	const Box& srcbox = tag.dbox;   // Note the switch of dst and src here.
	const Box& dstbox = tag.sbox;
        mf[tag.srcIndex].plus(mf[tag.dstIndex],srcbox,dstbox,scomp,scomp,ncomp);
    }

    //
    //  wait and unpack
    //
    const int N_rcvs = RcvTags.size();

    if (N_rcvs > 0)
    {
	Array<const FPB::FPBComTagsContainer*> recv_cctc;
	recv_cctc.reserve(N_rcvs);

        for (int k = 0; k < N_rcvs; k++)
        {
	    FPB::MapOfFPBComTagContainers::const_iterator m_it = RcvTags.find(recv_from[k]);
            BL_ASSERT(m_it != RcvTags.end());

	    recv_cctc.push_back(&(m_it->second));
	}

	stats.resize(N_rcvs);
	BL_MPI_REQUIRE( MPI_Waitall(N_rcvs, recv_reqs.dataPtr(), stats.dataPtr()) );

	// unsafe to do OMP
	{
	    FArrayBox fab;

	    for (int k = 0; k < N_rcvs; k++) 
	    {
		const Real* dptr = recv_data[k];
		BL_ASSERT(dptr != 0);
		
		const FPB::FPBComTagsContainer& cctc = *recv_cctc[k];
		
		for (FPB::FPBComTagsContainer::const_iterator it = cctc.begin();
		     it != cctc.end(); ++it)
		{
		    BL_ASSERT(mf.DistributionMap()[it->srcIndex] == ParallelDescriptor::MyProc());
		    const Box& srcbox = it->dbox;   // Note the switch of dst and src here.
		    const Box& dstbox = it->sbox;
		    fab.resize(srcbox,ncomp);
		    const int Cnt = dstbox.numPts()*ncomp;
		    memcpy(fab.dataPtr(), dptr, Cnt*sizeof(Real));
		    mf[it->srcIndex].plus(fab,srcbox,dstbox,0,scomp,ncomp);
		    dptr += Cnt;
		}
	    }
        }
    }

    BoxLib::The_Arena()->free(the_recv_data);

    if (FabArrayBase::do_async_sends && !SndTags.empty())
        FabArrayBase::GrokAsyncSends(N_snds,send_reqs,send_data,stats);

#endif

}

void
Geometry::PeriodicCopy (MultiFab&       dstmf,
			const MultiFab& srcmf) const
{
    PeriodicCopy(dstmf, srcmf, 0, 0, srcmf.nComp());
}

void 
Geometry::PeriodicCopy (MultiFab&       dstmf,
			const MultiFab& srcmf,
			int             dcomp,
			int             scomp,
			int             ncomp,
			int             dstng,
			int             srcng) const
{
    BoxLib::PeriodicCopy(*this, dstmf, srcmf, dcomp, scomp, ncomp, dstng, srcng,
			 FabArrayBase::COPY);
}

Geometry::Geometry () {}

Geometry::Geometry (const Box&     dom,
                    const RealBox* rb,
                    int            coord,
                    int*           is_per)
{
    define(dom,rb,coord,is_per);
}

Geometry::Geometry (const Geometry& g)
{
    ok     = g.ok;
    domain = g.domain;

    D_TERM(dx[0]=g.dx[0];,dx[1]=g.dx[1];,dx[2]=g.dx[2];)
}

Geometry::~Geometry() {}

void
Geometry::define (const Box&     dom,
                  const RealBox* rb,
                  int            coord,
                  int*           is_per)
{
    if (c_sys == undef)
        Setup(rb,coord,is_per);

    domain = dom;
    ok     = true;

    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        dx[k] = prob_domain.length(k)/(Real(domain.length(k)));
    }
    if (Geometry::spherical_origin_fix == 1)
    {
	if (c_sys == SPHERICAL && prob_domain.lo(0) == 0 && BL_SPACEDIM > 1)
        {
            prob_domain.setLo(0,2*dx[0]);

            for (int k = 0; k < BL_SPACEDIM; k++)
            {
                dx[k] = prob_domain.length(k)/(Real(domain.length(k)));
            }
	}
    } 
}

void
Geometry::Finalize ()
{
    c_sys = undef;

    Geometry::FlushPIRMCache();
}

void
Geometry::Setup (const RealBox* rb, int coord, int* isper)
{
    ParmParse pp("geometry");
    //
    // The default behavior is as before.  If rb and coord come
    // in with default values, we require that user set them through pp.
    // If not, use those coming in, and possibly override them w/pp
    //
    Array<Real> prob_lo(BL_SPACEDIM);
    Array<Real> prob_hi(BL_SPACEDIM);
    if (rb == 0  &&  coord==-1)
    {
        pp.get("coord_sys",coord);
        SetCoord( (CoordType) coord );
        pp.getarr("prob_lo",prob_lo,0,BL_SPACEDIM);
        BL_ASSERT(prob_lo.size() == BL_SPACEDIM);
        pp.getarr("prob_hi",prob_hi,0,BL_SPACEDIM);
        BL_ASSERT(prob_lo.size() == BL_SPACEDIM);
        prob_domain.setLo(prob_lo);
        prob_domain.setHi(prob_hi);
    }
    else
    {
        BL_ASSERT(rb != 0  &&  coord != -1);
        pp.query("coord_sys",coord);
        SetCoord( (CoordType) coord );
        prob_domain.setLo(rb->lo());
        prob_domain.setHi(rb->hi());

        if (pp.countval("prob_lo")>0)
        {
            pp.queryarr("prob_lo",prob_lo,0,BL_SPACEDIM);
            BL_ASSERT(prob_lo.size() == BL_SPACEDIM);
            prob_domain.setLo(prob_lo);
        }
        if (pp.countval("prob_hi")>0)
        {
            pp.queryarr("prob_hi",prob_hi,0,BL_SPACEDIM);
            BL_ASSERT(prob_hi.size() == BL_SPACEDIM);
            prob_domain.setHi(prob_hi);
        }
    }
    pp.query("spherical_origin_fix", Geometry::spherical_origin_fix);
    pp.query("fpb_cache_max_size",   Geometry::fpb_cache_max_size);
    //
    // Don't let the cache size get too small.  This simplifies some logic later.
    //
    if (Geometry::fpb_cache_max_size < 1)
        Geometry::fpb_cache_max_size = 1;
    //
    // Now get periodicity info.
    //
    if (isper == 0)
    {
        Array<int> is_per(BL_SPACEDIM);
        pp.queryarr("is_periodic",is_per,0,BL_SPACEDIM);
        for (int n = 0; n < BL_SPACEDIM; n++)  
            is_periodic[n] = is_per[n];
    }
    else
    {
        for (int n = 0; n < BL_SPACEDIM; n++)  
            is_periodic[n] = isper[n];
    }

#ifdef BL_MEM_PROFILING
    MemProfiler::add("Geometry", [] () -> MemProfiler::MemInfo {
	    return {fpb_cache_total_bytes, fpb_cache_total_bytes_hwm};
	});
#endif

    BoxLib::ExecOnFinalize(Geometry::Finalize);
}

void
Geometry::GetVolume (MultiFab&       vol,
                     const BoxArray& grds,
                     int             ngrow) const
{
    vol.define(grds,1,ngrow,Fab_allocate);
    GetVolume(vol);
}

void
Geometry::GetVolume (MultiFab&       vol) const
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(vol,true); mfi.isValid(); ++mfi)
    {
	CoordSys::SetVolume(vol[mfi], mfi.growntilebox());
    }
}

void
Geometry::GetVolume (FArrayBox&      vol,
                     const BoxArray& grds,
                     int             idx,
                     int             ngrow) const
{
    CoordSys::GetVolume(vol, BoxLib::grow(grds[idx],ngrow));
}

#if (BL_SPACEDIM <= 2)
void
Geometry::GetDLogA (MultiFab&       dloga,
                    const BoxArray& grds, 
                    int             dir,
                    int             ngrow) const
{
    dloga.define(grds,1,ngrow,Fab_allocate);
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(dloga,true); mfi.isValid(); ++mfi)
    {
	CoordSys::SetDLogA(dloga[mfi], mfi.growntilebox(), dir);
    }
}
#endif

void
Geometry::GetFaceArea (MultiFab&       area,
                       const BoxArray& grds,
                       int             dir,
                       int             ngrow) const
{
    BoxArray edge_boxes(grds);
    edge_boxes.surroundingNodes(dir);
    area.define(edge_boxes,1,ngrow,Fab_allocate);

    GetFaceArea(area, dir);
}

void
Geometry::GetFaceArea (MultiFab&       area,
                       int             dir) const
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(area,true); mfi.isValid(); ++mfi)
    {
	CoordSys::SetFaceArea(area[mfi],mfi.growntilebox(),dir);
    }
}

void
Geometry::GetFaceArea (FArrayBox&      area,
                       const BoxArray& grds,
                       int             idx,
                       int             dir,
                       int             ngrow) const
{
    CoordSys::GetFaceArea(area, BoxLib::grow(grds[idx],ngrow), dir);
}

void
Geometry::periodicShift (const Box&      target,
                         const Box&      src, 
                         Array<IntVect>& out) const
{
    out.resize(0);

    Box locsrc(src);

    int nist,njst,nkst;
    int niend,njend,nkend;
    nist = njst = nkst = 0;
    niend = njend = nkend = 0;
    D_TERM( nist , =njst , =nkst ) = -1;
    D_TERM( niend , =njend , =nkend ) = +1;

    int ri,rj,rk;
    for (ri = nist; ri <= niend; ri++)
    {
        if (ri != 0 && !is_periodic[0])
            continue;
        if (ri != 0 && is_periodic[0])
            locsrc.shift(0,ri*domain.length(0));

        for (rj = njst; rj <= njend; rj++)
        {
            if (rj != 0 && !is_periodic[1])
                continue;
            if (rj != 0 && is_periodic[1])
                locsrc.shift(1,rj*domain.length(1));

            for (rk = nkst; rk <= nkend; rk++)
            {
                if (rk!=0
#if (BL_SPACEDIM == 3)
                    && !is_periodic[2]
#endif
                    )
                {
                    continue;
                }
                if (rk!=0
#if (BL_SPACEDIM == 3)
                    && is_periodic[2]
#endif
                    )
                {
                    locsrc.shift(2,rk*domain.length(2));
                }

                if (ri == 0 && rj == 0 && rk == 0)
                    continue;
                //
                // If losrc intersects target, then add to "out".
                //
                if (target.intersects(locsrc))
                {
                    out.push_back(IntVect(D_DECL(ri*domain.length(0),
                                                 rj*domain.length(1),
                                                 rk*domain.length(2))));
                }
                if (rk != 0
#if (BL_SPACEDIM == 3)
                    && is_periodic[2]
#endif
                    )
                {
                    locsrc.shift(2,-rk*domain.length(2));
                }
            }
            if (rj != 0 && is_periodic[1])
                locsrc.shift(1,-rj*domain.length(1));
        }
        if (ri != 0 && is_periodic[0])
            locsrc.shift(0,-ri*domain.length(0));
    }
}

Geometry::FPBMMapIter
Geometry::GetFPB (const Geometry&      geom,
                  const Geometry::FPB& fpb,
                  const FabArrayBase&  mf)
{
    BL_PROFILE("Geometry::GetFPB");

    BL_ASSERT(fpb.m_ngrow > 0);
    BL_ASSERT(fpb.m_ba.size() > 0);
    BL_ASSERT(geom.isAnyPeriodic());

    const BoxArray&            ba     = fpb.m_ba;
    const DistributionMapping& dm     = fpb.m_dm;
    const int                  MyProc = ParallelDescriptor::MyProc();
    const IntVect&             Typ    = ba[0].type();
    const int                  Scale  = D_TERM(Typ[0],+3*Typ[1],+5*Typ[2]) + 11;
    const int                  Key    = ba.size() + ba[0].numPts() + Scale + fpb.m_ngrow;

    std::pair<Geometry::FPBMMapIter,Geometry::FPBMMapIter> er_it = m_FPBCache.equal_range(Key);
    
    for (Geometry::FPBMMapIter it = er_it.first; it != er_it.second; ++it)
    {
        if (it->second == fpb)
        {
            it->second.m_reused = true;

            return it;
        }
    }

    if (m_FPBCache.size() >= Geometry::fpb_cache_max_size)
    {
        //
        // Don't let the size of the cache get too big.
        // Get rid of entries with the biggest largest key that haven't been reused.
        // Otherwise just remove the entry with the largest key.
        //
        Geometry::FPBMMapIter End      = m_FPBCache.end();
        Geometry::FPBMMapIter last_it  = End;
        Geometry::FPBMMapIter erase_it = End;

        for (Geometry::FPBMMapIter it = m_FPBCache.begin(); it != End; ++it)
        {
            last_it = it;

            if (!it->second.m_reused)
                erase_it = it;
        }

        if (erase_it != End)
        {
#ifdef BL_MEM_PROFILING
	    fpb_cache_total_bytes -= erase_it->second.bytes();
#endif
            m_FPBCache.erase(erase_it);
        }
        else if (last_it != End)
        {
#ifdef BL_MEM_PROFILING
	    fpb_cache_total_bytes -= last_it->second.bytes();
#endif
            m_FPBCache.erase(last_it);
        }
    }
    //
    // Got to insert one & then build it.
    //
    Geometry::FPBMMapIter cache_it = m_FPBCache.insert(FPBMMap::value_type(Key,fpb));
    FPB&                  TheFPB   = cache_it->second;
    //
    // Here's where we allocate memory for the cache innards.
    // We do this so we don't have to build objects of these types
    // each time we search the cache.  Otherwise we'd be constructing
    // and destroying said objects quite frequently.
    //
    TheFPB.m_LocTags = new FPB::FPBComTagsContainer;
    TheFPB.m_SndTags = new FPB::MapOfFPBComTagContainers;
    TheFPB.m_RcvTags = new FPB::MapOfFPBComTagContainers;
    TheFPB.m_SndVols = new std::map<int,int>;
    TheFPB.m_RcvVols = new std::map<int,int>;

    const Array<int>& imap = mf.IndexMap();

    if (imap.empty()) {
        //
        // We don't own any of the relevant FABs so can't possibly have any work to do.
        //
#ifdef BL_MEM_PROFILING
	fpb_cache_total_bytes += TheFPB.bytes();
	fpb_cache_total_bytes_hwm = std::max(fpb_cache_total_bytes_hwm,
					     fpb_cache_total_bytes);
#endif
        return cache_it;
    }

    // All workers in the same team will have identical copies of tags for local copy
    // so that they can share work.  But for remote communication, they are all different.

    const int nlocal = imap.size();
    const int ng = fpb.m_ngrow;
    std::vector<std::pair<int,Box> > isects;
    Array<IntVect> pshifts(26);

    Box TheDomain = geom.Domain();
    TheDomain.convert(ba.ixType());
    const Box& GrownDomain = BoxLib::grow(TheDomain,ng);

    FPB::MapOfFPBComTagContainers send_tags; // temp copy

    for (int i = 0; i < nlocal; ++i)
    {
	const   int k_src = imap[i];
	const Box& bx_src = ba[k_src];

	if (TheDomain.contains(BoxLib::grow(bx_src,ng))) continue;

	geom.periodicShift(GrownDomain, bx_src, pshifts);

	for (Array<IntVect>::const_iterator pit = pshifts.begin(), pEnd = pshifts.end();
	     pit != pEnd; ++pit)
	{
	    const IntVect& iv   = *pit;
	    const Box&     shft = bx_src + iv;

	    ba.intersections(shft, isects, ng);

	    for (int j = 0, M = isects.size(); j < M; ++j)
	    {
		const int  k_dst     = isects[j].first;
		const Box& bx_dst    = ba[k_dst];
		Box        bx        = isects[j].second;
		const int  dst_owner = dm[k_dst];

		if (fpb.m_do_corners) {
		    for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
			if (!geom.isPeriodic(dir)) {
			    if (bx.smallEnd(dir) == TheDomain.smallEnd(dir) 
				&& bx_dst.smallEnd(dir) == TheDomain.smallEnd(dir) )
			    {
				bx.growLo(dir,ng);
			    }
			    if (bx.bigEnd(dir) == TheDomain.bigEnd(dir)
				&& bx_dst.bigEnd(dir) == TheDomain.bigEnd(dir) )
			    {
				bx.growHi(dir,ng);
			    }
			}
		    }
		}

		if (ParallelDescriptor::sameTeam(dst_owner)) {
		    continue; // local copy will be dealt with later
		} else if (MyProc == dm[k_src]) {
		    const BoxList& bl = BoxLib::boxDiff(bx, bx_dst);
		    for (BoxList::const_iterator lit = bl.begin(); lit != bl.end(); ++lit)
			send_tags[dst_owner].push_back(FPBComTag(*lit-iv, *lit, k_src, k_dst));
		}
	    }
	}
    }

    FPB::MapOfFPBComTagContainers recv_tags; // temp copy

    BaseFab<int> localtouch, remotetouch;
    bool check_local = false, check_remote = false;
#ifdef _OPENMP
    if (omp_get_max_threads() > 1) {
        check_local = true;
        check_remote = true;
    }
#endif    

    if (ParallelDescriptor::TeamSize() > 1) {
	check_local = true;
    }

    if ( ba.ixType().cellCentered() ) {
	TheFPB.m_threadsafe_loc = true;
	TheFPB.m_threadsafe_rcv = true;
        check_local = false;
        check_remote = false;
    }

    for (int i = 0; i < nlocal; ++i)
    {
	const int   k_dst   = imap[i];
	const Box& bx_dst   = ba[k_dst];
	const Box& bx_dst_g = BoxLib::grow(bx_dst, ng);

	if (TheDomain.contains(bx_dst_g)) continue;

	if (check_local) {
	    localtouch.resize(bx_dst_g);
	    localtouch.setVal(0);
	}

	if (check_remote) {
	    remotetouch.resize(bx_dst_g);
	    remotetouch.setVal(0);
	}

	geom.periodicShift(TheDomain, bx_dst_g, pshifts);

	for (Array<IntVect>::const_iterator pit = pshifts.begin(), pEnd = pshifts.end();
	     pit != pEnd; ++pit)
	{
	    const IntVect& iv   = *pit;
	    const Box&     shft = bx_dst_g + iv;

	    ba.intersections(shft, isects);

	    for (int j = 0, M = isects.size(); j < M; ++j)
	    {
		const int k_src     = isects[j].first;
		Box       bx        = isects[j].second;
		const int src_owner = dm[k_src];

		if (fpb.m_do_corners) {
		    for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
			if (!geom.isPeriodic(dir)) {
			    if (bx.smallEnd(dir) == TheDomain.smallEnd(dir)
				&& bx_dst.smallEnd(dir) == TheDomain.smallEnd(dir) )
			    {
				bx.growLo(dir,ng);
			    }
			    if (bx.bigEnd(dir) == TheDomain.bigEnd(dir)
				&& bx_dst.bigEnd(dir) == TheDomain.bigEnd(dir) )
			    {
				bx.growHi(dir,ng);
			    }
			}
		    }
		}

		const BoxList& bl = BoxLib::boxDiff(bx-iv, bx_dst); // destinatin boxes
		for (BoxList::const_iterator lit = bl.begin(); lit != bl.end(); ++lit)
		{
		    const Box& blbx = *lit;

		    if (ParallelDescriptor::sameTeam(src_owner)) { // local copy
			const BoxList tilelist(blbx, FabArrayBase::comm_tile_size);
			for (BoxList::const_iterator
				 it_tile  = tilelist.begin(),
				 End_tile = tilelist.end();   it_tile != End_tile; ++it_tile)
			{
			    TheFPB.m_LocTags->push_back(FPBComTag((*it_tile)+iv, *it_tile,
								  k_src, k_dst));
			}
			if (check_local) {
			    localtouch.plus(1, blbx);
			}
		    } else if (MyProc == dm[k_dst]) {
			recv_tags[src_owner].push_back(FPBComTag(blbx+iv, blbx, k_src, k_dst));
			if (check_remote) {
			    remotetouch.plus(1, blbx);
			}
		    }
		}
	    }
	}

	if (check_local) {  
	    // safe if a cell is touched no more than once 
	    // keep checking thread safety if it is safe so far
            check_local = TheFPB.m_threadsafe_loc = localtouch.max() <= 1;
        }

	if (check_remote) {
            check_remote = TheFPB.m_threadsafe_rcv = remotetouch.max() <= 1;
        }
    }

    for (int ipass = 0; ipass < 2; ++ipass) // pass 0: send; pass 1: recv
    {
	FPB::MapOfFPBComTagContainers & Tags
	    = (ipass == 0) ? *TheFPB.m_SndTags : *TheFPB.m_RcvTags;
	FPB::MapOfFPBComTagContainers & tmpTags
	    = (ipass == 0) ?         send_tags :         recv_tags;
	std::map<int,int> & Vols
	    = (ipass == 0) ? *TheFPB.m_SndVols : *TheFPB.m_RcvVols;

	for (FPB::MapOfFPBComTagContainers::iterator
		 it  = tmpTags.begin(),
		 End = tmpTags.end();   it != End; ++it)
	{
	    const int key = it->first;
	    std::vector<FPBComTag>& fctv = it->second;

	    // We need to fix the order so that the send and recv processes match.
	    std::sort(fctv.begin(), fctv.end());

	    std::vector<FPBComTag> new_fctv;
	    new_fctv.reserve(fctv.size());

	    for (std::vector<FPBComTag>::const_iterator
		     it2  = fctv.begin(),
		     End2 = fctv.end();   it2 != End2; ++it2)
	    {
		const Box& sbx = it2->sbox;
		const Box& dbx = it2->dbox;
		IntVect diff = sbx.smallEnd() - dbx.smallEnd();

		Vols[key] += sbx.numPts();
		
		const BoxList tilelist(sbx, FabArrayBase::comm_tile_size);
		for (BoxList::const_iterator 
			 it_tile  = tilelist.begin(), 
			 End_tile = tilelist.end();    it_tile != End_tile; ++it_tile)
                {
		    new_fctv.push_back(FPBComTag(*it_tile, (*it_tile)-diff,
						 it2->srcIndex, it2->dstIndex));
		}
	    }

	    Tags[key].swap(new_fctv);
	} 
    }

#ifdef BL_MEM_PROFILING
    fpb_cache_total_bytes += TheFPB.bytes();
    fpb_cache_total_bytes_hwm = std::max(fpb_cache_total_bytes_hwm,
					 fpb_cache_total_bytes);
#endif

    return cache_it;
}

long
Geometry::bytesOfFPBCache ()
{
    long r;
    if (m_FPBCache.empty()) {
	r = 0L;
    } else {
	r = sizeof(m_FPBCache);
	for (FPBMMapIter it = m_FPBCache.begin(), End = m_FPBCache.end(); it != End; ++it)
	{
	    r += sizeof(it->first) + it->second.bytes() 
		+ BoxLib::gcc_map_node_extra_bytes;
	}
    }
    return r;
}

void
Geometry::FlushPIRMCache ()
{
    long stats[3] = {0,0,0}; // size, reused, bytes

    stats[0] = m_FPBCache.size();

    for (FPBMMapIter it = m_FPBCache.begin(), End = m_FPBCache.end(); it != End; ++it)
    {
        if (it->second.m_reused)
            stats[1]++;
    }

    stats[2] = bytesOfFPBCache();

    if (BoxLib::verbose)
    {
#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceLongMax(&stats[0], 3, ParallelDescriptor::IOProcessorNumber());
        if (stats[0] > 0 && ParallelDescriptor::IOProcessor())
        {
            std::cout << "Geometry::TheFPBCache: max size: "
                      << stats[0]
                      << ", max # reused: "
                      << stats[1]
                      << ", max bytes used: "
                      << stats[2]
                      << std::endl;
        }
#ifdef BL_LAZY
	});
#endif
    }

    m_FPBCache.clear();

#ifdef BL_MEM_PROFILING
    fpb_cache_total_bytes = 0L;
#endif
}

#ifdef BL_USE_MPI
void
Geometry::SendGeometryToSidecars (Geometry *geom)
{
  const int MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;

  int coord;
  int is_periodic[BL_SPACEDIM];

  if (ParallelDescriptor::Communicator() == ParallelDescriptor::CommunicatorComp())
  {
    // Data to construct base Box
    const Box& baseBox = geom->Domain();
    const IntVect box_index_type = baseBox.type();
    const int* box_index_type_IV = box_index_type.getVect();
    const int *smallEnd = baseBox.smallEnd().getVect();
    const int *bigEnd = baseBox.bigEnd().getVect();

    // Data to construct RealBox
    const RealBox& realBox = geom->ProbDomain();
    const Real *realBox_lo = realBox.lo();
    const Real *realBox_hi = realBox.hi();

    CoordType coordtype = geom->Coord();
    // UGHHH STOP USING ENUMS IN MPI CODE. IT DESTROYS LIVES.
    int coord;
    switch (coordtype)
    {
        case undef:
            coord = -1;
            break;
        case cartesian:
            coord = 0;
            break;
        case RZ:
            coord = 1;
            break;
        case SPHERICAL:
            coord = 2;
            break;
    }

    for (unsigned int i = 0; i < BL_SPACEDIM; ++i)
    {
      is_periodic[i] = geom->isPeriodic(i);
    }

      // Step 1: send the base Box
      ParallelDescriptor::Bcast(const_cast<int*>(box_index_type_IV), BL_SPACEDIM, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
      ParallelDescriptor::Bcast(const_cast<int*>(smallEnd)      , BL_SPACEDIM, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
      ParallelDescriptor::Bcast(const_cast<int*>(bigEnd)        , BL_SPACEDIM, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());

      // Step 2: send the RealBox
      ParallelDescriptor::Bcast(const_cast<Real*>(realBox_lo), BL_SPACEDIM, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
      ParallelDescriptor::Bcast(const_cast<Real*>(realBox_hi), BL_SPACEDIM, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());

      // Step 3: send the coordinates
      ParallelDescriptor::Bcast(&coord, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());

      // Step 4: send the periodicity flags
      ParallelDescriptor::Bcast(is_periodic, BL_SPACEDIM, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
  }
  else
  {
    int box_index_type[BL_SPACEDIM];
    int smallEnd[BL_SPACEDIM];
    int bigEnd[BL_SPACEDIM];

    Real realBox_lo[BL_SPACEDIM];
    Real realBox_hi[BL_SPACEDIM];

    ParallelDescriptor::Bcast(box_index_type, BL_SPACEDIM, 0, ParallelDescriptor::CommunicatorInter());
    ParallelDescriptor::Bcast(smallEnd      , BL_SPACEDIM, 0, ParallelDescriptor::CommunicatorInter());
    ParallelDescriptor::Bcast(bigEnd        , BL_SPACEDIM, 0, ParallelDescriptor::CommunicatorInter());

    ParallelDescriptor::Bcast(realBox_lo, BL_SPACEDIM, 0, ParallelDescriptor::CommunicatorInter());
    ParallelDescriptor::Bcast(realBox_hi, BL_SPACEDIM, 0, ParallelDescriptor::CommunicatorInter());

    ParallelDescriptor::Bcast(&coord, 1, 0, ParallelDescriptor::CommunicatorInter());

    ParallelDescriptor::Bcast(is_periodic, BL_SPACEDIM, 0, ParallelDescriptor::CommunicatorInter());

    // Now reconstruct all the parts
    IntVect smallEnd_IV(smallEnd);
    IntVect bigEnd_IV(bigEnd);
    IntVect box_index_type_IV(box_index_type);
    Box baseBox(smallEnd_IV, bigEnd_IV, box_index_type_IV);

    RealBox realBox;
    for (int n = 0; n < BL_SPACEDIM; n++) {
      realBox.setLo(n, realBox_lo[n]);
      realBox.setHi(n, realBox_hi[n]);
    }

    geom->define(baseBox, &realBox, coord, is_periodic);
  }
}
#endif
