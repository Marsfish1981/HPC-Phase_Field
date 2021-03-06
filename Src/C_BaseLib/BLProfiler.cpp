#ifdef BL_PROFILING

#include <BLProfiler.H>
#include <REAL.H>
#include <Utility.H>
#include <ParallelDescriptor.H>
#include <Array.H>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <stack>
#include <algorithm>
#include <limits>
#include <stdlib.h>
#include <cmath>


bool BLProfiler::bWriteAll = true;
bool BLProfiler::bNoOutput = false;
bool BLProfiler::bWriteFabs = true;
bool BLProfiler::bFirstCommWriteH = true;  // header
bool BLProfiler::bFirstCommWriteD = true;  // data
bool BLProfiler::bInitialized = false;

int BLProfiler::currentStep = 0;
int BLProfiler::csFlushSize = 2000000;
int BLProfiler::traceFlushSize = 2000000;
int BLProfiler::nProfFiles  = 64;
int BLProfiler::finestLevel = -1;
int BLProfiler::maxLevel    = -1;

Real BLProfiler::pctTimeLimit = 5.0;
Real BLProfiler::calcRunTime  = 0.0;
Real BLProfiler::startTime    = 0.0;
Real BLProfiler::timerTime    = 0.0;

#ifndef BL_AMRPROF
Array<IntVect> BLProfiler::refRatio;
Array<Box> BLProfiler::probDomain;
#endif

std::stack<Real> BLProfiler::nestedTimeStack;
std::map<int, Real> BLProfiler::mStepMap;
std::map<std::string, BLProfiler::ProfStats> BLProfiler::mProfStats;
Array<BLProfiler::CommStats> BLProfiler::vCommStats;
std::map<std::string, BLProfiler *> BLProfiler::mFortProfs;
Array<std::string> BLProfiler::mFortProfsErrors;
const int mFortProfMaxErrors(32);
Array<BLProfiler *> BLProfiler::mFortProfsInt;
Array<std::string> BLProfiler::mFortProfsIntNames;
const int mFortProfsIntMaxFuncs(32);
std::map<std::string, BLProfiler::CommFuncType> BLProfiler::CommStats::cftNames;
std::set<BLProfiler::CommFuncType> BLProfiler::CommStats::cftExclude;
int BLProfiler::CommStats::barrierNumber(0);
int BLProfiler::CommStats::reductionNumber(0);
int BLProfiler::CommStats::tagWrapNumber(0);
int BLProfiler::CommStats::tagMin(0);
int BLProfiler::CommStats::tagMax(0);
int BLProfiler::CommStats::csVersion(1);
Array<std::pair<std::string,int> > BLProfiler::CommStats::barrierNames;
Array<std::pair<int,int> > BLProfiler::CommStats::nameTags;
Array<std::string> BLProfiler::CommStats::nameTagNames;
Array<int> BLProfiler::CommStats::reductions;
Array<int> BLProfiler::CommStats::tagWraps;

std::string BLProfiler::procName("NoProcName");
int BLProfiler::procNumber(-1);
bool BLProfiler::blProfDirCreated(false);
std::string BLProfiler::blProfDirName("bl_prof");
int BLProfiler::BLProfVersion(1);

std::map<std::string, int> BLProfiler::mFNameNumbers;
Array<BLProfiler::CallStats> BLProfiler::vCallTrace;

// Region support
std::map<std::string, int> BLProfiler::mRegionNameNumbers;
int BLProfiler::inNRegions(0);
Array<BLProfiler::RStartStop> BLProfiler::rStartStop;
const std::string BLProfiler::noRegionName("__NoRegion__");

bool BLProfiler::bFirstTraceWriteH(true);  // header
bool BLProfiler::bFirstTraceWriteD(true);  // data
int BLProfiler::CallStats::cstatsVersion(1);

Array<BLProfiler::CallStatsStack> BLProfiler::callIndexStack;
Array<BLProfiler::CallStatsPatch> BLProfiler::callIndexPatch;

#ifdef BL_TRACE_PROFILING
int BLProfiler::callStackDepth(-1);
int BLProfiler::prevCallStackDepth(0);
Real BLProfiler::CallStats::minCallTime(std::numeric_limits<Real>::max());
Real BLProfiler::CallStats::maxCallTime(std::numeric_limits<Real>::min());
#endif


BLProfiler::BLProfiler(const std::string &funcname)
    : bltstart(0.0), bltelapsed(0.0)
    , fname(funcname)
    , bRunning(false)
{
    start();
}


BLProfiler::BLProfiler(const std::string &funcname, bool bstart)
    : bltstart(0.0), bltelapsed(0.0)
    , fname(funcname)
    , bRunning(false)
{
    if(bstart) {
      start();
    }
}


BLProfiler::~BLProfiler() {
  if(bRunning) {
    stop();
  }
}


void BLProfiler::Initialize() {
  if(bInitialized) {
    return;
  }

  startTime = ParallelDescriptor::second();

  int resultLen(-1);
  char cProcName[MPI_MAX_PROCESSOR_NAME + 11];
#ifdef BL_USE_MPI
  MPI_Get_processor_name(cProcName, &resultLen);
#endif
  if(resultLen < 1) {
    procName   = "NoProcName";
    procNumber = ParallelDescriptor::MyProc();
  } else {
    procName = cProcName;
#ifdef BL_HOPPER
    procNumber = atoi(procName.substr(3, std::string::npos).c_str());
#else
#ifdef BL_SIM_HOPPER
    //procNumber = (100 * ParallelDescriptor::MyProc()) % 6527;
    procNumber = ParallelDescriptor::MyProc();
    std::stringstream pname;
    pname << "nid" << procNumber;
    procName = pname.str();
    std::cout << ParallelDescriptor::MyProc() << "::procNumber = " << procNumber
              << "  procName = " << procName << std::endl;
#else
    procNumber = ParallelDescriptor::MyProc();
#endif
#endif
  }
  //std::cout << myProc << ":::: " << procName << "  len =  " << resultLen << std::endl;

  Real t0, t1;
  int nTimerTimes(1000);
  for(int i(0); i < nTimerTimes; ++i) {  // ---- time the timer
    t0 = ParallelDescriptor::second();
    t1 = ParallelDescriptor::second();
    timerTime += t1 - t0;
  }
  timerTime /= static_cast<Real> (nTimerTimes);

#ifdef BL_COMM_PROFILING
  vCommStats.reserve(std::max(csFlushSize, 256000));
#endif

#ifdef BL_TRACE_PROFILING
  vCallTrace.reserve(std::max(traceFlushSize, 256000));
  // ---- make sure there is always at least one so we dont need to check in start()
  CallStats unusedCS(-1, -1, -1, -1.1, -1.2, -1.3);
  vCallTrace.push_back(unusedCS);
#endif
  BL_PROFILE_REGION_START(noRegionName);

  CommStats::cftExclude.insert(AllCFTypes);  // temporarily

  CommStats::cftNames["InvalidCFT"]     = InvalidCFT;
  CommStats::cftNames["AllReduceT"]     = AllReduceT; 
  CommStats::cftNames["AllReduceR"]     = AllReduceR; 
  CommStats::cftNames["AllReduceL"]     = AllReduceL; 
  CommStats::cftNames["AllReduceI"]     = AllReduceI; 
  CommStats::cftNames["AsendTsii"]      = AsendTsii;
  CommStats::cftNames["AsendTsiiM"]     = AsendTsiiM;
  CommStats::cftNames["AsendvTii"]      = AsendvTii;
  CommStats::cftNames["SendTsii"]       = SendTsii;
  CommStats::cftNames["SendvTii"]       = SendvTii;
  CommStats::cftNames["ArecvTsii"]      = ArecvTsii;
  CommStats::cftNames["ArecvTsiiM"]     = ArecvTsiiM;
  CommStats::cftNames["ArecvTii"]       = ArecvTii;
  CommStats::cftNames["ArecvvTii"]      = ArecvvTii;
  CommStats::cftNames["RecvTsii"]       = RecvTsii;
  CommStats::cftNames["RecvvTii"]       = RecvvTii;
  CommStats::cftNames["ReduceT"]        = ReduceT; 
  CommStats::cftNames["ReduceR"]        = ReduceR; 
  CommStats::cftNames["ReduceL"]        = ReduceL; 
  CommStats::cftNames["ReduceI"]        = ReduceI; 
  CommStats::cftNames["BCastTsi"]       = BCastTsi; 
  CommStats::cftNames["GatherTsT1Si"]   = GatherTsT1Si; 
  CommStats::cftNames["GatherTi"]       = GatherTi; 
  CommStats::cftNames["GatherRiRi"]     = GatherRiRi; 
  CommStats::cftNames["ScatterTsT1si"]  = ScatterTsT1si; 
  CommStats::cftNames["Barrier"]        = Barrier;
  CommStats::cftNames["Waitsome"]       = Waitsome;
  CommStats::cftNames["NameTag"]        = NameTag;
  CommStats::cftNames["AllCFTypes"]     = AllCFTypes;
  CommStats::cftNames["NoCFTypes"]      = NoCFTypes;
  CommStats::cftNames["IOStart"]        = IOStart;
  CommStats::cftNames["IOEnd"]          = IOEnd;
  CommStats::cftNames["TagWrap"]        = TagWrap;
  CommStats::cftNames["Allgather"]      = Allgather;
  CommStats::cftNames["Alltoall"]       = Alltoall;
  CommStats::cftNames["Alltoallv"]      = Alltoallv;
  CommStats::cftNames["Gatherv"]        = Gatherv;
  CommStats::cftNames["Get_count"]      = Get_count;
  CommStats::cftNames["Iprobe"]         = Iprobe;
  CommStats::cftNames["Test"]           = Test;
  CommStats::cftNames["Wait"]           = Wait;
  CommStats::cftNames["Waitall"]        = Waitall;

  // check for exclude file
  std::string exFile("CommFuncExclude.txt");
  Array<CommFuncType> vEx;

  Array<char> fileCharPtr;
  bool bExitOnError(false);  // in case the file does not exist
  //ParallelDescriptor::ReadAndBcastFile(exFile, fileCharPtr, bExitOnError, ParallelDescriptor::CommunicatorAll());
  ParallelDescriptor::ReadAndBcastFile(exFile, fileCharPtr, bExitOnError);

  CommStats::cftExclude.erase(AllCFTypes);

  if(fileCharPtr.size() > 0) {
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream cfex(fileCharPtrString, std::istringstream::in);

    while( ! cfex.eof()) {
        std::string cft;
        cfex >> cft;
        if( ! cfex.eof()) {
	  vEx.push_back(CommStats::StringToCFT(cft));
	}
    }
    for(int i(0); i < vEx.size(); ++i) {
      CommStats::cftExclude.insert(vEx[i]);
    }
  }

  // initialize fort int profilers
  mFortProfsInt.resize(mFortProfsIntMaxFuncs + 1);    // use 0 for undefined
  mFortProfsIntNames.resize(mFortProfsIntMaxFuncs + 1);  // use 0 for undefined
  for(int i(0); i < mFortProfsInt.size(); ++i) {
    std::ostringstream fname;
    fname << "FORTFUNC_" << i;
    mFortProfsIntNames[i] = fname.str();
#ifdef DEBUG
    mFortProfsInt[i] = 0;
#else
    mFortProfsInt[i] = new BLProfiler(mFortProfsIntNames[i], false);  // dont start
#endif
  }
  bInitialized = true;
}


void BLProfiler::ChangeFortIntName(const std::string &fname, int intname) {
#ifdef DEBUG
    mFortProfsIntNames[intname] = fname;
#else
    delete mFortProfsInt[intname];
    mFortProfsInt[intname] = new BLProfiler(fname, false);  // dont start
#endif
}


void BLProfiler::PStart() {
  bltelapsed = 0.0;
  start();
}


void BLProfiler::PStop() {
  if(bRunning) {
    stop();
  }
}


void BLProfiler::start() {
#ifdef _OPENMP
#pragma omp master
#endif
{
  bltelapsed = 0.0;
  bltstart = ParallelDescriptor::second();
  ++mProfStats[fname].nCalls;
  bRunning = true;
  nestedTimeStack.push(0.0);

#ifdef BL_TRACE_PROFILING
  int fnameNumber;
  std::map<std::string, int>::iterator it = BLProfiler::mFNameNumbers.find(fname);
  if(it == BLProfiler::mFNameNumbers.end()) {
    fnameNumber = BLProfiler::mFNameNumbers.size();
    BLProfiler::mFNameNumbers.insert(std::pair<std::string, int>(fname, fnameNumber));
  } else {
    fnameNumber = it->second;
  }
  ++callStackDepth;
  BL_ASSERT(vCallTrace.size() > 0);
  Real calltime(bltstart - startTime);
  vCallTrace.push_back(CallStats(callStackDepth, fnameNumber, 1, 0.0, 0.0, calltime));
  CallStats::minCallTime = std::min(CallStats::minCallTime, calltime);
  CallStats::maxCallTime = std::max(CallStats::maxCallTime, calltime);

  callIndexStack.push_back(CallStatsStack(vCallTrace.size() - 1));
  prevCallStackDepth = callStackDepth;

#endif
}
}

  
void BLProfiler::stop() {
#ifdef _OPENMP
#pragma omp master
#endif
{
  double tDiff(ParallelDescriptor::second() - bltstart);
  double nestedTime(0.0);
  bltelapsed += tDiff;
  bRunning = false;
  Real thisFuncTime(bltelapsed);
  if( ! nestedTimeStack.empty()) {
    nestedTime    = nestedTimeStack.top();
    thisFuncTime -= nestedTime;
    nestedTimeStack.pop();
  }
  if( ! nestedTimeStack.empty()) {
    nestedTimeStack.top() += bltelapsed;
  }
  mProfStats[fname].totalTime += thisFuncTime;

#ifdef BL_TRACE_PROFILING
  prevCallStackDepth = callStackDepth;
  --callStackDepth;
  BL_ASSERT(vCallTrace.size() > 0);
  if(vCallTrace.back().csFNameNumber == mFNameNumbers[fname]) {
    vCallTrace.back().totalTime = thisFuncTime + nestedTime;
    vCallTrace.back().stackTime = thisFuncTime;
  }
  if( ! callIndexStack.empty()) {
    CallStatsStack &cis(callIndexStack.back());
    if(cis.bFlushed) {
      callIndexPatch[cis.index].callStats.totalTime = thisFuncTime + nestedTime;
      callIndexPatch[cis.index].callStats.stackTime = thisFuncTime;
    } else {
      vCallTrace[cis.index].totalTime = thisFuncTime + nestedTime;
      vCallTrace[cis.index].stackTime = thisFuncTime;
    }
    callIndexStack.pop_back();
  }
#endif
}
}


void BLProfiler::InitParams(const Real ptl, const bool writeall, const bool writefabs) {
  pctTimeLimit = ptl;
  bWriteAll = writeall;
  bWriteFabs = writefabs;
}


#ifndef BL_AMRPROF
void BLProfiler::InitAMR(const int flev, const int mlev, const Array<IntVect> &rr,
                        const Array<Box> pd)
{
  finestLevel = flev;
  maxLevel    = mlev;
  refRatio.resize(rr.size());
  probDomain.resize(pd.size());
  for(int i(0); i < rr.size(); ++i) {
    refRatio[i] = rr[i];
  }
  for(int i(0); i < pd.size(); ++i) {
    probDomain[i] = pd[i];
  }
}
#endif


void BLProfiler::AddStep(const int snum) {
  currentStep = snum;
  mStepMap.insert(std::map<int, Real>::value_type(currentStep,
                                                  ParallelDescriptor::second()));
}


void BLProfiler::RegionStart(const std::string &rname) {
  Real rsTime(ParallelDescriptor::second() - startTime);

  if(rname != noRegionName) {
    ++inNRegions;
  }
  if(inNRegions == 1) {
    RegionStop(noRegionName);
  }

  int rnameNumber;
  std::map<std::string, int>::iterator it = BLProfiler::mRegionNameNumbers.find(rname);
  if(it == BLProfiler::mRegionNameNumbers.end()) {
    rnameNumber = BLProfiler::mRegionNameNumbers.size();
    BLProfiler::mRegionNameNumbers.insert(std::pair<std::string, int>(rname, rnameNumber));
  } else {
    rnameNumber = it->second;
  }
  rStartStop.push_back(RStartStop(true, rnameNumber, rsTime));
}


void BLProfiler::RegionStop(const std::string &rname) {
  Real rsTime(ParallelDescriptor::second() - startTime);

  int rnameNumber;
  std::map<std::string, int>::iterator it = BLProfiler::mRegionNameNumbers.find(rname);
  if(it == BLProfiler::mRegionNameNumbers.end()) {  // ---- error
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "-------- error in RegionStop:  region " << rname
                << " never started."  << std::endl;
    }
    rnameNumber = BLProfiler::mRegionNameNumbers.size();
    BLProfiler::mRegionNameNumbers.insert(std::pair<std::string, int>(rname, rnameNumber));
  } else {
    rnameNumber = it->second;
  }
  rStartStop.push_back(RStartStop(false, rnameNumber, rsTime));

  if(rname != noRegionName) {
    --inNRegions;
  }
  if(inNRegions == 0) {
    RegionStart(noRegionName);
  }
}


void BLProfiler::Finalize() {
  if( ! bInitialized) {
    return;
  }
  if(bNoOutput) {
    bInitialized = false;
    return;
  }

  // --------------------------------------- gather global stats
  Real finalizeStart(ParallelDescriptor::second());  // time the timer
  const int nProcs(ParallelDescriptor::NProcs());
  //const int myProc(ParallelDescriptor::MyProc());
  const int iopNum(ParallelDescriptor::IOProcessorNumber());

  // filter out profiler communications.
  CommStats::cftExclude.insert(AllCFTypes);


  // -------- make sure the set of profiled functions is the same on all processors
  Array<std::string> localStrings, syncedStrings;
  bool alreadySynced;

  for(std::map<std::string, ProfStats>::const_iterator it = mProfStats.begin();
      it != mProfStats.end(); ++it)
  {
    localStrings.push_back(it->first);
  }
  BoxLib::SyncStrings(localStrings, syncedStrings, alreadySynced);

  if( ! alreadySynced) {  // ---- add the new name
    for(int i(0); i < syncedStrings.size(); ++i) {
      std::map<std::string, ProfStats>::const_iterator it =
                                          mProfStats.find(syncedStrings[i]);
      if(it == mProfStats.end()) {
        ProfStats ps;
        mProfStats.insert(std::pair<std::string, ProfStats>(syncedStrings[i], ps));
      }
    }
  }

  // ---- add the following names if they have not been called already
  // ---- they will be called below to write the database and the names
  // ---- need to be in the database before it is written
  Array<std::string> addNames;
  addNames.push_back("ParallelDescriptor::Send(Tsii)i");
  addNames.push_back("ParallelDescriptor::Recv(Tsii)i");
  addNames.push_back("ParallelDescriptor::Gather(TsT1si)d");
  addNames.push_back("ParallelDescriptor::Gather(TsT1si)l");
  for(int iname(0); iname < addNames.size(); ++iname) {
    std::map<std::string, ProfStats>::iterator it = mProfStats.find(addNames[iname]);
    if(it == mProfStats.end()) {
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "BLProfiler::Finalize:  adding name:  " << addNames[iname] << std::endl;
      }
      ProfStats ps;
      mProfStats.insert(std::pair<std::string, ProfStats>(addNames[iname], ps));
    }
  }

  // ---------------------------------- now collect global data onto the ioproc
  int maxlen(0);
  Array<Real> gtimes(1);
  Array<long> ncalls(1);
  if(ParallelDescriptor::IOProcessor()) {
    gtimes.resize(nProcs);
    ncalls.resize(nProcs);
  }

  for(std::map<std::string, ProfStats>::const_iterator it = mProfStats.begin();
      it != mProfStats.end(); ++it)
  {
    std::string profName(it->first);
    int pnLen(profName.size());
    maxlen = std::max(maxlen, pnLen);

    ProfStats &pstats = mProfStats[profName];
    if(nProcs == 1) {
      gtimes[0] = pstats.totalTime;
      ncalls[0] = pstats.nCalls;
    } else {
      ParallelDescriptor::Gather(&pstats.totalTime, 1, gtimes.dataPtr(), 1, iopNum);
      ParallelDescriptor::Gather(&pstats.nCalls, 1, ncalls.dataPtr(), 1, iopNum);
    }
    Real tsum(0.0), tmin(gtimes[0]), tmax(gtimes[0]), tavg(0.0), variance(0.0);
    long ncsum(0);
    if(ParallelDescriptor::IOProcessor()) {
      for(int i(0); i < gtimes.size(); ++i) {
        tsum += gtimes[i];
        tmin = std::min(tmin, gtimes[i]);
        tmax = std::max(tmax, gtimes[i]);
      }
      tavg = tsum / static_cast<Real> (gtimes.size());
      for(int i(0); i < gtimes.size(); ++i) {
        variance += (gtimes[i] - tavg) * (gtimes[i] - tavg);
      }
      pstats.minTime = tmin;
      pstats.maxTime = tmax;
      pstats.avgTime = tavg;
      pstats.variance = variance / static_cast<Real> (gtimes.size());  // n - 1 for sample

      for(int i(0); i < ncalls.size(); ++i) {
        ncsum += ncalls[i];
      }
      // uncomment for reporting total calls summed over all procs
      //pstats.nCalls = ncsum;
    }
  }

  // --------------------------------------- print global stats to cout
  if(ParallelDescriptor::IOProcessor()) {
    bool bWriteAvg(true);
    if(nProcs == 1) {
      bWriteAvg = false;
    }
    BLProfilerUtils::WriteStats(std::cout, mProfStats, mFNameNumbers, vCallTrace, bWriteAvg);
  }


  // --------------------------------------- print all procs stats to a file
  if(bWriteAll) {
    // ----
    // ---- if we use an unordered_map for mProfStats, copy to a sorted container
    // ----
    ParallelDescriptor::Barrier();  // ---- wait for everyone (remove after adding filters)
    Array<long> nCallsOut(mProfStats.size(), 0);
    Array<Real> totalTimesOut(mProfStats.size(), 0.0);
    int count(0);
    for(std::map<std::string, ProfStats>::const_iterator phit = mProfStats.begin();
        phit != mProfStats.end(); ++phit)
    {
      nCallsOut[count] = phit->second.nCalls;
      totalTimesOut[count] = phit->second.totalTime;
      ++count;
    }


    // --------------------- start nfiles block
    std::string cdir(blProfDirName);
    if( ! blProfDirCreated) {
      BoxLib::UtilCreateCleanDirectory(cdir);
      blProfDirCreated = true;
    }

    const int   myProc    = ParallelDescriptor::MyProc();
    const int   nProcs    = ParallelDescriptor::NProcs();
    const int   nOutFiles = std::max(1, std::min(nProcs, nProfFiles));
    const int   nSets     = (nProcs + (nOutFiles - 1)) / nOutFiles;
    const int   mySet     = myProc/nOutFiles;
    std::string cFileName(cdir + '/' + cdir + "_D_");
    std::string FullName  = BoxLib::Concatenate(cFileName, myProc % nOutFiles, 4);

    long seekPos(0);

    for(int iSet = 0; iSet < nSets; ++iSet) {
      if(mySet == iSet) {
        {  // scope
          std::ofstream csFile;

          if(iSet == 0) {   // First set.
            csFile.open(FullName.c_str(),
                        std::ios::out|std::ios::trunc|std::ios::binary);
          } else {
            csFile.open(FullName.c_str(),
                        std::ios::out|std::ios::app|std::ios::binary);
            csFile.seekp(0, std::ios::end);   // set to eof
          }
          if( ! csFile.good()) {
            BoxLib::FileOpenFailed(FullName);
          }

          // ----------------------------- write to file here
	  seekPos = csFile.tellp();
	  if(nCallsOut.size() > 0) {
            csFile.write((char *) nCallsOut.dataPtr(),
	                 nCallsOut.size() * sizeof(long));
	  }
	  if(totalTimesOut.size() > 0) {
            csFile.write((char *) totalTimesOut.dataPtr(),
	                 totalTimesOut.size() * sizeof(Real));
	  }
          // ----------------------------- end write to file here

          csFile.flush();
          csFile.close();
        }  // end scope

        int iBuff     = 0;
        int wakeUpPID = (myProc + nOutFiles);
        int tag       = (myProc % nOutFiles);
        if(wakeUpPID < nProcs) {
          ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
        }
      }
      if(mySet == (iSet + 1)) {   // Next set waits.
        int iBuff;
        int waitForPID = (myProc - nOutFiles);
        int tag        = (myProc % nOutFiles);
        ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
      }
    }
    // --------------------- end nfiles block

    Array<long> seekPosOut(1);
    if(ParallelDescriptor::IOProcessor()) {
      seekPosOut.resize(nProcs, 0);
    }
    ParallelDescriptor::Gather(&seekPos, 1, seekPosOut.dataPtr(), 1, iopNum);

    if(ParallelDescriptor::IOProcessor()) {
      std::string phFilePrefix("bl_prof");
      std::string phFileName(cdir + '/' + phFilePrefix + "_H");
      std::ofstream phHeaderFile;
      phHeaderFile.open(phFileName.c_str(), std::ios::out | std::ios::trunc);
      phHeaderFile << "BLProfVersion " << BLProfVersion << '\n';
      phHeaderFile << "NProcs  " << nProcs << '\n';
      phHeaderFile << "NOutFiles  " << nOutFiles << '\n';
      for(std::map<std::string, ProfStats>::const_iterator phit = mProfStats.begin();
          phit != mProfStats.end(); ++phit)
      {
        phHeaderFile << "phFName " << '"' << phit->first << '"' << '\n';
      }

      std::string dFileName(cdir + "_D_");
      for(int p(0); p < nProcs; ++p) {
        std::string dFullName  = BoxLib::Concatenate(dFileName, p % nOutFiles, 4);
	phHeaderFile << "BLProfProc " << p << " datafile " << dFullName
	             << " seekpos " << seekPosOut[p] << '\n';
      }
      phHeaderFile << "calcEndTime " << std::setprecision(16)
                   << ParallelDescriptor::second() - startTime << '\n';
      phHeaderFile.close();
    }


    BL_PROFILE_REGION_STOP(noRegionName);

    WriteCallTrace();

    ParallelDescriptor::Barrier("BLProfiler::Finalize");
  }
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "BLProfiler::Finalize():  time:  "   // time the timer
              << ParallelDescriptor::second() - finalizeStart << std::endl;
  }

#ifdef BL_COMM_PROFILING
  WriteCommStats();
#endif

  WriteFortProfErrors();
#ifdef DEBUG
#else
  for(int i(0); i < mFortProfsInt.size(); ++i) {
    delete mFortProfsInt[i];
  }
#endif

  bInitialized = false;
}


namespace BLProfilerUtils {

void WriteHeader(std::ostream &ios, const int colWidth,
                 const Real maxlen, const bool bwriteavg)
{
  if(bwriteavg) {
    ios << std::setfill('-') << std::setw(maxlen+4 + 7 * (colWidth+2))
        << std::left << "Total times " << '\n';
    ios << std::right << std::setfill(' ');
    ios << std::setw(maxlen + 2) << "Function Name"
        << std::setw(colWidth + 2) << "NCalls"
        << std::setw(colWidth + 2) << "Min"
        << std::setw(colWidth + 2) << "Avg"
        << std::setw(colWidth + 2) << "Max"
        << std::setw(colWidth + 2) << "StdDev"
        << std::setw(colWidth + 2) << "CoeffVar"
        << std::setw(colWidth + 4) << "Percent %"
        << '\n';
  } else {
    ios << std::setfill('-') << std::setw(maxlen+4 + 3 * (colWidth+2))
        << std::left << "Total times " << '\n';
    ios << std::right << std::setfill(' ');
    ios << std::setw(maxlen + 2) << "Function Name"
        << std::setw(colWidth + 2) << "NCalls"
        << std::setw(colWidth + 2) << "Time"
        << std::setw(colWidth + 4) << "Percent %"
        << '\n';
  }
}


void WriteRow(std::ostream &ios, const std::string &fname,
              const BLProfiler::ProfStats &pstats, const Real percent,
	      const int colWidth, const Real maxlen,
	      const bool bwriteavg)
{
    int numPrec(4), pctPrec(2);
    Real stdDev(0.0), coeffVariation(0.0);
    if(pstats.variance > 0.0) {
      stdDev = std::sqrt(pstats.variance);
    }
    if(pstats.avgTime > 0.0) {
      coeffVariation = 100.0 * (stdDev / pstats.avgTime);  // ---- percent
    }

    if(bwriteavg) {
      ios << std::right;
      ios << std::setw(maxlen + 2) << fname << "  "
          << std::setw(colWidth) << pstats.nCalls << "  "
          << std::setprecision(numPrec) << std::fixed << std::setw(colWidth)
	  << pstats.minTime << "  "
          << std::setprecision(numPrec) << std::fixed << std::setw(colWidth)
	  << pstats.avgTime << "  "
          << std::setprecision(numPrec) << std::fixed << std::setw(colWidth)
	  << pstats.maxTime << "  "
          << std::setprecision(numPrec) << std::fixed << std::setw(colWidth)
	  << stdDev << "  "
          << std::setprecision(numPrec) << std::fixed << std::setw(colWidth)
	  << coeffVariation << "  "
          << std::setprecision(pctPrec) << std::fixed << std::setw(colWidth)
	  << percent << " %" << '\n';
    } else {
      ios << std::setw(maxlen + 2) << fname << "  "
          << std::setw(colWidth) << pstats.nCalls << "  "
          << std::setprecision(numPrec) << std::fixed << std::setw(colWidth)
	  << pstats.totalTime << "  "
          << std::setprecision(pctPrec) << std::fixed << std::setw(colWidth)
	  << percent << " %" << '\n';
    }
}


void WriteStats(std::ostream &ios,
                const std::map<std::string, BLProfiler::ProfStats> &mpStats,
		const std::map<std::string, int> &fnameNumbers,
		const Array<BLProfiler::CallStats> &callTraces,
		bool bwriteavg, bool bwriteinclusivetimes)
{
  const int myProc(ParallelDescriptor::MyProc());
  const int colWidth(10);
  const Real calcRunTime(BLProfiler::GetRunTime());

  std::map<Real, std::string, std::greater<Real> > mTimersTotalsSorted;

  Real totalTimers(0.0), percent(0.0);
  int maxlen(0);
  for(std::map<std::string, BLProfiler::ProfStats>::const_iterator it = mpStats.begin();
      it != mpStats.end(); ++it)
  {
    std::string profName(it->first);
    int pnLen(profName.size());
    maxlen = std::max(maxlen, pnLen);

    if(bwriteavg) {
      totalTimers += it->second.avgTime;
    } else {
      totalTimers += it->second.totalTime;
    }
  }
  Real pTimeTotal(totalTimers);
  if(calcRunTime > 0.0 && bwriteavg == false) {
    pTimeTotal = calcRunTime;
  }

  ios << '\n' << '\n';
  if( ! bwriteavg) {
    ios << std::setfill('*')
        << std::setw(maxlen + 2 + 3 * (colWidth + 2) - (colWidth+12)) << "";
    ios << std::setfill(' ');
    ios << "  Processor:  " << std::setw(colWidth) << myProc << '\n';
  }

  // -------- write timers sorted by name
  BLProfilerUtils::WriteHeader(ios, colWidth, maxlen, bwriteavg);
  for(std::map<std::string, BLProfiler::ProfStats>::const_iterator it = mpStats.begin();
      it != mpStats.end(); ++it)
  {
    if(pTimeTotal > 0.0) {
      if(bwriteavg) {
        percent = 100.0 * (it->second.avgTime / pTimeTotal);
      } else {
        percent = 100.0 * (it->second.totalTime / pTimeTotal);
      }
    } else {
      percent = 100.0;
    }
    std::string fname(it->first);
    const BLProfiler::ProfStats &pstats = it->second;
    BLProfilerUtils::WriteRow(ios, fname, pstats, percent, colWidth, maxlen, bwriteavg);
  }
  ios << '\n';
  ios << "Total Timers     = " << std::setw(colWidth) << totalTimers
      << " seconds." << '\n';
  if(calcRunTime > 0.0) {
    percent = 100.0 * totalTimers / calcRunTime;
    ios << "Calc Run Time    = " << std::setw(colWidth) << calcRunTime
        << " seconds." << '\n';
    ios << "Percent Coverage = " << std::setw(colWidth) << percent << " %" << '\n';
  }

  // -------- write timers sorted by percent
  ios << '\n' << '\n';
  BLProfilerUtils::WriteHeader(ios, colWidth, maxlen, bwriteavg);
  for(std::map<std::string, BLProfiler::ProfStats>::const_iterator it = mpStats.begin();
      it != mpStats.end(); ++it)
  {
    Real dsec;
    if(bwriteavg) {
      dsec = it->second.avgTime;
    } else {
      dsec = it->second.totalTime;
    }
    std::string sfir(it->first);
    mTimersTotalsSorted.insert(std::make_pair(dsec, sfir));
  }

  for(std::map<Real, std::string>::const_iterator it = mTimersTotalsSorted.begin();
      it != mTimersTotalsSorted.end(); ++it)
  {
    if(pTimeTotal > 0.0) {
      percent = 100.0 * (it->first / pTimeTotal);
    } else {
      percent = 100.0;
    }
    std::string fname(it->second);
    std::map<std::string, BLProfiler::ProfStats>::const_iterator mpsit = mpStats.find(fname);
    if(mpsit != mpStats.end()) {
      const BLProfiler::ProfStats &pstats = mpsit->second;
      BLProfilerUtils::WriteRow(ios, fname, pstats, percent, colWidth, maxlen, bwriteavg);
    } else {
      // error:  should not be able to get here if names are synced
    }
  }
  if(bwriteavg) {
    ios << std::setfill('=') << std::setw(maxlen+4 + 7 * (colWidth+2)) << ""
        << '\n';
  } else {
    ios << std::setfill('=') << std::setw(maxlen+4 + 3 * (colWidth+2)) << ""
        << '\n';
  }
  ios << std::setfill(' ');
  ios << std::endl;


#ifdef BL_TRACE_PROFILING
  // -------- write timers sorted by inclusive times
  Array<std::string> fNumberNames(fnameNumbers.size());
  for(std::map<std::string, int>::const_iterator it = fnameNumbers.begin();
      it != fnameNumbers.end(); ++it)
  {
    fNumberNames[it->second] = it->first;
  }

  // sort by total time
  Array<BLProfiler::RIpair> funcTotalTimes(fnameNumbers.size());
  for(int i(0); i < funcTotalTimes.size(); ++i) {
    funcTotalTimes[i].first  = 0.0;
    funcTotalTimes[i].second = i;
  }
  Array<int> callStack(64, -1);
  int maxCSD(0);
  std::set<int> recursiveFuncs;
  for(int i(0); i < callTraces.size(); ++i) {
    const BLProfiler::CallStats &cs = callTraces[i];
    if(cs.csFNameNumber < 0) {  // ---- an unused cs
      continue;
    }
    int depth(cs.callStackDepth);
    maxCSD = std::max(maxCSD, depth);
    if(depth >= callStack.size()) {
      callStack.resize(depth + 1);
    }
    callStack[depth] = cs.csFNameNumber;
    bool recursiveCall(false);
    for(int d(0); d <  depth; ++d) {
      if(cs.csFNameNumber == callStack[d]) {
	recursiveFuncs.insert(cs.csFNameNumber);
        recursiveCall = true;
      }
    }
    if( ! recursiveCall) {
      funcTotalTimes[cs.csFNameNumber].first += cs.totalTime;
    }
  }

  ios << " MaxCallStackDepth = " << maxCSD << '\n';
  for(std::set<int>::iterator rfi = recursiveFuncs.begin(); rfi != recursiveFuncs.end(); ++rfi) {
    ios << " Recursive function:  " << fNumberNames[*rfi] << '\n';
  }
  ios << '\n';

  if(bwriteinclusivetimes) {
    std::sort(funcTotalTimes.begin(), funcTotalTimes.end(), BLProfiler::fTTComp());

    int numPrec(4);
    ios << '\n' << '\n';
    ios << std::setfill('-') << std::setw(maxlen+4 + 1 * (colWidth+2))
        << std::left << "Inclusive times " << '\n';
    ios << std::right << std::setfill(' ');
    ios << std::setw(maxlen + 2) << "Function Name"
        << std::setw(colWidth + 4) << "Time s"
        << '\n';

    for(int i(0); i < funcTotalTimes.size(); ++i) {
      ios << std::setw(maxlen + 2) << fNumberNames[funcTotalTimes[i].second] << "  "
          << std::setprecision(numPrec) << std::fixed << std::setw(colWidth)
          << funcTotalTimes[i].first << " s"
          << '\n';
    }
    ios << std::setfill('=') << std::setw(maxlen+4 + 1 * (colWidth+2)) << ""
        << '\n';
    ios << std::setfill(' ');
    ios << std::endl;
  }

#endif
}


}  // end namespace BLProfilerUtils


void BLProfiler::WriteCallTrace(bool bFlushing) {   // ---- write call trace data

    if(bFlushing) {
      int nCT(vCallTrace.size());
      ParallelDescriptor::ReduceIntMax(nCT);
      if(nCT < traceFlushSize) {
        if(ParallelDescriptor::IOProcessor()) {
          std::cout << "Bypassing call trace flush, nCT < traceFlushSize:  " << nCT
                    << "  " << traceFlushSize << std::endl;
        }
        return;
      } else {
        if(ParallelDescriptor::IOProcessor()) {
          std::cout << "Flushing call traces:  nCT traceFlushSize = " << nCT
                    << "  " << traceFlushSize << std::endl;
        }
      }
    }

    std::string cdir(blProfDirName);
    const int   myProc    = ParallelDescriptor::MyProc();
    const int   nProcs    = ParallelDescriptor::NProcs();
    const int   nOutFiles = std::max(1, std::min(nProcs, nProfFiles));
    const int   nSets     = (nProcs + (nOutFiles - 1)) / nOutFiles;
    const int   mySet     = myProc/nOutFiles;
    std::string cFilePrefix("bl_call_stats");
    std::string cFileName(cdir + '/' + cFilePrefix + "_D_");
    std::string FullName  = BoxLib::Concatenate(cFileName, myProc % nOutFiles, 4);

    if( ! blProfDirCreated) {
      BoxLib::UtilCreateCleanDirectory(cdir);
      blProfDirCreated = true;
    }

    // -------- make sure the set of region names is the same on all processors
    Array<std::string> localStrings, syncedStrings;
    bool alreadySynced;
    for(std::map<std::string, int>::iterator it = mRegionNameNumbers.begin();
        it != mRegionNameNumbers.end(); ++it)
    {
      localStrings.push_back(it->first);
    }
    BoxLib::SyncStrings(localStrings, syncedStrings, alreadySynced);

    if( ! alreadySynced) {  // ---- need to remap names and numbers
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "**** Warning:  region names not synced:  unsupported."
                  << std::endl;
      }
      // unsupported for now
    }

    if(ParallelDescriptor::IOProcessor()) {
      std::string globalHeaderFileName(cdir + '/' + cFilePrefix + "_H");
      std::ofstream csGlobalHeaderFile;
      csGlobalHeaderFile.open(globalHeaderFileName.c_str(),
                              std::ios::out | std::ios::trunc);
      if( ! csGlobalHeaderFile.good()) {
        BoxLib::FileOpenFailed(globalHeaderFileName);
      }
      csGlobalHeaderFile << "CallStatsProfVersion  " << CallStats::cstatsVersion << '\n';
      csGlobalHeaderFile << "NProcs  " << nProcs << '\n';
      csGlobalHeaderFile << "NOutFiles  " << nOutFiles << '\n';

      for(std::map<std::string, int>::iterator it = mRegionNameNumbers.begin();
          it != mRegionNameNumbers.end(); ++it)
      {
        csGlobalHeaderFile << "RegionName " << '"' << it->first << '"'
	                   << ' ' << it->second << '\n';
      }
      for(int i(0); i < nOutFiles; ++i) {
        std::string headerName(cFilePrefix + "_H_");
        headerName = BoxLib::Concatenate(headerName, i, 4);
        csGlobalHeaderFile << "HeaderFile " << headerName << '\n';
      }
      csGlobalHeaderFile.flush();
      csGlobalHeaderFile.close();
    }

    std::string shortDFileName(cFilePrefix + "_D_");
    shortDFileName  = BoxLib::Concatenate(shortDFileName, myProc % nOutFiles, 4);
    std::string longDFileName(cdir + '/' + shortDFileName);

    std::string shortHeaderFileName(cFilePrefix + "_H_");
    shortHeaderFileName  = BoxLib::Concatenate(shortHeaderFileName, myProc % nOutFiles, 4);
    std::string longHeaderFileName(cdir + '/' + shortHeaderFileName);

    long baseSeekPos(-1);

    for(int iSet = 0; iSet < nSets; ++iSet) {
      if(mySet == iSet) {
        {  // scope
          std::ofstream csDFile, csHeaderFile;

          if(iSet == 0 && bFirstTraceWriteH) {   // First set.
	    bFirstTraceWriteH = false;
            csHeaderFile.open(longHeaderFileName.c_str(), std::ios::out | std::ios::trunc);
          } else {
            csHeaderFile.open(longHeaderFileName.c_str(), std::ios::out | std::ios::app);
            csHeaderFile.seekp(0, std::ios::end);   // set to eof
          }
          if( ! csHeaderFile.good()) {
            BoxLib::FileOpenFailed(longHeaderFileName);
          }

          if(iSet == 0 && bFirstTraceWriteD) {   // First set.
	    bFirstTraceWriteD = false;
            csDFile.open(longDFileName.c_str(),
	                 std::ios::out | std::ios::trunc | std::ios::binary);
          } else {
            csDFile.open(longDFileName.c_str(),
	                 std::ios::out | std::ios::app | std::ios::binary);
            csDFile.seekp(0, std::ios::end);   // set to eof
          }
          if( ! csDFile.good()) {
            BoxLib::FileOpenFailed(longDFileName);
          }

          // ----------------------------- write to file here
          csHeaderFile << "CallStatsProc " << myProc
		       << " nRSS " << rStartStop.size()
		       << " nTraceStats " << vCallTrace.size()
		       << "  datafile  " << shortDFileName
		       << "  seekpos  " << csDFile.tellp()
	               << '\n';
	  for(std::map<std::string, int>::iterator it = mFNameNumbers.begin();
	      it != mFNameNumbers.end(); ++it)
	  {
	    csHeaderFile << "fName " << '"' << it->first << '"'
	                 << ' ' << it->second << '\n';
	  }
#ifdef BL_TRACE_PROFILING
	  csHeaderFile << std::setprecision(16) << "timeMinMax  "
	               << CallStats::minCallTime << ' '
	               << CallStats::maxCallTime << '\n';
#endif

          csHeaderFile.flush();
          csHeaderFile.close();

	  if(rStartStop.size() > 0) {
	    csDFile.write((char *) rStartStop.dataPtr(),
	                  rStartStop.size() * sizeof(RStartStop));
	  }
	  if(vCallTrace.size() > 0) {
	    baseSeekPos = csDFile.tellp();
	    csDFile.write((char *) vCallTrace.dataPtr(),
	                  vCallTrace.size() * sizeof(CallStats));
	  }

	  csDFile.flush();
	  csDFile.close();
          // ----------------------------- end write to file here
        }  // end scope

        int iBuff     = 0;
        int wakeUpPID = (myProc + nOutFiles);
        int tag       = (myProc % nOutFiles);
        if(wakeUpPID < nProcs) {
          ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
        }
      }
      if(mySet == (iSet + 1)) {   // Next set waits.
        int iBuff;
        int waitForPID = (myProc - nOutFiles);
        int tag        = (myProc % nOutFiles);
        ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
      }
    }
    // --------------------- end nfiles block


    if(bFlushing) {  // ---- save stacked CallStats
      for(int ci(0); ci < callIndexStack.size(); ++ci) {
	CallStatsStack &csStack = callIndexStack[ci];
	if( ! csStack.bFlushed) {
	  if(baseSeekPos < 0) {
	    std::cout << "**** Error:  baseSeekPos = " << baseSeekPos << std::endl;
	    break;
	  }
	  long spos(baseSeekPos + csStack.index * sizeof(CallStats));
	  callIndexPatch.push_back(CallStatsPatch(spos, vCallTrace[csStack.index], longDFileName));
	  csStack.bFlushed = true;
	  csStack.index    = callIndexPatch.size() - 1;
	}
      }
    } else { // ---- patch the incomplete CallStats on disk
             // ---- probably should throttle these for large nprocs
      for(int ci(0); ci < callIndexPatch.size(); ++ci) {
	CallStatsPatch &csPatch = callIndexPatch[ci];

        std::fstream csDFile;
        CallStats csOnDisk;
        csDFile.open(csPatch.fileName.c_str(), std::ios::in | std::ios::out |
	                                       std::ios::binary);
        csDFile.seekg(csPatch.seekPos, std::ios::beg);
        csDFile.read((char *) &csOnDisk, sizeof(CallStats));
	bool bReportPatches(false);
	if(bReportPatches) {
	  std::cout << myProc << "::PATCH:  csOnDisk.st tt = " << csOnDisk.stackTime
	            << "  " << csOnDisk.totalTime << '\n'
		    << myProc << "::PATCH:  csPatch.st tt = " << csPatch.callStats.stackTime
		    << "  " << csPatch.callStats.totalTime << std::endl;
	}
	csOnDisk.totalTime = csPatch.callStats.totalTime;
	csOnDisk.stackTime = csPatch.callStats.stackTime;
        csDFile.seekp(csPatch.seekPos, std::ios::beg);
        csDFile.write((char *) &csOnDisk, sizeof(CallStats));
	csDFile.close();
      }
      callIndexPatch.clear();
      Array<CallStatsPatch>().swap(callIndexPatch);
    }

    // --------------------- delete the data
    rStartStop.clear();
    Array<RStartStop>().swap(rStartStop);
    vCallTrace.clear();
    Array<CallStats>().swap(vCallTrace);
    CallStats unusedCS(-1, -1, -1, -1.1, -1.2, -1.3);
    vCallTrace.push_back(unusedCS);
}



void BLProfiler::WriteCommStats(bool bFlushing) {

  Real wcsStart(ParallelDescriptor::second());
  bool bAllCFTypesExcluded(OnExcludeList(AllCFTypes));
  if( ! bAllCFTypesExcluded) {
    CommStats::cftExclude.insert(AllCFTypes);  // temporarily
  }

  if(bFlushing) {
    int nCS(vCommStats.size());
    ParallelDescriptor::ReduceIntMax(nCS);
    if(nCS < csFlushSize) {
      if( ! bAllCFTypesExcluded) {
        CommStats::cftExclude.erase(AllCFTypes);
      }
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "Bypassing comm stats flush, nCS < csFlushSize:  " << nCS
	          << "  " << csFlushSize << std::endl;
      }
      return;
    } else {
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "Flushing commstats:  nCSmax csFlushSize = " << nCS
	          << "  " << csFlushSize << std::endl;
      }
    }
  }

  std::string cdir(blProfDirName);
  std::string commprofPrefix("bl_comm_prof");
  if( ! blProfDirCreated) {
    BoxLib::UtilCreateCleanDirectory(cdir);
    blProfDirCreated = true;
  }

  bool bUseRelativeTimeStamp(true);
  if(bUseRelativeTimeStamp) {
    for(int ics(0); ics < BLProfiler::vCommStats.size(); ++ics) {
      CommStats &cs = BLProfiler::vCommStats[ics];
      cs.timeStamp -= startTime;
    }
  }


  // --------------------- start nfiles block
  const int   myProc    = ParallelDescriptor::MyProc();
  const int   nProcs    = ParallelDescriptor::NProcs();
  const int   nOutFiles = std::max(1, std::min(nProcs, nProfFiles));
  const int   nSets     = (nProcs + (nOutFiles - 1)) / nOutFiles;
  const int   mySet     = myProc/nOutFiles;

  std::string shortDFileName(commprofPrefix + "_D_");
  shortDFileName = BoxLib::Concatenate(shortDFileName, myProc % nOutFiles, 4);
  std::string longDFileName  = cdir + '/' + shortDFileName;

  std::string shortHeaderFileName(commprofPrefix + "_H_");
  shortHeaderFileName = BoxLib::Concatenate(shortHeaderFileName, myProc % nOutFiles, 4);
  std::string longHeaderFileName  = cdir + '/' + shortHeaderFileName;

  //double mpiWTick(MPI_Wtick());

  if(ParallelDescriptor::IOProcessor() && bFirstCommWriteH) {
    std::string globalHeaderFileName(cdir + '/' + commprofPrefix + "_H");
    std::ofstream csGlobalHeaderFile;
    csGlobalHeaderFile.open(globalHeaderFileName.c_str(), std::ios::out | std::ios::trunc);
    if( ! csGlobalHeaderFile.good()) {
      BoxLib::FileOpenFailed(globalHeaderFileName);
    }
    csGlobalHeaderFile << "CommProfVersion  " << CommStats::csVersion << '\n';
    csGlobalHeaderFile << "NProcs  " << nProcs << '\n';
    csGlobalHeaderFile << "CommStatsSize  " << sizeof(CommStats) << '\n';
    csGlobalHeaderFile << "NOutFiles  " << nOutFiles << '\n';
#ifndef BL_AMRPROF
    csGlobalHeaderFile << "FinestLevel  " << finestLevel << '\n';
    csGlobalHeaderFile << "MaxLevel  " << maxLevel << '\n';
    for(int i(0); i < refRatio.size(); ++i) {
      csGlobalHeaderFile << "RefRatio  " << i << "  " << refRatio[i] << '\n';
    }
    for(int i(0); i < probDomain.size(); ++i) {
      csGlobalHeaderFile << "ProbDomain  " << i << "  " << probDomain[i] << '\n';
    }
#endif
    for(int i(0); i < nOutFiles; ++i) {
      std::string headerName(commprofPrefix + "_H_");
      headerName = BoxLib::Concatenate(headerName, i, 4);
      csGlobalHeaderFile << "HeaderFile " << headerName << '\n';
    }

    csGlobalHeaderFile.flush();
    csGlobalHeaderFile.close();
  }


  for(int iSet = 0; iSet < nSets; ++iSet) {
    if(mySet == iSet) {
      {  // scope
        std::ofstream csDFile, csHeaderFile;

        if(iSet == 0 && bFirstCommWriteH) {   // First set.
	  bFirstCommWriteH = false;
          csHeaderFile.open(longHeaderFileName.c_str(), std::ios::out | std::ios::trunc);
        } else {
          csHeaderFile.open(longHeaderFileName.c_str(), std::ios::out | std::ios::app);
          csHeaderFile.seekp(0, std::ios::end);   // set to eof
        }
        if( ! csHeaderFile.good()) {
          BoxLib::FileOpenFailed(longHeaderFileName);
        }

        if(iSet == 0 && bFirstCommWriteD) {   // First set.
	  bFirstCommWriteD = false;
          csDFile.open(longDFileName.c_str(),
                      std::ios::out | std::ios::trunc | std::ios::binary);
        } else {
          csDFile.open(longDFileName.c_str(),
                      std::ios::out | std::ios::app | std::ios::binary);
          csDFile.seekp(0, std::ios::end);   // set to eof
        }
        if( ! csDFile.good()) {
          BoxLib::FileOpenFailed(longDFileName);
        }

        // ----------------------------- write to file here
        csHeaderFile << "CommProfProc  " << myProc
                     << "  nCommStats  " << vCommStats.size()
                     << "  datafile  " << shortDFileName
	             << "  seekpos  " << csDFile.tellp()
		     << "  " << procName << '\n';
        for(int ib(0); ib < CommStats::barrierNames.size(); ++ib) {
          int seekindex(CommStats::barrierNames[ib].second);
          CommStats &cs = vCommStats[seekindex];
          csHeaderFile << "bNum  " << cs.tag  // tag is used for barrier number
                       << ' ' << '"' << CommStats::barrierNames[ib].first << '"'
                       << ' ' << seekindex << '\n';
        }
        for(int ib(0); ib < CommStats::nameTags.size(); ++ib) {
          int seekindex(CommStats::nameTags[ib].second);
          csHeaderFile << "nTag  " << CommStats::nameTags[ib].first << ' '
                       << seekindex << '\n';
        }
        for(int ib(0); ib < CommStats::reductions.size(); ++ib) {
          int seekindex(CommStats::reductions[ib]);
          CommStats &cs = vCommStats[seekindex];
          csHeaderFile << "red  " << cs.tag  // tag is used for reduction number
	               << ' ' << seekindex << '\n';
        }
	if(vCommStats.size() > 0) {
	  csHeaderFile << std::setprecision(16) << std::setiosflags(std::ios::showpoint)
	               << "timeMinMax  " << vCommStats[0].timeStamp << ' '
	               << vCommStats[vCommStats.size()-1].timeStamp << '\n';
	  csHeaderFile << std::setprecision(16) << std::setiosflags(std::ios::showpoint)
	               << "timerTime  " << timerTime << '\n';
	} else {
	  csHeaderFile << "timeMinMax  0.0  0.0" << '\n';
	}
        for(int i(0); i < CommStats::nameTagNames.size(); ++i) {
          csHeaderFile << "nameTagNames  " << '"' << CommStats::nameTagNames[i]
                       << '"' << '\n';
        }
        csHeaderFile << "tagRange  " << CommStats::tagMin << ' '
	             << CommStats::tagMax << '\n';
        for(int i(0); i < CommStats::tagWraps.size(); ++i) {
          csHeaderFile << "tagWraps  " << CommStats::tagWraps[i] << '\n';
        }
	csHeaderFile.flush();
        csHeaderFile.close();

	if(vCommStats.size() > 0) {
	  csDFile.write((char *) vCommStats.dataPtr(), vCommStats.size() * sizeof(CommStats));
	}

        csDFile.flush();
        csDFile.close();
        // ----------------------------- end write to file here
      }  // end scope

      int iBuff     = 0;
      int wakeUpPID = (myProc + nOutFiles);
      int tag       = (myProc % nOutFiles);
      if(wakeUpPID < nProcs) {
        ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
      }
    }
    if(mySet == (iSet + 1)) {   // Next set waits.
      int iBuff;
      int waitForPID = (myProc - nOutFiles);
      int tag        = (myProc % nOutFiles);
      ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
    }
  }
  // --------------------- end nfiles block


  // --------------------- delete the data
  vCommStats.clear();
  Array<CommStats>().swap(vCommStats);
  CommStats::barrierNames.clear();
  Array<std::pair<std::string,int> >().swap(CommStats::barrierNames);
  CommStats::nameTags.clear();
  Array<std::pair<int,int> >().swap(CommStats::nameTags);
  CommStats::reductions.clear();
  Array<int>().swap(CommStats::reductions);
  if( ! bAllCFTypesExcluded) {
    CommStats::cftExclude.erase(AllCFTypes);
  }

  ParallelDescriptor::Barrier("BLProfiler::WriteCommStats::end");

  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "BLProfiler::WriteCommStats():  time:  "
              << ParallelDescriptor::second() - wcsStart << std::endl;
  }
}


void BLProfiler::WriteFortProfErrors() {
  // report any fortran errors.  should really check with all procs, just iop for now
  if(ParallelDescriptor::IOProcessor()) {
    if(BLProfiler::mFortProfs.size() > 0) {
      std::cout << "FFFFFFFF -------- FORTRAN PROFILING UNSTOPPED ERRORS" << std::endl;
      for(std::map<std::string, BLProfiler *>::iterator it = BLProfiler::mFortProfs.begin();
          it != BLProfiler::mFortProfs.end(); ++it)
      {
        std::cout << "FFFF function not stopped:  fname ptr = " << it->first
	          << "  ---->" << it->second << "<----" << std::endl;
      }
      std::cout << "FFFFFFFF -------- END FORTRAN PROFILING UNSTOPPED ERRORS" << std::endl;
    }
    if(BLProfiler::mFortProfsErrors.size() > 0) {
      std::cout << "FFFFFFFF FORTRAN PROFILING ERRORS" << std::endl;
      if(BLProfiler::mFortProfsErrors.size() >= mFortProfMaxErrors) {
        std::cout << "FFFFFFFF -------- MAX FORTRAN ERRORS EXCEEDED" << std::endl;
      }
      for(int i(0); i < BLProfiler::mFortProfsErrors.size(); ++i) {
        std::cout << "FFFF " << BLProfiler::mFortProfsErrors[i] << std::endl;
      }
      std::cout << "FFFFFFFF -------- END FORTRAN PROFILING ERRORS" << std::endl;
    }
  }
}


bool BLProfiler::OnExcludeList(CommFuncType cft) {
  // 
  // the idea for NoCFTypes is to allow local filtering/unfiltering
  // while preserving the users exclude list
  // possibly use a filter stack instead
  // might need caching if performance is a problem
  // what to do if both all and none are on the list?
  // 
  if(CommStats::cftExclude.empty()) {  // nothing on the exclude list
    return false;
  }
  std::set<CommFuncType>::iterator cfti;
  cfti = CommStats::cftExclude.find(NoCFTypes);
  if(cfti != CommStats::cftExclude.end()) {  // found, exclude nothing
    return false;
  }
  cfti = CommStats::cftExclude.find(cft);
  if(cfti != CommStats::cftExclude.end()) {  // found, exclude this cft
    return true;
  }
  cfti = CommStats::cftExclude.find(AllCFTypes);
  if(cfti != CommStats::cftExclude.end()) {  // found, exclude all types
    return true;
  }
  return false;
}


void BLProfiler::AddCommStat(const CommFuncType cft, const int size,
                           const int pid, const int tag)
{
  if(OnExcludeList(cft)) {
    return;
  }
  vCommStats.push_back(CommStats(cft, size, pid, tag, ParallelDescriptor::second()));
}


void BLProfiler::AddBarrier(const std::string &message, const bool beforecall) {
  const CommFuncType cft(BLProfiler::Barrier);
  if(OnExcludeList(cft)) {
    return;
  }
  if(beforecall) {
    int tag(CommStats::barrierNumber);
    vCommStats.push_back(CommStats(cft, 0, BeforeCall(), tag,
                                   ParallelDescriptor::second()));
    CommStats::barrierNames.push_back(std::make_pair(message, vCommStats.size() - 1));
    ++CommStats::barrierNumber;
  } else {
    int tag(CommStats::barrierNumber - 1);  // it was incremented before the call
    vCommStats.push_back(CommStats(cft, AfterCall(), AfterCall(), tag,
                                   ParallelDescriptor::second()));
  }
}


void BLProfiler::TagRange(const int mintag, const int maxtag) {
  CommStats::tagMin = mintag;
  CommStats::tagMax = maxtag;
}


void BLProfiler::AddTagWrap() {
  const CommFuncType cft(BLProfiler::TagWrap);
  if(OnExcludeList(cft)) {
    return;
  }
  ++CommStats::tagWrapNumber;
  int tag(CommStats::tagWrapNumber);
  int index(CommStats::nameTags.size());
  CommStats::tagWraps.push_back(index);
  vCommStats.push_back(CommStats(cft, index,  vCommStats.size(), tag,
                       ParallelDescriptor::second()));
}


void BLProfiler::AddAllReduce(const CommFuncType cft, const int size,
                            const bool beforecall)
{
  if(OnExcludeList(cft)) {
    return;
  }
  if(beforecall) {
    int tag(CommStats::reductionNumber);
    vCommStats.push_back(CommStats(cft, size, BeforeCall(), tag,
                                   ParallelDescriptor::second()));
    CommStats::reductions.push_back(vCommStats.size() - 1);
    ++CommStats::reductionNumber;
  } else {
    int tag(CommStats::reductionNumber - 1);
    vCommStats.push_back(CommStats(cft, size, AfterCall(), tag,
                                   ParallelDescriptor::second()));
  }
}


void BLProfiler::AddWaitsome(const CommFuncType cft, const Array<MPI_Request> &reqs,
                           const int completed, const Array<int> &indx,
			   const Array<MPI_Status> &status, const bool beforecall)
{
#ifdef BL_USE_MPI
  if(OnExcludeList(cft)) {
    return;
  }
  if(beforecall) {
    vCommStats.push_back(CommStats(cft, BeforeCall(), BeforeCall(), NoTag(),
                         ParallelDescriptor::second()));
  } else {
    for(int i(0); i < completed; ++i) {
      MPI_Status stat(status[i]);
      int c;
      BL_MPI_REQUIRE( MPI_Get_count(&stat, MPI_UNSIGNED_CHAR, &c) );
      vCommStats.push_back(CommStats(cft, c, stat.MPI_SOURCE, stat.MPI_TAG,
                           ParallelDescriptor::second()));
    }
  }
#endif
}


int BLProfiler::NameTagNameIndex(const std::string &name) {  // prob need to opt this
  for(int i(0); i < CommStats::nameTagNames.size(); ++i) {
    if(CommStats::nameTagNames[i] == name) {
      return i;
    }
  }
  CommStats::nameTagNames.push_back(name);
  return CommStats::nameTagNames.size() - 1;
}


void BLProfiler::AddNameTag(const std::string &name) {
  const CommFuncType cft(BLProfiler::NameTag);
  if(OnExcludeList(cft)) {
    return;
  }
  int tag(NameTagNameIndex(name));
  int index(CommStats::nameTags.size());
  vCommStats.push_back(CommStats(cft, index,  vCommStats.size(), tag,
                       ParallelDescriptor::second()));
  CommStats::nameTags.push_back(std::make_pair(tag, vCommStats.size() - 1));
}


BLProfiler::CommFuncType BLProfiler::CommStats::StringToCFT(const std::string &s) {
  return CommStats::cftNames[s];
}


void BLProfiler::CommStats::Filter(CommFuncType cft) {
  if( ! OnExcludeList(cft)) {
    CommStats::cftExclude.insert(cft);
  }
}


void BLProfiler::CommStats::UnFilter(CommFuncType cft) {
  if(OnExcludeList(cft)) {
    CommStats::cftExclude.erase(cft);
  }
}


namespace {
  const int EOS(-1);

  std::string Trim(const std::string &str) {
    int n;
    for(n = str.size(); --n >= 0; ) {
      if(str[n] != ' ' ) {
        break;
      }
    }
    std::string result;
    for(int i(0); i <= n; ++i) {
      result += str[i];
    }
    return result;
  }

  std::string Fint_2_string(const int *iarr, int nlen) {
    std::string res;
    for(int i(0); i < nlen && *iarr != EOS; ++i) {
      res += *iarr++;
    }
    return Trim(res);
  }
}

#include <BLFort.H>

BL_FORT_PROC_DECL(BL_PROFFORTFUNCSTART_CPP, bl_proffortfuncstart_cpp)
  (const int istr[], const int *NSTR)
{
  std::string fName(Fint_2_string(istr, *NSTR));
  std::map<std::string, BLProfiler *>::const_iterator it = BLProfiler::mFortProfs.find(fName);
  if(it == BLProfiler::mFortProfs.end()) {  // make a new profiler
    BLProfiler::mFortProfs.insert(std::pair<std::string, BLProfiler *>(fName, new BLProfiler(fName)));
  } else {  // error:  fname is already being profiled
    std::string estring("bl_proffortfuncstart error:  mFortProfs function already being profiled:  ");
    estring += fName;
    if(BLProfiler::mFortProfsErrors.size() < mFortProfMaxErrors) {
      BLProfiler::mFortProfsErrors.push_back(estring);
    }
  }
}

BL_FORT_PROC_DECL(BL_PROFFORTFUNCSTOP_CPP, bl_proffortfuncstop_cpp)
  (const int istr[], const int *NSTR)
{
  std::string fName(Fint_2_string(istr, *NSTR));
  std::map<std::string, BLProfiler *>::const_iterator it = BLProfiler::mFortProfs.find(fName);
  if(it == BLProfiler::mFortProfs.end()) {  // error:  fname not found
    std::string estring("bl_proffortfuncstop error:  mFortProfs function not started:  ");
    estring += fName;
    if(BLProfiler::mFortProfsErrors.size() < mFortProfMaxErrors) {
      BLProfiler::mFortProfsErrors.push_back(estring);
    }
  } else {  // delete the pointer and remove fname from map
    delete it->second;
    BLProfiler::mFortProfs.erase(fName);
  }
}


BL_FORT_PROC_DECL(BL_PROFFORTFUNCSTART_CPP_INT, bl_proffortfuncstart_cpp_int)
  (const int *i)
{
#ifdef DEBUG
  if(BLProfiler::mFortProfsInt[*i] == 0) {  // make a new profiler
    BLProfiler::mFortProfsInt[*i] =  new BLProfiler(BLProfiler::mFortProfsIntNames[*i]);
  } else {  // error:  fname is already being profiled
    std::string estring("bl_proffortfuncstart error:  mFortProfs function already being profiled:  ");
    estring += BLProfiler::mFortProfsIntNames[*i];
    if(BLProfiler::mFortProfsErrors.size() < mFortProfMaxErrors) {
      BLProfiler::mFortProfsErrors.push_back(estring);
    }
  }
#else
  BLProfiler::mFortProfsInt[*i]->PStart();
#endif
}


BL_FORT_PROC_DECL(BL_PROFFORTFUNCSTOP_CPP_INT, bl_proffortfuncstop_cpp_int)
  (const int *i)
{
#ifdef DEBUG
  if(BLProfiler::mFortProfsInt[*i] == 0) {  // error:  fname not found
    std::string estring("bl_proffortfuncstop error:  mFortProfs function not started:  ");
    estring += BLProfiler::mFortProfsIntNames[*i];
    if(BLProfiler::mFortProfsErrors.size() < mFortProfMaxErrors) {
      BLProfiler::mFortProfsErrors.push_back(estring);
    }
  } else {  // delete the pointer and remove fname from map
    delete BLProfiler::mFortProfsInt[*i];
    BLProfiler::mFortProfsInt[*i] = 0;
  }
#else
  BLProfiler::mFortProfsInt[*i]->PStop();
#endif
}


#else

#include <BLFort.H>

BL_FORT_PROC_DECL(BL_PROFFORTFUNCSTART_CPP,bl_proffortfuncstart_cpp)
  (
   const int istr[], const int *NSTR
   )
{
}

BL_FORT_PROC_DECL(BL_PROFFORTFUNCSTOP_CPP,bl_proffortfuncstop_cpp)
  (
   const int istr[], const int *NSTR
   )
{
}

BL_FORT_PROC_DECL(BL_PROFFORTFUNCSTART_CPP_INT,bl_proffortfuncstart_cpp_int)
  (
   int i
   )
{
}

BL_FORT_PROC_DECL(BL_PROFFORTFUNCSTOP_CPP_INT,bl_proffortfuncstop_cpp_int)
  (
   int i
   )
{
}


#endif





