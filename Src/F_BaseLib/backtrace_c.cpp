#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <limits>

#include <unistd.h>

#ifdef __linux__
#include <execinfo.h>
#endif

#include <signal.h>
#include <fenv.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef FORTRAN_BOXLIB
#include <BLBackTrace.H>
#endif

extern "C" {
    extern void abort_fortranboxlib ();
}

namespace
{
    static std::string fexename;
    static int myproc;

    void fflush_and_stderr (const char* str)
    {
	fflush(NULL);
	if (str)
	{
	    const char * const end = " !!!\n";
	    fwrite(str, strlen(str), 1, stderr);
	    fwrite(end, strlen(end), 1, stderr);
	}
    }

#ifdef __linux__
    void print_backtrace_info (FILE* f)
    {
	const int nbuf = 32;
	char **strings = NULL;
	void *buffer[nbuf];
	int nptrs = backtrace(buffer, nbuf);
	strings = backtrace_symbols(buffer, nptrs);
	if (strings != NULL) {
	    int have_addr2line = 0;
	    std::string cmd = "/usr/bin/addr2line";
	    if (FILE *fp = fopen(cmd.c_str(), "r")) {
		fclose(fp);
		have_addr2line = 1;
	    }
	    cmd += " -Cfie " + fexename; 
	    if (have_addr2line) {
		fprintf(f, "=== Please note that the line number reported by addr2line may not be accurate.\n");
		fprintf(f, "    If necessary, one can use 'readelf -wl my_exefile | grep my_line_address'\n");
		fprintf(f, "    to find out the offset for that line.\n");
	    }
	    for (int i = 0; i < nptrs; ++i) {
		std::string line = strings[i];
		line += "\n";
		if (have_addr2line) {
		    std::size_t found1 = line.rfind('[');
		    std::size_t found2 = line.rfind(']');
		    if (found1 != std::string::npos && found2 != std::string::npos) {
			std::string addr = line.substr(found1+1, found2-found1-1);
			std::string full_cmd = cmd + " " + addr;
			if (FILE * ps = popen(full_cmd.c_str(), "r")) {
			    char buff[512];
			    while (fgets(buff, sizeof(buff), ps)) {
				line += "    ";
				line += buff;
			    }
			    pclose(ps);
			}
		    }
		}
		fprintf(f, "%2d: %s\n", i, line.c_str());
	    }
	}
    }
#endif
}

extern "C"
{
#ifdef FORTRAN_BOXLIB
    void backtrace_handler (int s)
    {
	switch (s) {
	case SIGSEGV:
	    fflush_and_stderr("Segfault");
	    break;
	case SIGFPE:
	    fflush_and_stderr("Erroneous arithmetic operation");
	    break;
	case SIGINT:
	    fflush_and_stderr("SIGINT");
	    break;
	}

#ifdef __linux__
	std::string errfilename;
	{
	    std::ostringstream ss;
	    ss << "Backtrace." << myproc;
#ifdef _OPENMP
		ss << "." << omp_get_thread_num();
#endif
	    errfilename = ss.str();
	}
	if (FILE* p = fopen(errfilename.c_str(), "w")) {
	    print_backtrace_info(p);
	    fclose(p);
	}
	std::cerr << "See " << errfilename << " file for details" << std::endl;
	sleep(3);
#endif // __linux__
	
	abort_fortranboxlib();
    }
#else
    void backtrace_handler (int s)
    { 
	BLBackTrace::handler(s); 
    }
#endif

    void set_signal_handler (const char* ename, int rank)
    {
	if (ename[0] != '/') {
	    char temp[1024];
	    getcwd(temp,1024);
	    fexename = temp;
	    fexename += "/";
	}
	fexename += ename;
	
	signal(SIGSEGV, backtrace_handler); // catch seg falult
	signal(SIGINT,  backtrace_handler);
#if 0
#if (defined(BL_TESTING) || defined(DEBUG)) && defined(__linux__)
#if !defined(__PGI) || (__PGIC__ >= 16)
        // trap floating point exceptions
	feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
	signal(SIGFPE, backtrace_handler);
#endif
#endif
#endif

	myproc = rank;
    }

    int get_fpe_trap ()
    {
#if defined(__linux__) && (!defined(__PGI) || (__PGIC__ >= 16))
	static bool first = true;
	static int flags;
	if (first) {
	    first = false;
	    int all_enabled_exceptions = fegetexcept();
	    flags = (all_enabled_exceptions & FE_INVALID) 
		||  (all_enabled_exceptions & FE_DIVBYZERO) 
		||  (all_enabled_exceptions & FE_OVERFLOW);
		
	}
	return flags;
#else
	return 0;
#endif
    }

    double get_quiet_nan ()
    {
	return std::numeric_limits<double>::quiet_NaN();
    }
}
