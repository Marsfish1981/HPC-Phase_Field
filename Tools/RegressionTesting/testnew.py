#!/usr/bin/env python

"""
A simple regression test framework for a BoxLib-based code

There are several major sections to this source: the runtime parameter
routines, the test suite routines, and the report generation routines.
They are separated as such in this file.

This test framework understands source based out of the F_Src and
C_Src BoxLib frameworks.

"""

from __future__ import print_function

try: import ConfigParser as configparser
except ImportError:
    import configparser   # python 3

import datetime
import email
import argparse
import getpass
import json
import os
import shlex
import shutil
import smtplib
import socket
import string
import subprocess
import sys
import tarfile
import time

do_timings_plots = True

try: import numpy as np
except: do_timings_plots = False

try: import matplotlib
except: do_timings_plots = False
else:
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

try: import matplotlib.dates as dates
except: do_timings_plots = False

usage = """
The test suite and its tests are defined through an input file in an INI
configuration file format.

The "main" block specifies the global test suite parameters:

  [main]

  testTopDir     = < full path to test output directory >
  webTopDir      = < full path to test web output directory >

  sourceTree = < C_Src, F_Src, or BoxLib -- what type is it? >

  suiteName = < descriptive name (i.e. Castro) >

  reportActiveTestsOnly = <0: web shows every test ever run;
                           1: just current tests >

  goUpLink = <1: add "Go UP" link at top of the web page >

  FCOMP = < name of Fortran compiler >
  COMP  = < name of C/C++ compiler >

  add_to_f_make_command = < any additional defines to add to the make invocation for F_Src BoxLib >
  add_to_c_make_command = < any additional defines to add to the make invocation for C_Src BoxLib >

  purge_output = <0: leave all plotfiles in place;
                  1: delete plotfiles after compare >

  MAKE = < name of make >
  numMakeJobs = < number of make jobs >

  MPIcommand = < MPI run command, with holders for host, # of proc, command >

     This should look something like:

          mpiexec -host @host@ -n @nprocs@ @command@ >

  MPIhost = < host for MPI job -- depends on MPI implementation >

  sendEmailWhenFail = < 1: send email when any tests fail >

  emailTo = < list of email addresses separated by commas, such as,
              foo@example.com, bar@example.com >

  emailBody = < email body >


The source git repositories are defined in separate blocks.  There
will always be a "BoxLib" block, and usually a "source" block which is
the default directory used for compiling the tests.  Any extra repos
(including those where additional tests are to be build) are defined
in their own block starting with "extra-"

The general form is:

  [name]

  dir = < full path to git repo >

  branch = < desired branch in the git repo >

  build = < 1: this is a directory that tests will be compiled in >

  comp_string = < a string that is added to the make line >

      comp_string can refer to both the main source directory (as @source@)
      and its own directory (as @self@), for example:

      comp_string = CASTRO_DIR=@source@ WDMERGER_HOME=@self@


Each test is given its own block, with the general form:

  [Sod-x]

  buildDir = < relative path (from sourceDir) for this problem >

  inputFile = < input file name >
  probinFile = < probin file name >

  dim = < dimensionality: 1, 2, or 3 >

  aux?File = < name of additional file needed by the test >
  link?File = < name of additional file needed by the test >

      Here "?" is 1, 2, or 3, allowing for several files per test

  restartTest = < is this a restart test? 0 for no, 1 for yes >
  restartFileNum = < # of file to restart from (if restart test) >

  useMPI = <is this a parallel (MPI) job? 0 for no, 1 for yes) >
  numprocs = < # of processors to run on (if parallel job) >

  useOMP = <is this an OpenMP job? 0 for no, 1 for yes) >
  numthreads = < # of threads to us with OpenMP (if OpenMP job) >

  debug = < 0 for normal run, 1 if we want debugging options on >

  compileTest = < 0 for normal run, 1 if we just test compilation >

  selfTest = < 0 for normal run, 1 if test self-diagnoses if it succeeded >
  stSuccessString = < string to find in self-test output to determine success >

  doVis = < 0 for no visualization, 1 if we do visualization >
  visVar = < string of the variable to visualize >

  analysisRoutine = < name of the script to run on the output >

      The script is run as:

        analysisRoutine [options] plotfile

  analysisMainArgs = < commandline arguments to pass to the analysisRoutine --
                       these should refer to options from the [main] block >

  analysisOutputImage = < name on analysis result image to show on web page >

  compareFile = < explicit output file to do the comparison with -- this is
                  assumed to be prefixed with the test name when output by
                  the code at runtime, e.g. test_plt00100 >

  outputFile = < explicit output file to compare with -- exactly as it will
                 be written.  Not prefix of the test name will be done >

  diffDir = < directory/file to do a plain text diff on (recursive, if dir) >

  diffOpts = < options to use with the diff command for the diffDir comparison >


Getting started:

To set up a test suite, it is probably easiest to write the
testfile.ini as described above and then run the test routine with the
--make_benchmarks option to create the benchmark directory.
Subsequent runs can be done as usual, and will compare to the newly
created benchmarks.  If differences arise in the comparisons due to
(desired) code changes, the benchmarks can be updated using
--make_benchmarks to reflect the new ``correct'' solution.

"""


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# T E S T   C L A S S E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
class Test(object):

    def __init__ (self, name):

        self.name = name

        self.log = None

        self.buildDir = ""

        self.extra_build_dir = ""

        self.testSrcTree = ""

        self.inputFile = ""
        self.probinFile = ""
        self.auxFiles = []
        self.linkFiles = []

        self.dim = -1

        self.restartTest = 0
        self.restartFileNum = -1

        self.compileTest = 0

        self.selfTest = 0
        self.stSuccessString = ""

        self.debug = 0

        self.useMPI = 0
        self.numprocs = -1

        self.useOMP = 0
        self.numthreads = -1

        self.doVis = 0
        self.visVar = ""

        self.analysisRoutine = ""
        self.analysisMainArgs = ""
        self.analysisOutputImage = ""

        self.png_file = None

        self.outputFile = ""
        self.compareFile = ""

        self.compare_file_used = ""

        self.diffDir = ""
        self.diffOpts = ""

        self.addToCompileString = ""

        self.reClean = 0    # set automatically, not by users

        self.wall_time = 0   # set automatically, not by users

        self.nlevels = None  # set but running fboxinfo on the output

        self.comp_string = None  # set automatically
        self.run_command = None  # set automatically

        self.job_info_field1 = ""
        self.job_info_field2 = ""
        self.job_info_field3 = ""
        
        self.has_jobinfo = 0  # filled automatically
        
        self.backtrace = []   # filled automatically


    def __lt__(self, other):
        return self.value() < other.value()

    def value(self):
        return self.name

    def find_backtrace(self):
        """ find any backtrace files produced """        
        return [ft for ft in os.listdir(self.output_dir)
                if (os.path.isfile(ft) and ft.startswith("Backtrace."))]

    
    def get_last_plotfile(self, output_dir=None):
        """ Find the last plotfile written.  Note: we give an error if the
            last plotfile is 0.  If output_dir is specified, then we use
            that instead of the default
        """

        if output_dir is None:
            output_dir = self.output_dir   # not yet implemented

        plts = [d for d in os.listdir(output_dir) if \
                (os.path.isdir(d) and d.startswith("{}_plt".format(self.name))) or \
                (os.path.isfile(d) and d.startswith("{}_plt".format(self.name)) and d.endswith(".tgz"))]

        if len(plts) == 0:
            self.log.warn("WARNING: test did not produce any output")
            return ""

        plts.sort()
        last_plot = plts.pop()

        if last_plot.endswith("00000"):
            self.log.warn("WARNING: only plotfile 0 was output -- skipping comparison")
            return ""

        return last_plot


class Suite(object):

    def __init__ (self, args):

        self.args = args

        # this will hold all of the Repo() objects for the BoxLib, source,
        # and build directories
        self.repos = {}

        self.test_file_path = os.getcwd() + '/' + self.args.input_file[0]

        self.suiteName = "testDefault"
        self.sub_title = ""

        self.sourceTree = ""
        self.testTopDir = ""
        self.webTopDir = ""

        # set automatically
        self.source_dir = ""
        self.boxlib_dir = ""

        self.MPIcommand = ""
        self.MPIhost = ""

        self.FCOMP = "gfortran"
        self.COMP = "g++"

        self.add_to_f_make_command = ""
        self.add_to_c_make_command = ""

        self.summary_job_info_field1 = ""
        self.summary_job_info_field2 = ""
        self.summary_job_info_field3 = ""
        
        self.MAKE = "gmake"
        self.numMakeJobs = 1

        self.reportActiveTestsOnly = 0
        self.goUpLink = 0
        self.lenTestName = 0

        self.sendEmailWhenFail = 0
        self.emailFrom = ""
        self.emailTo = []
        self.emailSubject = ""
        self.emailBody = ""

        self.slack_post = 0
        self.slack_webhookfile = ""
        self.slack_webhook_url = None
        self.slack_channel = ""
        self.slack_username = ""

        self.globalAddToExecString = ""

        # delete all plot/checkfiles but the plotfile used for comparison upon
        # completion
        self.purge_output = 0

        self.log = None

    def check_test_dir(self, dir_name):
        """ given a string representing a directory, check if it points to
            a valid directory.  If so, return the directory name """

        dir_name = os.path.normpath(dir_name) + "/"

        if not os.path.isdir(dir_name):
            self.log.fail("ERROR: {} is not a valid directory".format(dir_name))

        return dir_name

    def get_tests_to_run(self, test_list_old):
        """ perform various tests based on the runtime options to determine
            which of the tests in the input file we run """

        # if we only want to run the tests that failed previously,
        # remove the others
        if self.args.redo_failed or not self.args.copy_benchmarks is None:
            last_run = self.get_last_run()
            failed = self.get_test_failures(last_run)

            test_list = [t for t in test_list_old if t.name in failed]
        else:
            test_list = test_list_old[:]

        # if we only want to run tests of a certain dimensionality, remove
        # the others
        if self.args.d in [1, 2, 3]:
            test_list = [t for t in test_list_old if t.dim == self.args.d]

        # if we are doing a single test, remove all other tests; if we
        # specified a list of tests, check each one; if we did both
        # --single_test and --tests, complain
        if not self.args.single_test == "" and not self.args.tests == "":
            self.log.fail("ERROR: specify tests either by --single_test or --tests, not both")

        if not self.args.single_test == "":
            tests_find = [self.args.single_test]
        elif not self.args.tests == "":
            tests_find = self.args.tests.split()
        else:
            tests_find = []

        if len(tests_find) > 0:
            new_test_list = []
            for test in tests_find:
                _tmp = [o for o in test_list if o.name == test]
                if len(_tmp) == 1:
                    new_test_list += _tmp
                else:
                    self.log.fail("ERROR: {} is not a valid test".format(test))

            test_list = new_test_list

        if len(test_list) == 0:
            self.log.fail("No valid tests defined")

        return test_list

    def get_bench_dir(self):
        bench_dir = self.testTopDir + self.suiteName + "-benchmarks/"
        if not os.path.isdir(bench_dir):
            if not self.args.make_benchmarks == None:
                os.mkdir(bench_dir)
            else:
                self.log.fail("ERROR: benchmark directory, %s, does not exist" % (bench_dir))
        return bench_dir

    def make_test_dirs(self):
        os.chdir(self.testTopDir)

        today_date = datetime.date.today()
        today = today_date.__str__()

        # figure out what the current output directory should be
        maxRuns = 100      # maximum number of tests in a given day

        test_dir = today + "/"

        # test output stored in a directory suiteName-tests/2007-XX-XX/
        # make sure that the suiteName-tests directory exists
        if not os.path.isdir(self.testTopDir + self.suiteName + "-tests/"):
            os.mkdir(self.testTopDir + self.suiteName + "-tests/")

        full_test_dir = self.testTopDir + self.suiteName + "-tests/" + test_dir

        if self.args.do_temp_run:
            test_dir = "TEMP_RUN/"
            full_test_dir = self.testTopDir + self.suiteName + "-tests/" + test_dir
            if os.path.isdir(full_test_dir):
                shutil.rmtree(full_test_dir)
        else:
            for i in range(1, maxRuns):
                if not os.path.isdir(full_test_dir): break
                test_dir = today + "-{:03d}/".format(i)
                full_test_dir = self.testTopDir + self.suiteName + "-tests/" + test_dir

        self.log.skip()
        self.log.bold("testing directory is: " + test_dir)
        os.mkdir(full_test_dir)

        # make the web directory -- this is where all the output and HTML will be
        # put, so it is easy to move the entire test website to a different disk
        full_web_dir = "%s/%s/"  % (self.webTopDir, test_dir)

        if self.args.do_temp_run:
            if os.path.isdir(full_web_dir):
                shutil.rmtree(full_web_dir)

        os.mkdir(full_web_dir)

        # copy the test file into the web output directory
        shutil.copy(self.test_file_path, full_web_dir)

        self.test_dir = test_dir
        self.full_test_dir = full_test_dir
        self.full_web_dir = full_web_dir

    def get_run_history(self, active_test_list):
        """ return the list of output directories run over the
            history of the suite and a separate list of the tests
            run (unique names) """

        valid_dirs = []
        all_tests = []

        # start by finding the list of valid test directories
        for file in os.listdir(self.webTopDir):

            # look for a directory of the form 20* (this will work up until 2099
            if file.startswith("20") and os.path.isdir(file):

                # look for the status file
                status_file = file + '/' + file + '.status'
                if os.path.isfile(status_file):
                    valid_dirs.append(file)

        valid_dirs.sort()
        valid_dirs.reverse()

        # now find all of the unique problems in the test directories
        for dir in valid_dirs:

            for file in os.listdir(self.webTopDir + dir):
                if file.endswith(".status") and not file.startswith("20"):
                    index = string.rfind(file, ".status")
                    test_name = file[0:index]

                    if all_tests.count(test_name) == 0:
                        if (not self.reportActiveTestsOnly) or (test_name in active_test_list):
                            all_tests.append(test_name)

        all_tests.sort()

        return valid_dirs, all_tests

    def make_timing_plots(self, active_test_list):
        """ plot the wallclock time history for all the valid tests """

        valid_dirs, all_tests = self.get_run_history(active_test_list)

        # store the timings in NumPy arrays in a dictionary
        timings = {}
        N = len(valid_dirs)
        for t in all_tests:
            timings[t] = np.zeros(N, dtype=np.float64)

        # now get the timings from the web output
        for n, d in enumerate(valid_dirs):
            for t in all_tests:
                ofile = "{}/{}/{}.html".format(self.webTopDir, d, t)
                try: f = open(ofile)
                except:
                    timings[t][n] = 0.0
                    continue

                found = False
                for line in f:
                    if "Execution time" in line:
                        found = True
                        # this is of the form: <li>Execution time: 412.930 s
                        timings[t][n] = float(line.split(":")[1].strip().split(" ")[0])
                        break

                    elif "(seconds)" in line:
                        found = True
                        # this is the older form -- split on "="
                        # form: <p><b>Execution Time</b> (seconds) = 399.414828
                        timings[t][n] = float(line.split("=")[1])
                        break

                f.close()
                if not found:
                    timings[t][n] = 0.0

        # make the plots
        for t in all_tests:
            _d = valid_dirs[:]
            _t = list(timings[t])

            days = []
            times = []
            for n, ttime in enumerate(_t):
                if not ttime == 0.0:
                    # sometimes the date is of the form YYYY-MM-DD-NNN, where NNN
                    # is the run -- remove that
                    date = _d[n]
                    if len(date) > 10:
                        date = date[:date.rfind("-")]

                    days.append(dates.datestr2num(date))
                    times.append(ttime)


            if len(times) == 0: continue

            plt.clf()
            plt.plot_date(days, times, "o", xdate=True)

            years = dates.YearLocator()   # every year
            months = dates.MonthLocator()
            yearsFmt = dates.DateFormatter('%Y')

            ax = plt.gca()
            ax.xaxis.set_major_locator(years)
            ax.xaxis.set_major_formatter(yearsFmt)
            ax.xaxis.set_minor_locator(months)

            plt.ylabel("time (seconds)")
            plt.title(t)

            if max(times)/min(times) > 10.0:
                ax.set_yscale("log")

            fig = plt.gcf()
            fig.autofmt_xdate()

            plt.savefig("{}/{}-timings.png".format(self.webTopDir, t))


    def get_last_run(self):
        """ return the name of the directory corresponding to the previous
            run of the test suite """

        outdir = self.testTopDir + self.suiteName + "-tests/"

        # this will work through 2099
        if os.path.isdir(outdir):
            dirs = [d for d in os.listdir(outdir) if (os.path.isdir(outdir + d) and
                                                      d.startswith("20"))]
            dirs.sort()

            return dirs[-1]
        else:
            return None

    def get_test_failures(self, test_dir):
        """ look at the test run in test_dir and return the list of tests that
            failed """

        cwd = os.getcwd()

        outdir = self.testTopDir + self.suiteName + "-tests/"

        os.chdir(outdir + test_dir)

        failed = []

        for test in os.listdir("."):
            if not os.path.isdir(test): continue

            # the status files are in the web dir
            status_file = "{}/{}/{}.status".format(self.webTopDir, test_dir, test)
            with open(status_file, "r") as sf:
                for line in sf:
                    if line.find("FAILED") >= 0:
                        failed.append(test)

        os.chdir(cwd)
        return failed

    def make_realclean(self, repo="source"):
        build_comp_string = ""
        if self.repos[repo].build == 1:
            if not self.repos[repo].comp_string is None:
                build_comp_string = self.repos[repo].comp_string

        extra_src_comp_string = ""
        if not self.extra_src_comp_string is None:
            extra_src_comp_string = self.extra_src_comp_string

        cmd = "{} BOXLIB_HOME={} {} {} realclean".format(
            self.MAKE, self.boxlib_dir,
            extra_src_comp_string, build_comp_string)

        run(cmd)

    def build_f(self, opts="", target="", outfile=None):
        comp_string = "{} -j{} BOXLIB_HOME={} COMP={} {} {} {}".format(
            self.MAKE, self.numMakeJobs, self.boxlib_dir,
            self.FCOMP, self.add_to_f_make_command, opts, target)
        self.log.log(comp_string)
        run(comp_string, outfile=outfile)
        return comp_string

    def run_test(self, test, base_command):
        test_env = None
        if test.useOMP:
            test_env = dict(os.environ, OMP_NUM_THREADS="{}".format(test.numthreads))

        if test.useMPI:
            test_run_command = self.MPIcommand
            test_run_command = test_run_command.replace("@host@", self.MPIhost)
            test_run_command = test_run_command.replace("@nprocs@", "{}".format(test.numprocs))
            test_run_command = test_run_command.replace("@command@", base_command)
        else:
            test_run_command = base_command

        self.log.log(test_run_command)
        sout, serr, ierr = run(test_run_command, stdin=True, outfile="{}.run.out".format(test.name), env=test_env)
        test.run_command = test_run_command

    def copy_backtrace(self, test):
        """
        if any backtrace files were output (because the run crashed), find them
        and copy them to the web directory
        """
        backtrace = test.find_backtrace()

        for btf in backtrace:
            ofile = "{}/{}.{}".format(self.full_web_dir, test.name, btf)
            shutil.copy(btf, ofile)
            test.backtrace.append("{}.{}".format(test.name, btf))


    def build_tools(self, test_list):

        self.compare_tool_dir = "{}/Tools/Postprocessing/F_Src/".format(
            os.path.normpath(self.boxlib_dir))

        os.chdir(self.compare_tool_dir)

        self.make_realclean(repo="BoxLib")

        tools = ["fcompare", "fboxinfo"]
        if any([t for t in test_list if t.dim == 2]): tools.append("fsnapshot2d")
        if any([t for t in test_list if t.dim == 3]): tools.append("fsnapshot3d")

        self.tools = {}

        self.log.skip()
        self.log.bold("building tools...")
        self.log.indent()

        for t in tools:
            self.log.log("building {}...".format(t))
            self.build_f(target="programs={}".format(t), opts="NDEBUG=t MPI= ")
            exe = get_recent_filename(self.compare_tool_dir, t, ".exe")
            self.tools[t] = "{}/{}".format(self.compare_tool_dir, exe)

        self.log.outdent()

    def slack_post_it(self, message):

        payload = {}

        # make sure there are no quotes in the strings
        payload["channel"] = self.slack_channel.replace('"', '')
        payload["username"] = self.slack_username.replace('"', '')
        payload["text"] = message

        s = json.dumps(payload)
        cmd = "curl -X POST --data-urlencode 'payload={}' {}".format(s, self.slack_webhook_url)
        print(cmd)
        run(cmd)


class Repo(object):
    """ a simple class to manage our git operations """
    def __init__(self, suite, directory, name,
                 branch_wanted=None, hash_wanted=None,
                 build=0, comp_string=None):

        self.suite = suite
        self.dir = directory
        self.name = name
        self.branch_wanted = branch_wanted
        self.hash_wanted = hash_wanted

        self.build = build   # does this repo contain build directories?
        self.comp_string = comp_string   # environment vars needed to build

        # for storage
        self.branch_orig = None
        self.hash_current = None

        self.update = True
        if hash_wanted:
            self.update = False

    def git_update(self):
        """ Do a git update of the repository.  If githash is not empty, then
            we will check out that version instead of git-pulling. """

        os.chdir(self.dir)

        # find out current branch so that we can go back later if we need.
        stdout0, stderr0, rc = run("git rev-parse --abbrev-ref HEAD")
        self.branch_orig = stdout0.rstrip('\n')

        if self.branch_orig != self.branch_wanted:
            self.suite.log.log("git checkout {} in {}".format(self.branch_wanted, self.dir))
            stdout, stderr, rc = run("git checkout {}".format(self.branch_wanted),
                                     stdin=True)
        else:
            self.branch_wanted = self.branch_orig

        if self.hash_wanted == "" or self.hash_wanted == None:
            self.suite.log.log("'git pull' in {}".format(self.dir))

            # we need to be tricky here to make sure that the stdin is
            # presented to the user to get the password.
            stdout, stderr, rc = run("git pull", stdin=True,
                                     outfile="git.{}.out".format(self.name))

        else:
            stdout, stderr, rc = run("git checkout {}".format(self.hash_wanted),
                                     outfile="git.{}.out".format(self.name))

        if not rc == 0:
            self.suite.log.fail("ERROR: git update was unsuccessful")

        shutil.copy("git.{}.out".format(self.name), self.suite.full_web_dir)

    def save_head(self):

        os.chdir(self.dir)

        self.suite.log.log("saving git HEAD for {}/".format(self.name))

        stdout, stderr, rc = run("git rev-parse HEAD",
                                 outfile="git.{}.HEAD".format(self.name) )

        self.hash_current = stdout
        shutil.copy("git.{}.HEAD".format(self.name), self.suite.full_web_dir)

    def make_changelog(self):
        """ generate a ChangeLog git repository, and copy it to the
            web directory"""

        os.chdir(self.dir)

        self.suite.log.log("generating ChangeLog for {}/".format(self.name))

        run("git log --name-only", outfile="ChangeLog.{}".format(self.name), outfile_mode="w")
        shutil.copy("ChangeLog.{}".format(self.name), self.suite.full_web_dir)

    def git_back(self):
        """ switch the repo back to its original branch """

        os.chdir(self.dir)
        self.suite.log.log("git checkout {} in {}".format(self.branch_orig, self.dir))

        stdout, stderr, rc = run("git checkout {}".format(self.branch_orig),
                                 stdin=True,
                                 outfile="git.{}.out".format(self.name))

        if not rc == 0:
            self.suite.log.fail("ERROR: git checkout was unsuccessful")


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# R U N T I M E   P A R A M E T E R   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

def convert_type(string):
    """ return an integer, float, or string from the input string """
    if string is None:
        return None
        
    try: int(string)
    except: pass
    else: return int(string)

    try: float(string)
    except: pass
    else: return float(string)

    return string.strip()

def safe_get(cp, sec, opt, default=None):
    try: v = cp.get(sec, opt)
    except: v = default
    return v

def load_params(args):
    """
    reads the parameter file and creates as list of test objects as well as
    the suite object
    """

    test_list = []

    cp = configparser.ConfigParser()    # note, need strict=False for Python3
    cp.optionxform = str

    log = Log()

    log.bold("loading " + args.input_file[0])

    try: cp.read(args.input_file[0])
    except:
        log.fail("ERROR: unable to read parameter file {}".format(file))

    # "main" is a special section containing the global suite parameters.
    mysuite = Suite(args)

    mysuite.log = log

    valid_options = list(mysuite.__dict__.keys())

    for opt in cp.options("main"):

        # get the value of the current option
        value = convert_type(cp.get("main", opt))

        if opt in valid_options:

            if opt == "sourceTree":
                if not value in ["C_Src", "F_Src", "BoxLib"]:
                    mysuite.log.fail("ERROR: invalid sourceTree")
                else:
                    mysuite.sourceTree = value

            elif opt == "testTopDir": mysuite.testTopDir = mysuite.check_test_dir(value)
            elif opt == "webTopDir": mysuite.webTopDir = os.path.normpath(value) + "/"

            elif opt == "emailTo": mysuite.emailTo = value.split(",")

            else:
                # generic setting of the object attribute
                setattr(mysuite, opt, value)

        else:
            mysuite.log.warn("WARNING: suite parameter %s not valid" % (opt))


    # BoxLib -- this will always be defined
    rdir = mysuite.check_test_dir(safe_get(cp, "BoxLib", "dir"))

    branch = convert_type(safe_get(cp, "BoxLib", "branch"))
    rhash = convert_type(safe_get(cp, "BoxLib", "hash"))

    mysuite.repos["BoxLib"] = Repo(mysuite, rdir, "BoxLib",
                                   branch_wanted=branch, hash_wanted=rhash)


    # now all the other build and source directories
    other_srcs = [s for s in cp.sections() if s.startswith("extra-")]
    if not mysuite.sourceTree == "BoxLib": other_srcs.append("source")

    for s in other_srcs:
        if s.startswith("extra-"):
            k = s.split("-")[1]
        else:
            k = "source"

        rdir = mysuite.check_test_dir(safe_get(cp, s, "dir"))
        branch = convert_type(safe_get(cp, s, "branch"))
        rhash = convert_type(safe_get(cp, s, "hash"))

        build = convert_type(safe_get(cp, s, "build", default=0))
        if s == "source": build = 1

        comp_string = safe_get(cp, s, "comp_string")

        name = os.path.basename(os.path.normpath(rdir))

        mysuite.repos[k] = Repo(mysuite, rdir, name,
                                branch_wanted=branch, hash_wanted=rhash,
                                build=build, comp_string=comp_string)


    # BoxLib-only tests don't have a sourceDir
    mysuite.boxlib_dir = mysuite.repos["BoxLib"].dir

    if mysuite.sourceTree == "BoxLib":
        mysuite.source_dir = mysuite.repos["BoxLib"].dir
    else:
        mysuite.source_dir = mysuite.repos["source"].dir


    # now flesh out the compile strings -- they may refer to either themselves
    # or the source dir
    for r in mysuite.repos.keys():
        s = mysuite.repos[r].comp_string
        if not s is None:
            mysuite.repos[r].comp_string = \
                s.replace("@self@", mysuite.repos[r].dir).replace("@source@", mysuite.repos["source"].dir)


    # the suite needs to know both ext_src_comp_string
    mysuite.extra_src_comp_string = ""
    for r in mysuite.repos.keys():
        if not mysuite.repos[r].build == 1:
            if not mysuite.repos[r].comp_string is None:
                mysuite.extra_src_comp_string += "{} ".format(mysuite.repos[r].comp_string)

    # checks
    if mysuite.sendEmailWhenFail and not args.send_no_email:
        if mysuite.emailTo == [] or mysuite.emailBody == "":
            mysuite.log.fail("ERROR: when sendEmailWhenFail = 1, you must specify emailTo and emailBody\n")

        if mysuite.emailFrom == "":
            mysuite.emailFrom = '@'.join((getpass.getuser(), socket.getfqdn()))

        if mysuite.emailSubject == "":
            mysuite.emailSubject = mysuite.suiteName+" Regression Test Failed"

    if mysuite.slack_post:
        if not os.path.isfile(mysuite.slack_webhookfile):
            mysuite.log.warn("ERROR: slack_webhookfile invalid")
            mysuite.slack_post = 0
        else:
            print(mysuite.slack_webhookfile)
            try: f = open(mysuite.slack_webhookfile)
            except:
                mysuite.log.warn("ERROR: unable to open webhook file")
                mysuite.slack_post = 0
            else:
                mysuite.slack_webhook_url = str(f.readline())
                f.close()

    if (mysuite.sourceTree == "" or mysuite.boxlib_dir == "" or
        mysuite.source_dir == "" or mysuite.testTopDir == ""):
        mysuite.log.fail("ERROR: required suite-wide directory not specified\n" + \
                         "(sourceTree, boxLibDir, sourceDir, testTopDir)")

    # Make sure the web dir is valid (or use the default is none specified)
    if mysuite.webTopDir == "":
        mysuite.webTopDir = "{}/{}-web/".format(mysuite.testTopDir, mysuite.suiteName)

    if not os.path.isdir(mysuite.webTopDir):
        try: os.mkdir(mysuite.webTopDir)
        except:
            mysuite.log.fail("ERROR: unable to create the web directory: {}\n".format(
                mysuite.webTopDir))

    # all other sections are tests
    mysuite.log.skip()
    mysuite.log.bold("finding tests and checking parameters...")

    for sec in cp.sections():

        if sec in ["main", "BoxLib", "source"] or sec.startswith("extra-"): continue

        # maximum test name length -- used for HTML formatting
        mysuite.lenTestName = max(mysuite.lenTestName, len(sec))

        # create the test object for this test
        mytest = Test(sec)
        mytest.log = log
        invalid = 0

        # set the test object data by looking at all the options in
        # the current section of the parameter file
        valid_options = list(mytest.__dict__.keys())
        valid_options += ["aux1File", "aux2File", "aux3File"]
        valid_options += ["link1File", "link2File", "link3File"]

        for opt in cp.options(sec):

            # get the value of the current option
            value = convert_type(cp.get(sec, opt))

            if opt in valid_options:

                if opt in ["aux1File", "aux2File", "aux3File"]:
                    mytest.auxFiles.append(value)

                elif opt in ["link1File", "link2File", "link3File"]:
                    mytest.linkFiles.append(value)

                else:
                    # generic setting of the object attribute
                    setattr(mytest, opt, value)

            else:
                mysuite.log.warn("WARNING: unrecognized parameter {} for test {}".format(opt, sec))


        # make sure that the build directory actually exists
        if not mytest.extra_build_dir == "":
            bDir = mysuite.repos[mytest.extra_build_dir].dir + mytest.buildDir
        else:
            bDir = mysuite.source_dir + mytest.buildDir

        if not os.path.isdir(bDir):
            mysuite.log.warn("WARNING: invalid build directory: {}".format(bDir))
            invalid = 1


        # make sure all the require parameters are present
        if mytest.compileTest:
            if mytest.buildDir == "":
                mysuite.log.warn("WARNING: mandatory parameters for test {} not set".format(sec))
                invalid = 1

        else:
            if (mytest.buildDir == "" or mytest.inputFile == "" or
                (mysuite.sourceTree == "C_Src" and mytest.probinFile == "") or
                mytest.dim == -1):
                mysuite.log.warn("WARNING: mandatory parameters for test {} not set".format(sec))
                mysuite.log.warn("         buildDir = {}".format(mytest.buildDir))
                mysuite.log.warn("         inputFile = {}".format(mytest.inputFile))
                if mysuite.sourceTree == "C_Src":
                    mysuite.log.warn("         probinFile = {}".format(mytest.probinFile))
                mysuite.log.warn("            dim = {}".format(mytest.dim))

                invalid = 1

        # check the optional parameters
        if mytest.restartTest and mytest.restartFileNum == -1:
            mysuite.log.warn("WARNING: restart-test {} needs a restartFileNum".format(sec))
            invalid = 1

        if mytest.selfTest and mytest.stSuccessString == "":
            mysuite.log.warn("WARNING: self-test {} needs a stSuccessString".format(sec))
            invalid = 1

        if mytest.useMPI and mytest.numprocs == -1:
            mysuite.log.warn("WARNING: MPI parallel test {} needs numprocs".format(sec))
            invalid = 1

        if mytest.useOMP and mytest.numthreads == -1:
            mysuite.log.warn("WARNING: OpenMP parallel test {} needs numthreads".format(sec))
            invalid = 1

        if mytest.doVis and mytest.visVar == "":
            mysuite.log.warn("WARNING: test {} has visualization, needs visVar".format(sec))
            invalid = 1

        if mysuite.sourceTree == "BoxLib" and mytest.testSrcTree == "":
            mysuite.log.warn("WARNING: test {} is a BoxLib test but testSrcTree not set".format(sec))
            invalid = 1


        # add the current test object to the master list
        if not invalid:
            test_list.append(mytest)
        else:
            mysuite.log.warn("WARNING: test {} will be skipped".format(sec))


    # if any runs are parallel, make sure that the MPIcommand is defined
    anyMPI = any([t.useMPI for t in test_list])

    if anyMPI and mysuite.MPIcommand == "":
        mysuite.log.fail("ERROR: some tests are MPI parallel, but MPIcommand not defined")

    test_list.sort()

    return mysuite, test_list


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# O U T P U T   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
class Log(object):
    def __init__(self, output_file=None):

        # http://stackoverflow.com/questions/287871/print-in-terminal-with-colors-using-python
        # which in-turn cites the blender build scripts
        self.warn_color = '\033[33m'
        self.success_color = '\033[32m'
        self.fail_color = '\033[31m'
        self.bold_color = '\033[1m'
        self.end_color = '\033[0m'

        self.current_indent = 0
        self.indent_str = ""

        if not output_file is None:
            self.of = output_file
        else:
            self.of = None

    def indent(self):
        self.current_indent += 1
        self.indent_str = self.current_indent*"   "

    def outdent(self):
        self.current_indent -= 1
        self.current_indent = max(0, self.current_indent)
        self.indent_str = self.current_indent*"   "

    def fail(self, string):
        nstr = self.fail_color + string + self.end_color
        print("{}{}".format(self.indent_str, nstr))
        self.close_log()
        sys.exit()

    def testfail(self, string):
        nstr = self.fail_color + string + self.end_color
        print("{}{}".format(self.indent_str, nstr))

    def warn(self, string):
        nstr = self.warn_color + string + self.end_color
        print("{}{}".format(self.indent_str, nstr))

    def success(self, string):
        nstr = self.success_color + string + self.end_color
        print("{}{}".format(self.indent_str, nstr))

    def log(self, string):
        print("{}{}".format(self.indent_str, string))

    def skip(self):
        print("")

    def bold(self, string):
        nstr = self.bold_color + string + self.end_color
        print("{}{}".format(self.indent_str, nstr))

    def close_log(self):
        if not self.of is None: self.of.close()



#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# S Y S T E M   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
def get_args(arg_string=None):
    """ parse the commandline arguments.  If arg_string is present, we
        parse from there, otherwise we use the default (sys.argv) """

    parser = argparse.ArgumentParser(description=usage, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-d", type=int, default=-1,
                        help="restrict tests to a particular dimensionality")
    parser.add_argument("--make_benchmarks", type=str, default=None, metavar="comment",
                        help="make new benchmarks? (must provide a comment)")
    parser.add_argument("--copy_benchmarks", type=str, default=None, metavar="comment",
                        help="copy the last plotfiles from the failed tests of the most recent run as the new benchmarks.  No git pull is done and no new runs are performed (must provide a comment)")
    parser.add_argument("--no_update", type=str, default="None", metavar="name",
                        help="which codes to exclude from the git update? (None, All, or a comma-separated list of codes)")
    parser.add_argument("--single_test", type=str, default="", metavar="test-name",
                        help="name of a single test to run")
    parser.add_argument("--tests", type=str, default="", metavar="'test1 test2 test3'",
                        help="a space-separated list of tests to run")
    parser.add_argument("--do_temp_run", action="store_true",
                        help="is this a temporary run? (output not stored or logged)")
    parser.add_argument("--send_no_email", action="store_true",
                        help="do not send emails when tests fail")
    parser.add_argument("--with_valgrind", action="store_true",
                        help="run with valgrind")
    parser.add_argument("--valgrind_options", type=str, default="--leak-check=yes --log-file=vallog.%p",
                        help="valgrind options")
    parser.add_argument("--boxLibGitHash", type=str, default=None, metavar="hash",
                        help="git hash of a version of BoxLib.  If provided, this version will be used to run tests.")
    parser.add_argument("--sourceGitHash", type=str, default=None, metavar="hash",
                        help="git hash of a version of the source code.  For BoxLib tests, this will be ignored.")
    parser.add_argument("--extSrcGitHash", type=str, default=None, metavar="hash",
                        help="git hash of a version of the source code.  For BoxLib tests, this will be ignored.")
    parser.add_argument("--note", type=str, default="",
                        help="a note on the resulting test webpages")
    parser.add_argument("--complete_report_from_crash", type=str, default="", metavar="testdir",
                        help="complete the generation of the report from a crashed test suite run named testdir")
    parser.add_argument("--redo_failed", action="store_true",
                        help="only run the tests that failed last time")
    parser.add_argument("input_file", metavar="input-file", type=str, nargs=1,
                        help="the input file (INI format) containing the suite and test parameters")

    if not arg_string is None:
        args = parser.parse_args(arg_string)
    else:
        args = parser.parse_args()

    return args


def run(string, stdin=False, outfile=None, store_command=False, env=None, outfile_mode="a", log=None):

    # shlex.split will preserve inner quotes
    prog = shlex.split(string)
    if stdin:
        p0 = subprocess.Popen(prog, stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, env=env)
    else:
        p0 = subprocess.Popen(prog, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, env=env)

    stdout0, stderr0 = p0.communicate()
    if stdin: p0.stdin.close()
    rc = p0.returncode
    p0.stdout.close()

    if not outfile == None:
        try: cf = open(outfile, outfile_mode)
        except IOError:
            log.fail("  ERROR: unable to open file for writing")
        else:
            if store_command:
                cf.write(string)
            for line in stdout0:
                cf.write(line)
            cf.close()

    return stdout0, stderr0, rc


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# T E S T   S U I T E   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

def find_build_dirs(tests):
    """ given the list of test objects, find the set of UNIQUE build
        directories.  Note if we have the useExtraBuildDir flag set """

    build_dirs = []
    reClean = []

    for obj in tests:

        # keep track of the build directory and which source tree it is
        # in (e.g. the extra build dir)

        # first find the list of unique build directories
        dir_pair = (obj.buildDir, obj.extra_build_dir)
        if build_dirs.count(dir_pair) == 0:
            build_dirs.append(dir_pair)


        # re-make all problems that specify an extra compile argument,
        # just to make sure that any unique build commands are seen.
        if not obj.addToCompileString == "":
            reClean.append(dir_pair)

    for bdir, _ in reClean:
        for obj in tests:
            if obj.buildDir == bdir:
                obj.reClean = 1

    return build_dirs

def copy_benchmarks(old_full_test_dir, full_web_dir, test_list, bench_dir, log):
    """ copy the last plotfile output from each test in test_list
        into the benchmark directory.  Also copy the diffDir, if
        it exists """
    td = os.getcwd()

    for t in test_list:
        wd = "{}/{}".format(old_full_test_dir, t.name)
        os.chdir(wd)

        if t.compareFile == "" and t.outputFile == "":
            p = t.get_last_plotfile(output_dir=wd)
        elif not t.outputFile == "":
            if not os.path.isdir(t.outputFile):
                p = get_recent_filename(wd, t.outputFile, ".tgz")
            else:
                p = t.outputFile
        else:
            if not os.path.isdir(t.compareFile):
                p = get_recent_filename(wd, t.compareFile, ".tgz")
            else:
                p = t.compareFile

        if not p == "":
            if p.endswith(".tgz"):
                try:
                    tg = tarfile.open(name=p, mode="r:gz")
                    tg.extractall()
                except:
                    log.fail("ERROR extracting tarfile")
                idx = p.rfind(".tgz")
                p = p[:idx]

            store_file = p
            if not t.outputFile == "":
                store_file = "{}_{}".format(t.name, p)

            try: shutil.rmtree("{}/{}".format(bench_dir, store_file))
            except: pass
            shutil.copytree(p, "{}/{}".format(bench_dir, store_file))

            with open("{}/{}.status".format(full_web_dir, t.name), 'w') as cf:
                cf.write("benchmarks updated.  New file:  {}\n".format(store_file) )

        else:   # no benchmark exists
            with open("{}/{}.status".format(full_web_dir, t.name), 'w') as cf:
                cf.write("benchmarks update failed")

        # is there a diffDir to copy too?
        if not t.diffDir == "":
            diff_dir_bench = "{}/{}_{}".format(bench_dir, t.name, t.diffDir)
            if os.path.isdir(diff_dir_bench):
                shutil.rmtree(diff_dir_bench)
                shutil.copytree(t.diffDir, diff_dir_bench)
            else:
                shutil.copy(t.diffDir, diff_dir_bench)
            log.log("new diffDir: {}_{}".format(t.name, t.diffDir))

        os.chdir(td)

def get_recent_filename(dir, base, extension):
    """ find the most recent file matching the base and extension """

    files = [f for f in os.listdir(dir) if (f.startswith(base) and
                                            f.endswith(extension))]

    files.sort(key=lambda x: os.path.getmtime(x))

    try: return files.pop()
    except: return None

def convert_to_f_make_flag(opt, test_not=False):
    if test_not:
        if opt: return " "
        else: return "t"
    else:
        if opt: return "t"
        else: return " "


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# test
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
def test_suite(argv):


    #--------------------------------------------------------------------------
    # parse the commandline arguments
    #--------------------------------------------------------------------------
    args=get_args(arg_string=argv)

    #--------------------------------------------------------------------------
    # read in the test information
    #--------------------------------------------------------------------------
    suite, test_list = load_params(args)

    active_test_list = [t.name for t in test_list]

    test_list = suite.get_tests_to_run(test_list)

    suite.log.skip()
    suite.log.bold("running tests: ")
    suite.log.indent()
    for obj in test_list:
        suite.log.log(obj.name)
    suite.log.outdent()

    if not args.complete_report_from_crash == "":

        # make sure the web directory from the crash run exists
        suite.full_web_dir = "%s/%s/"  % (suite.webTopDir, args.complete_report_from_crash)
        if not os.path.isdir(suite.full_web_dir):
            suite.log.fail("Crash directory does not exist")

        suite.test_dir = args.complete_report_from_crash

        # find all the tests that completed in that web directory
        tests = []
        test_file = ""
        was_benchmark_run = 0
        for file in os.listdir(suite.full_web_dir):
            if os.path.isfile(file) and file.endswith(".status"):
                index = string.rfind(file, ".status")
                tests.append(file[:index])

                with open(suite.full_web_dir + file, "r") as f:
                    for line in f:
                        if line.find("benchmarks updated") > 0:
                            was_benchmark_run = 1

            if os.path.isfile(file) and file.endswith(".ini"):
                test_file = file


        # create the report for this test run
        num_failed = report_this_test_run(suite, was_benchmark_run,
                                          "recreated report after crash of suite",
                                          "", tests, test_file)

        # create the suite report
        suite.log.bold("creating suite report...")
        report_all_runs(suite, active_test_list)
        suite.log.close_log()
        sys.exit("done")


    #--------------------------------------------------------------------------
    # check bench dir and create output directories
    #--------------------------------------------------------------------------
    all_compile = all([t.compileTest == 1 for t in test_list])

    if not all_compile:
        bench_dir = suite.get_bench_dir()

    if not args.copy_benchmarks is None:
        last_run = suite.get_last_run()

    suite.make_test_dirs()

    if suite.slack_post:
        msg = "{} ({}) test suite started, id: {}\n{}".format(suite.suiteName, suite.sub_title, suite.test_dir, args.note)
        suite.slack_post_it(msg)

    if not args.copy_benchmarks is None:
        old_full_test_dir = suite.testTopDir + suite.suiteName + "-tests/" + last_run
        copy_benchmarks(old_full_test_dir, suite.full_web_dir, test_list, bench_dir, suite.log)

        num_failed = report_this_test_run(suite, args.copy_benchmarks,   # plays the role of make_benchmarks here
                                          "copy_benchmarks used -- no new tests run",
                                          "",
                                          test_list, args.input_file[0])
        report_all_runs(suite, active_test_list)

        if suite.slack_post:
            msg = "copied benchmarks\n{}".format(args.copy_benchmarks)
            suite.slack_post_it(msg)

        sys.exit("done")


    #--------------------------------------------------------------------------
    # figure out what needs updating and do the git updates, save the
    # current hash / HEAD, and make a ChangeLog
    # --------------------------------------------------------------------------
    now = time.localtime(time.time())
    updateTime = time.strftime("%Y-%m-%d %H:%M:%S %Z", now)

    no_update = args.no_update.lower()
    if not args.copy_benchmarks is None:
        no_update = "all"

    # the default is to update everything, unless we specified a hash
    # when constructing the Repo object
    if no_update == "none":
        pass

    elif no_update == "all":
        for k in suite.repos:
            suite.repos[k].update = False

    else:
        nouplist = no_update.split(",")

        if "boxlib" in nouplist: suite.repos["BoxLib"].update = False
        if suite.srcName.lower() in nouplist: suite.repos["source"].update = False
        if suite.extSrcName.lower() in nouplist: suite.repos["extra_source"].update = False

        # each extra build directory has its own update flag
        for n, e in enumerate(suite.extra_build_names):
            if e.lower() in nouplist:
                suite.repos["extra_build-{}".format(n)].update = False

    os.chdir(suite.testTopDir)

    for k in suite.repos:
        suite.log.skip()
        suite.log.bold("repo: {}".format(suite.repos[k].name))
        suite.log.indent()

        if suite.repos[k].update or suite.repos[k].hash_wanted:
            suite.repos[k].git_update()

        suite.repos[k].save_head()

        if suite.repos[k].update:
            suite.repos[k].make_changelog()

        suite.log.outdent()

    #--------------------------------------------------------------------------
    # build the tools and do a make clean, only once per build directory
    #--------------------------------------------------------------------------
    suite.build_tools(test_list)

    all_build_dirs = find_build_dirs(test_list)

    suite.log.skip()
    suite.log.bold("make clean in...")

    for dir, source_tree in all_build_dirs:

        if not source_tree == "":
            suite.log.log("{} in {}".format(dir, source_tree))
            os.chdir(suite.repos[source_tree].dir + dir)
            suite.make_realclean(repo=source_tree)
        else:
            suite.log.log("{}".format(dir))
            os.chdir(suite.source_dir + dir)
            if suite.sourceTree == "BoxLib":
                suite.make_realclean(repo="BoxLib")
            else:
                suite.make_realclean()

    os.chdir(suite.testTopDir)


    #--------------------------------------------------------------------------
    # main loop over tests
    #--------------------------------------------------------------------------
    for test in test_list:

        suite.log.outdent()  # just to make sure we have no indentation
        suite.log.skip()
        suite.log.bold("working on test: {}".format(test.name))
        suite.log.indent()

        if not args.make_benchmarks == None and (test.restartTest or test.compileTest or
                                                 test.selfTest):
            suite.log.warn("  WARNING: test {} doesn't need benchmarks... skipping".format(test.name))
            continue

        output_dir = suite.full_test_dir + test.name + '/'
        os.mkdir(output_dir)
        test.output_dir = output_dir


        #----------------------------------------------------------------------
        # compile the code
        #----------------------------------------------------------------------
        if not test.extra_build_dir == "":
            bdir = suite.repos[test.extra_build_dir].dir + test.buildDir
        else:
            bdir = suite.source_dir + test.buildDir

        os.chdir(bdir)

        if test.reClean == 1:
            # for one reason or another, multiple tests use different
            # build options, make clean again to be safe
            suite.log.log("re-making clean...")
            if not test.extra_build_dir == "":
                suite.make_realclean(repo=test.extra_build_dir)
            else:
                suite.make_realclean()

        suite.log.log("building...")

        if suite.sourceTree == "C_Src" or test.testSrcTree == "C_Src":

            buildOptions = ""

            if test.debug:
                buildOptions += "DEBUG=TRUE "
            else:
                buildOptions += "DEBUG=FALSE "

            if test.useMPI:
                buildOptions += "USE_MPI=TRUE "
            else:
                buildOptions += "USE_MPI=FALSE "

            if test.useOMP:
                buildOptions += "USE_OMP=TRUE "
            else:
                buildOptions += "USE_OMP=FALSE "

            if not test.extra_build_dir == "":
                buildOptions += suite.repos[test.extra_build_dir].comp_string + " "

            comp_string = "{} -j{} BOXLIB_HOME={} {} {} DIM={} {} COMP={} FCOMP={} {}".format(
                suite.MAKE, suite.numMakeJobs, suite.boxlib_dir,
                suite.extra_src_comp_string, test.addToCompileString,
                test.dim, buildOptions, suite.COMP, suite.FCOMP,
                suite.add_to_c_make_command)

            suite.log.log(comp_string)
            so, se, r = run(comp_string,
                            outfile="{}/{}.make.out".format(output_dir, test.name))

            # get the executable
            executable = get_recent_filename(bdir, "", ".ex")
            
        elif suite.sourceTree == "F_Src" or test.testSrcTree == "F_Src":

            build_options = ""
            build_options += "NDEBUG={} ".format(convert_to_f_make_flag(test.debug, test_not=True))
            build_options += "MPI={} ".format(convert_to_f_make_flag(test.useMPI))
            build_options += "OMP={} ".format(convert_to_f_make_flag(test.useOMP))

            if not test.extra_build_dir == "":
                build_options += suite.repos[test.extra_build_dir].comp_string + " "

            comp_string = suite.build_f(opts="{} {} {}".format(
                suite.extra_src_comp_string, test.addToCompileString, build_options),
                                        outfile="{}/{}.make.out".format(output_dir, test.name))

            # get the executable
            executable = get_recent_filename(bdir, "main", ".exe")

        test.comp_string = comp_string

        if test.compileTest:

            # compilation tests are done now -- just make the report and ...
            shutil.copy("%s/%s.make.out"    % (output_dir, test.name), suite.full_web_dir)

            suite.log.log("creating problem test report ...")
            report_single_test(suite, test, test_list)

            # ... skip to the next test in the loop
            continue


        #----------------------------------------------------------------------
        # copy the necessary files over to the run directory
        #----------------------------------------------------------------------
        suite.log.log("copying files to run directory...")

        try: shutil.copy(executable, output_dir)
        except (IOError, AttributeError):

            # compilation failed.  First copy the make.out into the
            # web directory and then report
            shutil.copy("{}/{}.make.out".format(output_dir, test.name), suite.full_web_dir)

            error_msg = "ERROR: compilation failed"
            report_single_test(suite, test, test_list, failure_msg=error_msg)
            continue

        try: shutil.copy(test.inputFile, output_dir)
        except IOError:
            error_msg = "ERROR: unable to copy input file: {}".format(test.inputFile)
            report_single_test(suite, test, test_list, failure_msg=error_msg)
            continue

        # sometimes the input file was in a subdirectory under the
        # build directory.  Keep only the input file for latter
        index = string.rfind(test.inputFile, "/")
        if index > 0:
            test.inputFile = test.inputFile[index+1:]

        # if we are a "C_Src" build, we need the probin file
        if (suite.sourceTree == "C_Src" or \
                (test.testSrcTree == "C_Src" and test.probinFile != "")):
            try: shutil.copy(test.probinFile, output_dir)
            except IOError:
                error_msg = "ERROR: unable to copy probin file: {}".format(test.probinFile)
                report_single_test(suite, test, test_list, failure_msg=error_msg)
                continue

            # sometimes the probin file was in a subdirectory under the
            # build directory.  Keep only the probin file for latter
            index = string.rfind(test.probinFile, "/")
            if index > 0:
               test.probinFile = test.probinFile[index+1:]


        # python doesn't allow labelled continue statements, so we
        # use skip_to_next_test to decide if we need to skip to
        # the next test
        skip_to_next_test = 0
        for file in test.auxFiles:
            try: shutil.copy(file, output_dir)
            except IOError:
                error_msg = "ERROR: unable to copy aux file: {}".format(file)
                report_single_test(suite, test, test_list, failure_msg=error_msg)
                skip_to_next_test = 1
                break

        if skip_to_next_test: continue

        # python doesn't allow labelled continue statements, so we
        # use skip_to_next_test to decide if we need to skip to
        # the next test
        skip_to_next_test = 0
        for file in test.linkFiles:
            if not os.path.exists(file):
                error_msg = "ERROR: link file {} does not exist".format(file)
                report_single_test(suite, test, test_list, failure_msg=error_msg)
                skip_to_next_test = 1
                break

            else:
                if os.path.isabs(file):
                    link_source = file
                    link_name = output_dir + os.path.basename(file)
                else:
                    link_source = os.path.abspath(file)
                    link_name = output_dir + file
                try: os.symlink(link_source, link_name)
                except IOError:
                    error_msg = "ERROR: unable to symlink link file: {}".format(file)
                    report_single_test(suite, test, test_list, failure_msg=error_msg)
                    skip_to_next_test = 1
                    break

        if skip_to_next_test: continue


        #----------------------------------------------------------------------
        # run the test
        #----------------------------------------------------------------------
        suite.log.log("running the test...")

        os.chdir(output_dir)

        test.wall_time = time.time()

        if suite.sourceTree == "C_Src" or test.testSrcTree == "C_Src":

            base_command = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk" % \
                           (executable, test.inputFile, test.name, test.name)

            # keep around the checkpoint files only for the restart runs
            if test.restartTest:
                base_command += " amr.checkpoint_files_output=1 amr.check_int=%d" % \
                                (test.restartFileNum)
            else:
                base_command += " amr.checkpoint_files_output=0"

        elif suite.sourceTree == "F_Src" or test.testSrcTree == "F_Src":

            base_command = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk " % \
                           (executable, test.inputFile, test.name, test.name)

            # keep around the checkpoint files only for the restart runs
            if not test.restartTest: base_command += " --chk_int 0 "

            base_command += "{}".format(suite.globalAddToExecString)

        if args.with_valgrind:
            base_command = "valgrind " + args.valgrind_options + " " + base_command

        suite.run_test(test, base_command)


        # if it is a restart test, then rename the final output file and
        # restart the test
        if test.restartTest:
            skip_restart = False
            
            last_file = test.get_last_plotfile(output_dir=output_dir)

            if last_file == "":
                error_msg = "ERROR: test did not produce output.  Restart test not possible"
                skip_restart = True

            if len(test.find_backtrace()) > 0:
                error_msg = "ERROR: test produced backtraces.  Restart test not possible"
                skip_restart = True

            if skip_restart:
                # copy what we can
                test.wall_time = time.time() - test.wall_time
                shutil.copy("{}.run.out".format(test.name), suite.full_web_dir)
                shutil.copy("{}.make.out".format(test.name), suite.full_web_dir)
                suite.copy_backtrace(test)
                report_single_test(suite, test, test_list, failure_msg=error_msg)
                continue
            
            orig_last_file = "orig_%s" % (last_file)
            shutil.move(last_file, orig_last_file)

            if test.diffDir:
                origDiffDir = "orig_%s" % (test.diffDir)
                shutil.move(test.diffDir, origDiffDir)

            # get the file number to restart from
            restartFile = "%s_chk%5.5d" % (test.name, test.restartFileNum)

            suite.log.log("restarting from {} ... ".format(restartFile))

            if suite.sourceTree == "C_Src" or test.testSrcTree == "C_Src":

                base_command = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk amr.checkpoint_files_output=0 amr.restart=%s" % \
                        (executable, test.inputFile, test.name, test.name, restartFile)

            elif suite.sourceTree == "F_Src" or test.testSrcTree == "F_Src":

                base_command = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk --chk_int 0 --restart %d %s" % \
                        (executable, test.inputFile, test.name, test.name, test.restartFileNum, suite.globalAddToExecString)

            suite.run_test(test, base_command)

        test.wall_time = time.time() - test.wall_time


        #----------------------------------------------------------------------
        # do the comparison
        #----------------------------------------------------------------------
        if not test.selfTest:

            if test.outputFile == "":
                if test.compareFile == "":
                    compare_file = test.get_last_plotfile(output_dir=output_dir)
                else:
                    compare_file = test.compareFile
                output_file = compare_file
            else:
                output_file = test.outputFile
                compare_file = test.name+'_'+output_file


            # get the number of levels for reporting
            prog = "{} -l {}".format(suite.tools["fboxinfo"], output_file)
            stdout0, stderr0, rc = run(prog)
            test.nlevels = stdout0.rstrip('\n')
            if not type(convert_type(test.nlevels)) is int:
                test.nlevels = ""

            if args.make_benchmarks == None:

                suite.log.log("doing the comparison...")
                suite.log.indent()
                suite.log.log("comparison file: {}".format(output_file))

                test.compare_file_used = output_file

                if not test.restartTest:
                    bench_file = bench_dir + compare_file
                else:
                    bench_file = orig_last_file

                # see if it exists
                # note, with BoxLib, the plotfiles are actually directories

                if not os.path.isdir(bench_file):
                    suite.log.warn("WARNING: no corresponding benchmark found")
                    bench_file = ""

                    cf = open("%s.compare.out" % (test.name), 'w')
                    cf.write("WARNING: no corresponding benchmark found\n")
                    cf.write("         unable to do a comparison\n")
                    cf.close()

                else:
                    if not compare_file == "":

                        suite.log.log("benchmark file: {}".format(bench_file))

                        command = "{} -n 0 --infile1 {} --infile2 {}".format(
                            suite.tools["fcompare"], bench_file, output_file)
                        sout, serr, ierr = run(command, outfile="{}.compare.out".format(test.name), store_command=True)

                    else:
                        suite.log.warn("WARNING: unable to do a comparison")

                        cf = open("%s.compare.out" % (test.name), 'w')
                        cf.write("WARNING: run did not produce any output\n")
                        cf.write("         unable to do a comparison\n")
                        cf.close()

                suite.log.outdent()

                if not test.diffDir == "":
                    if not test.restartTest:
                        diff_dir_bench = bench_dir + '/' + test.name + '_' + test.diffDir
                    else:
                        diff_dir_bench = origDiffDir

                    suite.log.log("doing the diff...")
                    suite.log.log("diff dir: {}".format(test.diffDir))

                    command = "diff %s -r %s %s" \
                        % (test.diffOpts, diff_dir_bench, test.diffDir)

                    outfile = "{}.compare.out".format(test.name)
                    sout, serr, diff_status = run(command, outfile=outfile, store_command=True)

                    if diff_status == 0:
                        with open("{}.compare.out".format(test.name), 'a') as cf:
                            cf.write("\ndiff was SUCCESSFUL\n")

            else:   # make_benchmarks

                suite.log.log("storing output of {} as the new benchmark...".format(test.name))
                suite.log.indent()
                suite.log.warn("new benchmark file: {}".format(compare_file))
                suite.log.outdent()

                if not compare_file == "":
                    if not output_file == compare_file:
                        source_file = output_file
                    else:
                        source_file = compare_file

                    try: shutil.rmtree("{}/{}".format(bench_dir, compare_file))
                    except: pass
                    shutil.copytree(source_file, "{}/{}".format(bench_dir, compare_file))

                    with open("%s.status" % (test.name), 'w') as cf:
                        cf.write("benchmarks updated.  New file:  %s\n" % (compare_file) )

                else:
                    with open("%s.status" % (test.name), 'w') as cf:
                        cf.write("benchmarks failed")

                    # copy what we can
                    shutil.copy("{}.run.out".format(test.name), suite.full_web_dir)
                    shutil.copy("{}.make.out".format(test.name), suite.full_web_dir)
                    suite.copy_backtrace(test)
                    error_msg = "ERROR: runtime failure during benchmark creation"
                    report_single_test(suite, test, test_list, failure_msg=error_msg)


                if not test.diffDir == "":
                    diff_dir_bench = "{}/{}_{}".format(bench_dir, test.name, test.diffDir)
                    if os.path.isdir(diff_dir_bench):
                        shutil.rmtree(diff_dir_bench)
                        shutil.copytree(test.diffDir, diff_dir_bench)
                    else:
                        shutil.copy(test.diffDir, diff_dir_bench)
                    suite.log.log("new diffDir: {}_{}".format(test.name, test.diffDir))

        else:   # selfTest

            if args.make_benchmarks == None:

                suite.log.log("looking for selfTest success string: {} ...".format(test.stSuccessString))

                try: of = open("{}.run.out".format(test.name), 'r')
                except IOError:
                    suite.log.warn("WARNING: no output file found")
                    compare_successful = 0
                    outLines = ['']
                else:
                    outLines = of.readlines()

                    # successful comparison is indicated by PLOTFILES AGREE
                    compare_successful = 0

                    for line in outLines:
                        if line.find(test.stSuccessString) >= 0:
                            compare_successful = 1
                            break

                    of.close()

                with open("%s.compare.out" % (test.name), 'w') as cf:
                    if compare_successful:
                        cf.write("SELF TEST SUCCESSFUL\n")
                    else:
                        cf.write("SELF TEST FAILED\n")


        #----------------------------------------------------------------------
        # do any requested visualization (2- and 3-d only) and analysis
        #----------------------------------------------------------------------
        if output_file != "":
            if args.make_benchmarks == None:

                # get any parameters for the summary table
                job_info_file = "{}/job_info".format(output_file)
                if os.path.isfile(job_info_file):
                    test.has_jobinfo = 1
                
                try: jif = open(job_info_file, "r")
                except:
                    suite.log.warn("unable to open the job_info file")
                else:
                    job_file_lines = jif.readlines()

                    if suite.summary_job_info_field1 is not "":
                        for l in job_file_lines:
                            if l.find(suite.summary_job_info_field1) >= 0 and l.find(":") >= 0:
                                _tmp = l.split(":")[1]
                                idx = _tmp.rfind("/") + 1
                                test.job_info_field1 = _tmp[idx:]
                                break

                    if suite.summary_job_info_field2 is not "":
                        for l in job_file_lines:
                            if l.find(suite.summary_job_info_field2) >= 0 and l.find(":") >= 0:
                                _tmp = l.split(":")[1]
                                idx = _tmp.rfind("/") + 1
                                test.job_info_field2 = _tmp[idx:]
                                break

                    if suite.summary_job_info_field3 is not "":
                        for l in job_file_lines:
                            if l.find(suite.summary_job_info_field3) >= 0 and l.find(":") >= 0:
                                _tmp = l.split(":")[1]
                                idx = _tmp.rfind("/") + 1
                                test.job_info_field3 = _tmp[idx:]
                                break

                # visualization
                if test.doVis:

                    if test.dim == 1:
                        suite.log.log("Visualization not supported for dim = {}".format(test.dim))
                    else:
                        suite.log.log("doing the visualization...")
                        tool = suite.tools["fsnapshot{}d".format(test.dim)]
                        run('{} --palette {}/Palette -cname "{}" -p "{}"'.format(
                            tool, suite.compare_tool_dir, test.visVar, output_file))

                        # convert the .ppm files into .png files
                        ppm_file = get_recent_filename(output_dir, "", ".ppm")
                        if not ppm_file is None:
                            png_file = ppm_file.replace(".ppm", ".png")
                            run("convert {} {}".format(ppm_file, png_file))
                            test.png_file = png_file

                # analysis
                if not test.analysisRoutine == "":

                    suite.log.log("doing the analysis...")
                    if not test.extra_build_dir == "":
                        tool = "{}/{}".format(suite.repos[test.extra_build_dir].dir, test.analysisRoutine)
                    else:
                        tool = "{}/{}".format(suite.source_dir, test.analysisRoutine)

                    shutil.copy(tool, os.getcwd())

                    option = eval("suite.{}".format(test.analysisMainArgs))
                    run("{} {} {}".format(os.path.basename(test.analysisRoutine),
                                          option, output_file))

        else:
            if test.doVis or test.analysisRoutine != "":
                suite.log.warn("WARNING: no output file.  Skipping visualization")


        #----------------------------------------------------------------------
        # move the output files into the web directory
        #----------------------------------------------------------------------
        if args.make_benchmarks == None:
            shutil.copy("{}.run.out".format(test.name), suite.full_web_dir)
            shutil.copy("{}.make.out".format(test.name), suite.full_web_dir)
            shutil.copy("{}.compare.out".format(test.name), suite.full_web_dir)

            shutil.copy(test.inputFile, "{}/{}.{}".format(
                suite.full_web_dir, test.name, test.inputFile) )

            if test.has_jobinfo:
                shutil.copy(job_info_file, "{}/{}.job_info".format(
                    suite.full_web_dir, test.name))
                
            if suite.sourceTree == "C_Src":
                shutil.copy(test.probinFile, "{}/{}.{}".format(
                    suite.full_web_dir, test.name, test.probinFile) )

            for af in test.auxFiles:

                # sometimes the auxFile was in a subdirectory under the
                # build directory.
                index = string.rfind(af, "/")
                if index > 0:
                    af = af[index+1:]

                shutil.copy(af, "{}/{}.{}".format(suite.full_web_dir, test.name, af) )

            if not test.png_file is None:
                try: shutil.copy(test.png_file, suite.full_web_dir)
                except IOError:
                    # visualization was not successful.  Reset image
                    test.png_file = None

            if not test.analysisRoutine == "":
                try: shutil.copy(test.analysisOutputImage, suite.full_web_dir)
                except IOError:
                    # analysis was not successful.  Reset the output image
                    test.analysisOutputImage = ""

            # were any Backtrace files output (indicating a crash)
            suite.copy_backtrace(test)

        else:
            shutil.copy("%s.status" % (test.name), suite.full_web_dir)


        #----------------------------------------------------------------------
        # archive (or delete) the output
        #----------------------------------------------------------------------
        suite.log.log("archiving the output...")
        for file in os.listdir(output_dir):
            if (os.path.isdir(file) and
                (file.startswith("%s_plt" % (test.name)) or
                 file.startswith("%s_chk" % (test.name)) ) ):

                if suite.purge_output == 1 and not file == output_file:
                    # delete the plt/chk file
                    if os.path.isdir(file):
                        try: shutil.rmtree(file)
                        except:
                            suite.log.warn("WARNING: unable to remove {}".format(file))

                else:
                    # tar it up
                    try:
                        tar = tarfile.open("%s.tgz" % (file), "w:gz")
                        tar.add("%s" % (file))
                        tar.close()

                    except:
                        suite.log.warn("WARNING: unable to tar output file %s" % (file))

                    else:
                        shutil.rmtree(file)


        #----------------------------------------------------------------------
        # write the report for this test
        #----------------------------------------------------------------------
        if args.make_benchmarks == None:
            suite.log.log("creating problem test report ...")
            report_single_test(suite, test, test_list)


    #--------------------------------------------------------------------------
    # write the report for this instance of the test suite
    #--------------------------------------------------------------------------
    suite.log.outdent()
    suite.log.skip()
    suite.log.bold("creating new test report...")
    num_failed = report_this_test_run(suite, args.make_benchmarks, args.note,
                                      updateTime,
                                      test_list, args.input_file[0])


    # make sure that all of the files in the web directory are world readable
    for file in os.listdir(suite.full_web_dir):
       current_file = suite.full_web_dir + file

       if os.path.isfile(current_file):
          os.chmod(current_file, 0o644)

    # reset the branch to what it was originally
    suite.log.skip()
    suite.log.bold("reverting git branches/hashes")
    suite.log.indent()

    for k in suite.repos:
        if suite.repos[k].update or suite.repos[k].hash_wanted:
            suite.repos[k].git_back()

    suite.log.outdent()

    # For temporary run, return now without creating suote report.
    if args.do_temp_run:
        return num_failed


    # store an output file in the web directory that can be parsed easily by
    # external program
    name = "source"
    if suite.sourceTree == "BoxLib": name = "BoxLib"
    branch = suite.repos[name].branch_wanted.strip("\"")

    with open("{}/suite.{}.status".format(suite.webTopDir, branch), "w") as f:
        f.write("{}; num failed: {}; source hash: {}".format(
            suite.repos[name].name, num_failed, suite.repos[name].hash_current))


    #--------------------------------------------------------------------------
    # generate the master report for all test instances
    #--------------------------------------------------------------------------
    suite.log.skip()
    suite.log.bold("creating suite report...")
    report_all_runs(suite, active_test_list)

    def emailDevelopers():
        msg = email.message_from_string(suite.emailBody)
        msg['From'] = suite.emailFrom
        msg['To'] = ",".join(suite.emailTo)
        msg['Subject'] = suite.emailSubject

        server = smtplib.SMTP('localhost')
        server.sendmail(suite.emailFrom, suite.emailTo, msg.as_string())
        server.quit()

    if num_failed > 0 and suite.sendEmailWhenFail and not args.send_no_email:
        suite.log.skip()
        suite.log.bold("sending email...")
        emailDevelopers()


    if suite.slack_post:
        suite.slack_post_it("test complete, num failed = {}\n{}".format(num_failed, suite.emailBody))

    return num_failed


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# R E P O R T   W R I T I N G   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

cssContents = \
r"""
body {font-family: "Arial", san-serif;}

h1 {font-family: "Tahoma","Arial", sans-serif;
    color: #333333;}

h3 {display: inline;}

h3.passed {text-decoration: none; display: inline;
           color: black; background-color: lime; padding: 2px;}

a.passed:link {color: black; text-decoration: none;}
a.passed:visited {color: black; text-decoration: none;}
a.passed:hover {color: #ee00ee; text-decoration: underline;}

h3.failed {text-decoration: none; display: inline;
           color: yellow; background-color: red; padding: 2px;}

a.failed:link {color: yellow; text-decoration: none;}
a.failed:visited {color: yellow; text-decoration: none;}
a.failed:hover {color: #00ffff; text-decoration: underline;}

h3.benchmade {text-decoration: none; display: inline;
              color: black; background-color: orange; padding: 2px;}

a.benchmade:link {color: black; text-decoration: none;}
a.benchmade:visited {color: black; text-decoration: none;}
a.benchmade:hover {color: #00ffff; text-decoration: underline;}


span.nobreak {white-space: nowrap;}

a.main:link {color: yellow; text-decoration: none;}
a.main:visited {color: yellow; text-decoration: none;}
a.main:hover {color: #00ffff; text-decoration: underline;}

td {border-width: 0px;
    padding: 5px;
    background-color: white;
    vertical-align: middle;}

td.passed {background-color: lime; opacity: 0.8;}
td.failed {background-color: red; color: yellow; opacity: 0.8;}
td.benchmade {background-color: orange; opacity: 0.8;}
td.date {background-color: #666666; color: white; opacity: 0.8; font-weight: bold;}

.maintable tr:hover {background-color: blue;}


table {border-collapse: separate;
       border-spacing: 2px;
       margin-left: auto;
       margin-right: auto;
       border-width: 1px;
       border-color: gray;
       border-style: solid;
       box-shadow: 10px 10px 5px #888888;}

table.head {border-collapse: separate;
       border-spacing: 0px;
       margin-left: auto;
       margin-right: auto;
       border-width: 0px;
       border-style: solid;
       box-shadow: none;}

/* http://blog.petermares.com/2010/10/27/vertical-text-in-html-table-headers-for-webkitmozilla-browsers-without-using-images/ */

div.verticaltext {text-align: center;
                  vertical-align: middle;
                  width: 20px;
                  margin: 0px;
                  padding: 0px;
                  padding-left: 3px;
                  padding-right: 3px;
                  padding-top: 10px;
                  white-space: nowrap;
                  -webkit-transform: rotate(-90deg);
                  -moz-transform: rotate(-90deg);}

#summary th {background-color: grey;
    color: yellow;
    text-align: center;
    height: 2em;
    padding-bottom: 3px;
    padding-left: 5px;
    padding-right: 5px;}


#summary td {background: transparent;}

#summary tr:nth-child(even) {background: #dddddd;}
#summary tr:nth-child(odd) {background: #eeeeee;}

#summary tr.special {background: #ccccff;}
#summary td.highlight {color: red;}

#summary td.passed {background-color: lime;}
#summary td.failed {background-color: red;}
#summary td.benchmade {background-color: orange;}

div.small {font-size: 75%;}

th {background-color: grey;
    color: yellow;
    text-align: center;
    vertical-align: bottom;
    height: @TABLEHEIGHT@;
    padding-bottom: 3px;
    padding-left: 5px;
    padding-right: 5px;}

li {padding-top: 0.5em;}

ul li {color: blue;
       font-weight: bold;}
ul li ul li {color: black;
             font-weight: normal;}

ul li h3 {border: 1px solid black;}

#compare td {font-family: "Lucida Console", Monaco, monospace;
             font-size: 80%;}

#box {  width: 900px;
  margin: 0 auto;
  padding: 1em;
  background: #ffffff;
}

.alignright {
   text-align: right;
}

"""

HTMLHeader = \
r"""
<HTML>
<HEAD>
<TITLE>@TESTDIR@ / @TESTNAME@</TITLE>
<META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=ISO-8859-1">
<LINK REL="stylesheet" TYPE="text/css" HREF="tests.css">
</HEAD>
<BODY>
<div id="box">
"""

MainHeader = \
r"""
<HTML>
<HEAD>
<TITLE>@TITLE@</TITLE>
<META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=ISO-8859-1">
<LINK REL="stylesheet" TYPE="text/css" HREF="tests.css">
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.5.0/css/font-awesome.min.css">
</HEAD>
<BODY>
<!--GOUPLINK-->
<CENTER><H1>@TITLE@</H1></CENTER>
<CENTER><H2>@SUBTITLE@</H2></CENTER>
"""


#==============================================================================
# create_css
#==============================================================================
def create_css(tableHeight=16):
    """ write the css file for the webpages """

    cssC = cssContents.replace("@TABLEHEIGHT@", "%sem" % (tableHeight))

    with open("tests.css", 'w') as cf:
        cf.write(cssC)

class HTMLList(object):
    def __init__(self, of=None):
        self.list_items = []
        self.current_indent = 0
        self.of = of

    def item(self, content):
        self.list_items.append((self.current_indent, content))

    def indent(self):
        self.current_indent += 1

    def outdent(self):
        self.current_indent -= 1

    def write_list(self):
        self.of.write("<ul>\n")
        current_indent = -1
        for i, c in self.list_items:
            if current_indent == -1:
                current_indent = i
            else:
                if i < current_indent:
                    self.of.write("</li></ul></li>\n")
                elif i > current_indent:
                    self.of.write("<ul>\n")
                else:
                    self.of.write("</li>\n")

            current_indent = i
            self.of.write("<li>{}\n".format(c))

        # finish current item
        self.of.write("</li>")

        # finish nesting
        for n in range(0,current_indent):
            self.of.write("</ul></li>\n")

        self.of.write("</ul>\n")

class HTMLTable(object):
    def __init__(self, out_file, columns=1, divs=None):
        self.hf = out_file
        self.columns = columns
        if not divs is None:
            self.divs = list(divs)
        else:
            self.divs = None

    def start_table(self):
        if not self.divs is None:
            for d in self.divs:
                self.hf.write("<div id='{}'>\n".format(d))
        self.hf.write("<p><table>\n")

    def header(self, header_list):
        n = len(header_list)
        line = "<tr>"+n*"<th>{}</th>"+"</tr>\n"
        self.hf.write(line.format(*header_list))

    def print_single_row(self, row):
        self.hf.write("<tr class='special'><td colspan={}>".format(self.columns)+row+"</td></tr>\n")

    def print_row(self, row_list, highlight=False):
        """ row_list are the individual table elements.  Note that if
        a list item is a tuple, then the first element is assumed to
        be the cell data and the second element is an html tag that
        goes in the <td >, e.g. to set the class or colspan"""

        n = len(row_list)
        if highlight:
            line = "<tr>"+n*"<td class='highlight'>{}</td>"+"</tr>\n"
        else:
            line = "<tr>"
            for d in row_list:
                if isinstance(d, tuple):
                    line += "<td {}>{}</td>".format(d[1], d[0])
                else:
                    line += "<td>{}</td>".format(d)
            line += "</tr>\n"
        self.hf.write(line.format(*row_list))

    def end_table(self):
        self.hf.write("</table>\n")
        if not self.divs is None:
            for n in range(len(self.divs)):
                self.hf.write("</div>\n")


#==============================================================================
# REPORT ROUTINES
#==============================================================================
def report_single_test(suite, test, tests, failure_msg=None):
    """ generate a single problem's test result page.  If
        failure_msg is set to a string, then it is assumed
        that the test did not complete.  The string will
        be reported on the test page as the error. """

    # for navigation
    tnames = [t.name for t in tests]
    current_index = tnames.index(test.name)

    if not failure_msg is None:
        suite.log.testfail("aborting test")
        suite.log.testfail(failure_msg)

    # get the current directory
    current_dir = os.getcwd()

    # switch to the web directory and open the report file
    os.chdir(suite.full_web_dir)

    #--------------------------------------------------------------------------
    # parse the compilation report and determine if we compiled
    #--------------------------------------------------------------------------
    compileFile = "%s.make.out" % (test.name)

    try: cf = open(compileFile, 'r')
    except IOError:
        suite.log.warn("WARNING: no compilation file found")
        compile_successful = 0
    else:
        # successful compilation be indicated by SUCCESS or
        # Nothing to be done for `all'.  Look for both
        compile_successful = 0

        for line in cf:
            if (line.find("SUCCESS") >= 0 or
                line.find("is up to date.") >= 0 or
                line.find("Nothing to be done") >= 0):
                compile_successful = 1
                break

        cf.close()


    #--------------------------------------------------------------------------
    # parse the compare report and determine if we passed
    #--------------------------------------------------------------------------
    if failure_msg is None:
        if not test.compileTest:
            compare_file = "{}.compare.out".format(test.name)

            try: cf = open(compare_file, 'r')
            except IOError:
                suite.log.warn("WARNING: no comparison file found")
                compare_successful = 0
                diff_lines = ['']
            else:
                diff_lines = cf.readlines()

                # successful comparison is indicated by PLOTFILES AGREE
                compare_successful = 0
                for line in diff_lines:
                    if (line.find("PLOTFILES AGREE") >= 0 or
                        line.find("SELF TEST SUCCESSFUL") >= 0):
                        compare_successful = 1
                        break
                    
                if compare_successful:
                    if not test.diffDir == "":
                        compare_successful = 0
                        for line in diff_lines:
                            if line.find("diff was SUCCESSFUL") >= 0:
                                compare_successful = 1
                                break

                cf.close()

                # last check: did we produce any backtrace files?
                if len(test.backtrace) > 0: compare_successful = 0

        # write out the status file for this problem, with either
        # PASSED or FAILED
        status_file = "{}.status".format(test.name)
        with open(status_file, 'w') as sf:
            if (compile_successful and
                (test.compileTest or (not test.compileTest and compare_successful))):
                sf.write("PASSED\n")
                suite.log.success("{} PASSED".format(test.name))
            else:
                sf.write("FAILED\n")
                suite.log.testfail("{} FAILED".format(test.name))

    else:
        # we came in already admitting we failed...
        status_file = "{}.status".format(test.name)
        with open(status_file, 'w') as sf:
            sf.write("FAILED\n")
        suite.log.testfail("{} FAILED".format(test.name))


    #--------------------------------------------------------------------------
    # generate the HTML page for this test
    #--------------------------------------------------------------------------

    # write the css file
    create_css()

    html_file = "{}.html".format(test.name)
    hf = open(html_file, 'w')

    new_head = HTMLHeader

    # arrows for previous and next test
    new_head += r"""<table style="width: 100%" class="head"><br><tr>"""
    if current_index > 0:
        new_head += r"""<td><< <a href="{}.html">previous test</td>""".format(tests[current_index-1].name)
    else:
        new_head += r"""<td>&nbsp;</td>"""

    if current_index < len(tests)-1:
        new_head += r"""<td class="alignright"><a href="{}.html">next test >></td>""".format(tests[current_index+1].name)
    else:
        new_head += r"""<td>&nbsp;</td>"""

    new_head += r"</tr></table>" + "\n"


    new_head += r"""<center><h1><a href="index.html">@TESTDIR@</a> / @TESTNAME@</h1></center>"""

    new_head = new_head.replace("@TESTDIR@", os.path.normpath(suite.test_dir))
    new_head = new_head.replace("@TESTNAME@", test.name)

    hf.write(new_head)


    ll = HTMLList(of=hf)

    if not failure_msg is None:
        ll.item("Test error: ")
        ll.indent()

        ll.item("<h3 class=\"failed\">Failed</h3>")
        ll.item("{}".format(failure_msg))

        ll.outdent()

    # build summary
    ll.item("Build/Test information:")
    ll.indent()

    ll.item("Build directory: {}".format(test.buildDir))

    if not test.extra_build_dir == "":
        ll.indent()
        ll.item("in {}".format(suite.repos[test.extra_build_dir].dir))
        ll.outdent()

    if not test.compileTest:

        if test.debug:
            ll.item("Debug test")

        if test.useMPI or test.useOMP:
            ll.item("Parallel run")
            ll.indent()
            if test.useMPI:
                ll.item("MPI numprocs = {}".format(test.numprocs))
            if test.useOMP:
                ll.item("OpenMP numthreads = {}".format(test.numthreads))
            ll.outdent()

        if test.restartTest:

            ll.item("Restart test")
            ll.indent()
            ll.item("Job was run as normal and then restarted from checkpoint # {}, and the two final outputs were compared".format(test.restartFileNum))
            ll.outdent()


        ll.item("Files:")
        ll.indent()

        ll.item("input file: <a href=\"{}.{}\">{}</a>".format(test.name, test.inputFile, test.inputFile))

        if suite.sourceTree == "C_Src":
            ll.item("probin file: <a href=\"{}.{}\">{}</a>".format(test.name, test.probinFile, test.probinFile))

        for i, afile in enumerate(test.auxFiles):
            # sometimes the auxFile was in a subdirectory under the
            # build directory.
            index = string.rfind(afile, "/")
            if index > 0:
                root_file = afile[index+1:]
            else:
                root_file = afile

            ll.item("auxillary file {}: <a href=\"{}.{}\">{}</a>".format(i+1, test.name, root_file, afile))

        ll.outdent()

        ll.item("Dimensionality: {}".format(test.dim))

    ll.outdent()   # end of build information

    # compilation summary
    ll.item("Compilation:")
    ll.indent()

    if compile_successful:
        ll.item("<h3 class=\"passed\">Successful</h3>")
    else:
        ll.item("<h3 class=\"failed\">Failed</h3>")

    ll.item("Compilation command:<br><tt>{}</tt>".format(test.comp_string))
    ll.item("<a href=\"{}.make.out\">make output</a>".format(test.name))

    ll.outdent()


    if not test.compileTest:

        # execution summary
        ll.item("Execution:")
        ll.indent()
        ll.item("Execution time: {:.3f} s".format(test.wall_time))
        ll.item("Execution command:<br><tt>{}</tt>".format(test.run_command))
        ll.item("<a href=\"{}.run.out\">execution output</a>".format(test.name))
        if test.has_jobinfo:
            ll.item("<a href=\"{}.job_info\">job_info</a>".format(test.name))
        ll.outdent()


        # were there backtrace files?
        if len(test.backtrace) > 0:
            ll.item("Backtraces:")
            ll.indent()
            for bt in test.backtrace:
                ll.item("<a href=\"{}\">{}</a>".format(bt, bt))
            ll.outdent()

        # comparison summary
        if failure_msg is None:
            ll.item("Comparison: ")
            ll.indent()

            if compare_successful:
                ll.item("<h3 class=\"passed\">Successful</h3>")
            else:
                ll.item("<h3 class=\"failed\">Failed</h3>")

    ll.write_list()

    if (not test.compileTest) and failure_msg is None:

        # parse the compare output and make an HTML table
        ht = HTMLTable(hf, columns=3, divs=["summary", "compare"])
        in_diff_region = False

        box_error = False
        grid_error = False
        
        for line in diff_lines:

            if "number of boxes do not match" in line:
                box_error = True
                break

            if "grids do not match" in line:
                grid_error = True
                break

            if not in_diff_region:
                if line.find("fcompare") > 1:
                    hf.write("<tt>"+line+"</tt>\n")

                    ht.start_table()
                    continue

                if line.strip().startswith("diff"):
                    ht.end_table()
                    hf.write("<pre>\n")

                    hf.write(line)
                    in_diff_region = True
                    continue

                if line.strip().startswith("level"):
                    ht.print_single_row(line.strip())
                    continue

                if line.strip().startswith("-----"):
                    continue

                if line.strip().startswith("<<<"):
                    ht.print_single_row(line.strip().replace('<','&lt;').replace('>','&gt;'))
                    continue

                fields = [q.strip() for q in line.split("  ") if not q == ""]

                if fields[0].startswith("variable"):
                    ht.header(fields)
                    continue

                if len(fields) == 2:
                    if "NaN present" in line:
                        ht.print_row([fields[0], (fields[1], "colspan='2'")])
                        continue
                    elif "variable not present" in line:
                        ht.print_row([fields[0], (fields[1], "colspan='2'")])
                        continue
                    else:
                        ht.header([" "] + fields)
                        continue

                if len(fields) == 1:
                    continue

                else:
                    abs_err = float(fields[1])
                    rel_err = float(fields[2])
                    if abs(rel_err) > 1.e-6:
                        ht.print_row([fields[0], abs_err, rel_err], highlight=True)
                    else:
                        ht.print_row([fields[0], abs_err, rel_err])

            else:
                # diff region
                hf.write(line)

        if in_diff_region:
            hf.write("</pre>\n")
        else:
            ht.end_table()

        if box_error:
            hf.write("<p>number of boxes do not match</p>\n")

        if grid_error:
            hf.write("<p>grids do not match</p>\n")            

        # show any visualizations
        if test.doVis:
            if not test.png_file is None:
                hf.write("<P>&nbsp;\n")
                hf.write("<P><IMG SRC='{}' BORDER=0>".format(test.png_file))

        # show any analysis
        if not test.analysisOutputImage == "":
            hf.write("<P>&nbsp;\n")
            hf.write("<P><IMG SRC='%s' BORDER=0>" % (test.analysisOutputImage) )


    # close
    hf.write("</div></body>\n")
    hf.write("</html>\n")

    hf.close()


    # switch back to the original directory
    os.chdir(current_dir)



def report_this_test_run(suite, make_benchmarks, note, update_time,
                         test_list, test_file):
    """ generate the master page for a single run of the test suite """

    # get the current directory
    current_dir = os.getcwd()

    # switch to the web directory and open the report file
    os.chdir(suite.full_web_dir)


    # keep track of the number of tests that passed and the number that failed
    num_failed = 0
    numPassed = 0


    #--------------------------------------------------------------------------
    # generate the HTML page for this run of the test suite
    #--------------------------------------------------------------------------

    # always create the css (in case it changes)
    create_css()

    # create the master filename
    htmlFile = "index.html"

    hf = open(htmlFile, 'w')

    new_head = HTMLHeader + r"""<CENTER><H1><A HREF="../">@TESTDIR@</A> / @TESTNAME@</H1></CENTER>"""

    new_head = new_head.replace("@TESTDIR@", suite.suiteName)
    new_head = new_head.replace("@TESTNAME@", suite.test_dir)

    hf.write(new_head)

    if not note == "":
       hf.write("<p><b>Test run note:</b><br><font color=\"gray\">%s</font>\n" % (note) )

    if not make_benchmarks == None:
       hf.write("<p><b>Benchmarks updated</b><br>comment: <font color=\"gray\">{}</font>\n".format(make_benchmarks) )
       hf.write("<p>&nbsp;\n")


    hf.write("<p><b>test input parameter file:</b> <A HREF=\"%s\">%s</A>\n" %
             (test_file, test_file) )

    # git info lists
    any_update = any([suite.repos[t].update for t in suite.repos])

    if any_update and not update_time == "":
        hf.write("<p>&nbsp;\n")
        hf.write("<p><b>Git update was done at: </b>%s\n" % (update_time) )

        hf.write("<ul>\n")
        code_str = "<li><b>{}</b><ul><li><b>branch:</b> {}; <b>hash:</b> {}</li><li><b>changelog:</b> <a href=\"{}\">{}</a></li></ul></li>"

        for k, r in suite.repos.items():
            if r.update:
                hf.write(code_str.format(r.name, r.branch_wanted, r.hash_current,
                                         "ChangeLog.{}".format(r.name),
                                         "ChangeLog.{}".format(r.name)))

        hf.write("</ul>")

    else:
        hf.write("<p>No git update done\n")

    hf.write("<p>&nbsp;\n")

    # summary table
    if make_benchmarks == None:
        special_cols = []
        if suite.summary_job_info_field1 is not "":
            special_cols.append(suite.summary_job_info_field1)
        if suite.summary_job_info_field2 is not "":
            special_cols.append(suite.summary_job_info_field2)            
        if suite.summary_job_info_field3 is not "":
            special_cols.append(suite.summary_job_info_field3)            

        cols = ["test name", "dim", "compare plotfile",
                "# levels", "MPI procs", "OMP threads", "debug",
                "compile", "restart"] + special_cols + ["wall time", "result"]
        ht = HTMLTable(hf, columns=len(cols), divs=["summary"])
        ht.start_table()
        ht.header(cols)

    else:
        ht = HTMLTable(hf, columns=3, divs=["summary"])
        ht.start_table()
        ht.header(["test name", "result", "comment"])

    # loop over the tests and add a line for each
    for test in test_list:

        if make_benchmarks == None:

            # check if it passed or failed
            status_file = "%s.status" % (test.name)

            testPassed = 0

            with open(status_file, 'r') as sf:
                for line in sf:
                    if line.find("PASSED") >= 0:
                        testPassed = 1
                        numPassed += 1
                        break

                if not testPassed:
                    num_failed += 1

            row_info = []
            row_info.append("<a href=\"{}.html\">{}</a>".format(test.name, test.name))
            row_info.append(test.dim)
            row_info.append("<div class='small'>{}</div>".format(test.compare_file_used))

            if not test.nlevels == None:
                row_info.append(test.nlevels)
            else:
                row_info.append("")

            if test.useMPI:
                row_info.append("&check; ({})".format(test.numprocs))
            else:
                row_info.append("")

            # OMP ?
            if test.useOMP:
                row_info.append("&check; ({})".format(test.numthreads))
            else:
                row_info.append("")

            # debug ?
            if test.debug:
                row_info.append("&check;")
            else:
                row_info.append("")

            # compile ?
            if test.compileTest:
                row_info.append("&check;")
            else:
                row_info.append("")

            # restart ?
            if test.restartTest:
                row_info.append("&check;")
            else:
                row_info.append("")


            # special columns
            if suite.summary_job_info_field1 is not "":
                row_info.append("<div class='small'>{}</div>".format(test.job_info_field1))

            if suite.summary_job_info_field2 is not "":
                row_info.append("<div class='small'>{}</div>".format(test.job_info_field2))

            if suite.summary_job_info_field3 is not "":
                row_info.append("<div class='small'>{}</div>".format(test.job_info_field3))

                
            # wallclock time
            row_info.append("{:.3f}&nbsp;s".format(test.wall_time))

            if testPassed:
                row_info.append(("PASSED", "class='passed'"))
            else:
                row_info.append(("FAILED", "class='failed'"))

            ht.print_row(row_info)

        else:
            if test.restartTest: continue
            if test.compileTest: continue
            if test.selfTest: continue

            # the benchmark was updated -- find the name of the new benchmark file
            benchStatusFile = "%s.status" % (test.name)

            bench_file = "none"

            with open(benchStatusFile, 'r') as bf:
                for line in bf:
                    index = line.find("file:")
                    if index >= 0:
                        bench_file = line[index+5:]
                        break

            row_info = []
            row_info.append("{}".format(test.name))
            if not bench_file == "none":
                row_info.append(("BENCHMARK UPDATED", "class='benchmade'"))
                row_info.append("new benchmark file is {}".format(bench_file))
            else:
                row_info.append(("BENCHMARK NOT UPDATED", "class='failed'"))
                row_info.append("compilation or execution failed")

            ht.print_row(row_info)

    ht.end_table()

    # close
    hf.write("</div></body>\n")
    hf.write("</html>\n")
    hf.close()


    #--------------------------------------------------------------------------
    # write out a status file for all the tests
    #--------------------------------------------------------------------------

    status_file = os.path.normpath(suite.test_dir) + ".status"
    with open(status_file, 'w') as sf:

        if make_benchmarks == None:
            if num_failed == 0:
                sf.write("ALL PASSED\n")
            elif num_failed > 0 and numPassed > 0:
                sf.write("SOME FAILED\n")
            else:
                sf.write("ALL FAILED\n")

        else:
            sf.write("BENCHMARKS UPDATED\n")

    # switch back to the original directory
    os.chdir(current_dir)

    return num_failed


def report_all_runs(suite, active_test_list):

    tableHeight = min(max(suite.lenTestName, 4), 18)

    os.chdir(suite.webTopDir)

    create_css(tableHeight=tableHeight)

    valid_dirs, all_tests = suite.get_run_history(active_test_list)

    if do_timings_plots: suite.make_timing_plots(active_test_list)

    #--------------------------------------------------------------------------
    # generate the HTML
    #--------------------------------------------------------------------------
    htmlFile = "index.html"

    title = "%s regression tests" % (suite.suiteName)

    hf = open(htmlFile, 'w')

    header = MainHeader.replace("@TITLE@", title).replace("@SUBTITLE@", suite.sub_title)

    if suite.goUpLink:
        header2 = header.replace("<!--GOUPLINK-->", '<a href="../">GO UP</a>')
        hf.write(header2)
    else:
        hf.write(header)

    hf.write("<P><TABLE class='maintable'>\n")

    # write out the header
    hf.write("<TR><TH ALIGN=CENTER>date</TH>\n")
    for test in all_tests:
        hf.write("<TH><div class='verticaltext'>%s</div></TH>\n" % (test))

    hf.write("</TR>\n")


    if do_timings_plots:
        hf.write("<tr><td class='date'>plots</td>")
        for t in all_tests:
            plot_file = "{}-timings.png".format(t)
            if os.path.isfile(plot_file):
                hf.write("<TD ALIGN=CENTER title=\"{} timings plot\"><H3><a href=\"{}\"><i class=\"fa fa-line-chart\"></i></a></H3></TD>\n".format(t, plot_file))
            else:
                hf.write("<TD ALIGN=CENTER><H3>&nbsp;</H3></TD>\n")

        hf.write("</TR>\n")

    # loop over all the test runs
    for dir in valid_dirs:

        # first look to see if there are any valid tests at all --
        # otherwise we don't do anything for this date
        valid = 0
        for test in all_tests:
            status_file = "%s/%s/%s.status" % (suite.webTopDir, dir, test)
            if os.path.isfile(status_file):
                valid = 1
                break

        if not valid: continue

        # write out the directory (date)
        hf.write("<TR><TD class='date'><SPAN CLASS='nobreak'><A class='main' HREF=\"%s/index.html\">%s&nbsp;</A></SPAN></TD>\n" %
                 (dir, dir) )

        for test in all_tests:

            # look to see if the current test was part of this suite run
            status_file = "%s/%s/%s.status" % (suite.webTopDir, dir, test)
            status = 0

            if os.path.isfile(status_file):

                with open(status_file, 'r') as sf:

                    # status = -1 (failed); 1 (passed); 10 (benchmark update)
                    status = -1
                    for line in sf:
                        if line.find("PASSED") >= 0:
                            status = 1
                            break
                        elif line.find("FAILED") >= 0:
                            status = -1
                            break
                        elif line.find("benchmarks updated") >= 0:
                            status = 10
                            break

            # write out this test's status
            if status == 1:
                hf.write("<TD ALIGN=CENTER title=\"%s\" class=\"passed\"><H3><a href=\"%s/%s.html\" class=\"passed\">:)</a></H3></TD>\n" % (test, dir, test))
            elif status == -1:
                hf.write("<TD ALIGN=CENTER title=\"%s\" class=\"failed\"><H3><a href=\"%s/%s.html\" class=\"failed\">&nbsp;!&nbsp;</a></H3></TD>\n" % (test, dir, test))
            elif status == 10:
                hf.write("<TD ALIGN=CENTER title=\"%s\" class=\"benchmade\"><H3>U</H3></TD>\n" % (test))
            else:
                hf.write("<TD>&nbsp;</TD>\n")


        hf.write("</TR>\n\n")

    hf.write("</TABLE>\n")

    # close
    hf.write("</BODY>\n")
    hf.write("</HTML>\n")

    hf.close()


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# m a i n
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
if __name__== "__main__":
    test_suite(sys.argv[1:])
