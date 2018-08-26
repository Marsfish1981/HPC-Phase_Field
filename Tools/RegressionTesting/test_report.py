import os


CSS_CONTENTS = \
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

HTML_HEADER = \
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

MAIN_HEADER = \
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

def create_css(table_height=16):
    """ write the css file for the webpages """

    css = CSS_CONTENTS.replace("@TABLEHEIGHT@", "{}em".format(table_height))

    with open("tests.css", 'w') as cf:
        cf.write(css)

class HTMLList(object):
    """ a simple class for managing nested HTML lists """

    def __init__(self, of=None):
        # items will hold tuples: (indent, string), where indent
        # specifies how deeply nested we are
        self.list_items = []
        self.current_indent = 0
        self.of = of

    def item(self, content):
        # add an item to the list
        self.list_items.append((self.current_indent, content))

    def indent(self):
        # indent (nest a new list)
        self.current_indent += 1

    def outdent(self):
        # close the current nest level
        self.current_indent -= 1

    def write_list(self):
        # output the list to the outfile, of, specified at creation
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
        for n in range(0, current_indent):
            self.of.write("</ul></li>\n")

        self.of.write("</ul>\n")

class HTMLTable(object):
    """ a simple class for creating an HTML table """

    def __init__(self, out_file, columns=1, divs=None):
        """ create the table object.  Here divs is the name of 
            any HTML div(s) we want to wrap the table with """

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
        """ write the table header """
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
            if any(sstr in line for sstr in ["SUCCESS", "is up to date.",
                                             "Nothing to be done"]):
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

    new_head = HTML_HEADER

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
            root_file = os.path.basename(afile)
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

    new_head = HTML_HEADER + r"""<CENTER><H1><A HREF="../">@TESTDIR@</A> / @TESTNAME@</H1></CENTER>"""

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

    table_height = min(max(suite.lenTestName, 4), 16)

    os.chdir(suite.webTopDir)

    create_css(table_height=table_height)

    valid_dirs, all_tests = suite.get_run_history(active_test_list)

    if suite.do_timings_plots: suite.make_timing_plots(active_test_list)

    #--------------------------------------------------------------------------
    # generate the HTML
    #--------------------------------------------------------------------------
    htmlFile = "index.html"

    title = "%s regression tests" % (suite.suiteName)

    hf = open(htmlFile, 'w')

    header = MAIN_HEADER.replace("@TITLE@", title).replace("@SUBTITLE@", suite.sub_title)

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


    if suite.do_timings_plots:
        hf.write("<tr><td class='date'>plots</td>")
        for t in all_tests:
            plot_file = "{}-timings.png".format(t)
            if os.path.isfile(plot_file):
                hf.write("<TD ALIGN=CENTER title=\"{} timings plot\"><H3><a href=\"{}\"><i class=\"fa fa-line-chart\"></i></a></H3></TD>\n".format(t, plot_file))
            else:
                hf.write("<TD ALIGN=CENTER><H3>&nbsp;</H3></TD>\n")

        hf.write("</TR>\n")

    # loop over all the test runs
    for tdir in valid_dirs:

        # first look to see if there are any valid tests at all --
        # otherwise we don't do anything for this date
        valid = 0
        for test in all_tests:
            status_file = "%s/%s/%s.status" % (suite.webTopDir, tdir, test)
            if os.path.isfile(status_file):
                valid = 1
                break

        if not valid: continue

        # write out the directory (date)
        hf.write("<TR><TD class='date'><SPAN CLASS='nobreak'><A class='main' HREF=\"%s/index.html\">%s&nbsp;</A></SPAN></TD>\n" %
                 (tdir, tdir) )

        for test in all_tests:

            # look to see if the current test was part of this suite run
            status_file = "%s/%s/%s.status" % (suite.webTopDir, tdir, test)
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
                hf.write("<TD ALIGN=CENTER title=\"%s\" class=\"passed\"><H3><a href=\"%s/%s.html\" class=\"passed\">:)</a></H3></TD>\n" % (test, tdir, test))
            elif status == -1:
                hf.write("<TD ALIGN=CENTER title=\"%s\" class=\"failed\"><H3><a href=\"%s/%s.html\" class=\"failed\">&nbsp;!&nbsp;</a></H3></TD>\n" % (test, tdir, test))
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
