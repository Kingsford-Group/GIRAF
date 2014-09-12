GIRAF
=====

GIRAF is a command-line utility that will read collections of trees to attempt
to find influenza reassortments.

Contents
1. INSTALLATION:
2. TESTING INSTALLATION
3. TYPICAL USAGE
4. GIRAF INPUT FILES
5. OPTIONS:
6. COMPILING, IF NECESSARY
7. ADVANCED MODE


1) INSTALLATION:
----------------

The heart of GIRAF is a command line program. Versions of the program for Apple
Mac, Windows, and Linux are in the bin/ directory of the distributed .zip file.
(The Linux version was built on RedHat Enterprise Linux 5 --- it may not work
with other linux versions; if it doesn't work on your machine, trying compiling
it.) 


2) TESTING INSTALLATION
-----------------------

You can test your installation by running the following commands (of course,
"giraf" should be replaced by perhaps giraf_mac or giraf_windows, and should
either be in your PATH or you should give the full path to it):

    cd testdata/h5n1
    giraf h5n1.giraf

This will run a small example (should take ~ a minute). You can compare the
output in "left-right_report" with "report.correct". They should basically
match.  (The number following "ID:" may be different.)

You can run a slightly larger test (that requires the MrBayes command "mb" be
installed on your path) with the following commands:

    cd testdata/h5n1-catalog/alignments
    run_mrbayes.sh *-aln.nex
    giraf in.giraf

You can then compare the file "catalog" produced by the above commands to the
file "catalog.correct" in the testdata/h5n1-catalog directory.


3) TYPICAL USAGE
----------------

GIRAF requires some way to produce .nex files that contain many trees. The
easiest way is to use MrBayes. You should install MrBayes and put the "mb"
executable file someplace in your PATH.

GIRAF comes with a script to run MrBayes to produce the required trees: 

    bin/run_mrbayes.sh

This script works only on UNIX, or Mac. On Windows, you have to run MrBayes by
hand.  NOTE: if you produce your own .nex files without using the
run_mrbayes.sh script, you should be aware that giraf will ignore the first 500
trees in each file by default. This is to give the MCMC process used by MrBayes
as chance to "burn-in". If you want to change the number of trees ignored by
MrBayes use the --burnin=N option to GIRAF. See OPTIONS below.

To find reassortments among the HA, NA, and MP segments:

    run_mrbayes.sh HA.nex NA.nex MP.nex
    giraf in.giraf

The first command will run MrBayes on each of the alignment NEXUS files
provided. It will produce many *.run1.t and *.run2.t files containing trees.
GIRAF will read these files.  The run_mrbayes.sh command also produces an
"in.giraf" file that serves can be passed to giraf (see GIRAF FILES below).

This will produce:

    SEG1-SEG2_report  - reassortments between SEG1 and SEG2
    catalog           - reassortment architectures

It will also produce possibly many intermediate files with names ending with
_splits, _trees, _graph.labels, _graph.labelled, _bicliques. These are files
that GIRAF writes and reads as it is working. In most cases, they can be safely
ignored or deleted after GIRAF finishes.

If you are only considering 2 segments, GIRAF will not produce the "catalog"
file (since "architectures" do not make sense for only 2 segments).

While running GIRAF, status messages will be output to standard output (the
screen).  The messages are preceded by a "word:" that indicates the "phase" of
giraf that is currently happening (mcmc_split_info = parsing trees;
build_incompat_graph = finding incompatible splits; extract_reassortments =
analyze the splits for reassortments; the names correspond to the programs in
advanced mode (see below)).


4) GIRAF INPUT FILES
--------------------

In typical usage, the giraf command is "giraf in.giraf". Here, "in.giraf" is
the name of a file that contains the information about the segments. The format
for this file is:

    name1 treefile1.nex [treefile2.nex...]
    name2 treefile1.nex [treefile2.nex...]

where "name1" and "name2" are names you want to give to the segments (they can
be anything you want without spaces). The nexus files contain estimates for the
trees for that segment.  The .nex files don't have to have a .nex extension,
and the ones produced by MrBayes will end with .t. 

The run_mrybayes.sh script will produce a file called "in.giraf" of the
required format. An example file is:

    HA left.nex.run1.t left.nex.run2.t
    NA right.nex.run1.t right.nex.run2.t


5) OPTIONS:
-----------

Typical Options:

These are options you might actually want to change:

   --burnin=N   (default N=500)
        Ignore the first N trees from each input tree file.

   --cull=F  (default F=0.05)
        Remove from considerations splits that happen in fewer than F fraction
        of the trees.

   --threshold=F (default=0.70) : confidence cutoff
        The confidence treshold for reporting a reassortment.  Higher values
        mean GIRAF will be more strict when outputing reassortments.

Advanced Options:

These are options that fundementally change how GIRAF works. Almost certainly
you do not want to change them. Understanding their meanings requires knowing
how GIRAF works in some detail.

   --use-dist=[0,1]   (default 1)
        If 0, distances in the trees will not be used to filter candidate sets
        of reassortments. This may produce more putative reassortments.

   --test-all-candidates=[0,1] (default 0)
        If 1, even large candidate sets will be tested as potential reassortant
        sets.

   --non-star=[0,1]   (default 1)
        Restrict what types of bicliques are considered. By default, only
        non-star bicliques (biclique (U,V) where both U and V contain at least
        2 splits) are considered. If non-star=0, all bicliques count.

   --single
        Don't require that a candidate be supported by multiple lines of
        biclique evidence. This may produce more putative reassortments.

Debugging options:

Unless you are trying to extend or debug GIRAF itself, you probably do not need
to use these options.

   --debug-out-pairs  : output result of statistical tests (debugging only)
   --debug-out-unlabeled : output unlabeled incompat graph (debugging only)

   --version-0.9-compat
        Try to be as similar to version 0.9 of GIRAF as possible. Version 0.9
        was primarily written in PERL, and the are some small differences in
        how it sorts and filters putative reassortments. With this option, an
        effort is made to match the version 0.9 behavior as much as possible.
        You should not need to use this option unless you are trying to
        compare with results that were generated using version 0.9 of GIRAF.

        All new work should be done with GIRAF 1.0 or later.

6) COMPILING, IF NECESSARY
--------------------------

If you change GIRAF, or want to get it to work on a different type of machine,
GIRAF must be compiled, and probably requires the GNU g++ compiler.  To compile
it, type the commands:

    cd src
    make

You should see number of lines starting with "g++" appear, along with with a
few others. There should be no warnings or errors.

This will produce a program file called "giraf" in the src/ directory. You can
move the giraf file anyplace you want.

You can use "make clean" to reset the compilation process.

7) ADVANCED MODE

You probably do not need to use advanced mode. The usage described in (2)
TYPICAL USAGE above is easier, faster, and in some ways more powerful.

The advanced mode of GIRAF is intend to make running multiple commands in
parallel easier. It can also be used if you want to replace one part of GIRAF
with your own analysis.  The "advanced mode" is almost never what you want, and
the standard "giraf" command is probably sufficient for nearly all users.  At
present, the advance mode cannot compute reassortment architectures among > 2
segments.

To create the advanced version of GIRAF, run the command:

    make advanced

The advanced version of GIRAF is actually 3 separate commands. The reason for
this is to help run multiple analyses in parallel and to restart in case a step
fails. Step 1 is to use MrBayes to create a collection of trees. Then you would
execute something like the following commands:

	mcmc_split_info left left.nex.run*.t
	mcmc_split_info right right.nex.run*.t
	build_incompat_graph left right LR
	extract_reassortments left right LR final

Details about the above commands are given below.

1. The command:

	mcmc_split_info right right.nex.run*.t

will read the output of MrBayes and produce "right_splits", "right_dist", and
"right_trees" just exactly as parse_mcmc.pl would. The _splits file will
contain the splits that are found in the trees, and _dist will contain the
distances between isolates.

This command takes somewhere between 1 and 3 minutes to run on 137 genomes on
my macbook. 

2. The command:

	build_incompat_graph left right OUT

will read the _splits & _dist files (with prefixes "left" and "right") and create: 

	OUT_graph.labelled
	OUT_graph.labels

This command takes about 1 minute to run on 137 genomes on my macbook.

3. The command:

	extract_reassortments left right left-right finalout

It will produce:

	finalout_bicliques
	finalout_report

The _report file contains the putative reassortant sets.

Each of the commands can take options, which you can find by entering the
command with no parameters. The options behave as described in section (4)
OPTIONS above.


