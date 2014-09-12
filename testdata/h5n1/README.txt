A test case to check your installation of GIRAF.

left.nex - HA segments from 35 H5N1 isolates
right.nex - NA segments from 35 H5N1 isolates
names.txt - mapping between taxa names, genbank identifiers, and isolate name
*.run*.t - Trees estimated with MrBayes
h5n1.giraf - the graph input file

To Run:

    giraf h5n1.giraf

This will produce HA-NA_report which you can compare with the expected output
in "report.correct".


If you want to re-run MrBayes, you can use the commands:

    run_mrbayes.sh left.nex right.nex
    giraf in.giraf

This will produce the file left-right_report can can be compared with
report.correct.

Note: because the MrBayes sampling process is randomized, there is some small
probability that you will not get the correct output. If several repetitions
of the command do not produce the expected output, something may be wrong.


