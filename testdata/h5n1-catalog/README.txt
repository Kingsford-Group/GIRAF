Small test case for finding reassortment architectures.

catalog.correct - the expected output
alignments/*.nex - the input fasta files for each segment

To run,

    cd alignments
    run_mrbayes.sh *.nex
    giraf in.giraf

The "mb" MrBayes executable must be in your path. This will produce
SEG1-SEG2_report files for every pair of segments, and also will produce a
file called "catalog" that contains predictions for the architectures.

the run_mrbayes.sh command may take an hour or two to run. (It has to estimate
trees for 8 segments).

The "giraf" command should take a few minutes at most.
