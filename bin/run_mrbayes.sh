#!/bin/sh

# make sure the given command exists
check_exists () 
{
    if which $1 > /dev/null ; then
        return 0
    fi
    echo "error: GIRAF requires $1 ($2) which could not be found." > /dev/stderr
    return 3
}

check_exists "mb" "MrBayes http://mrbayes.csit.fsu.edu/"

if [ ${#} -lt 1 ] ; then 
    echo "Usage: run_mrbayes alignment1.nex [alignment2.nex...]" >> /dev/stderr
    echo "The mb command must be in your PATH. It will be run on each" >>/dev/stderr
    echo "alignment, producing trees and also in.giraf to supply to giraf." >/dev/stderr
    exit 3
fi

rm -f in.giraf

for nex in $* ; do

    # run mrbayes
    cat <<EOT | mb
execute $nex;
lset nst=6 rates=invgamma;
mcmc ngen=200000 samplefreq=200;
sump burnin=500;
sumt burnin=500;
quit;
EOT

    # for some reason, we have to run it twice to really get the
    # output
    cat <<EOT | mb
execute $nex;
sump burnin=500;
sumt burnin=500;
quit;
EOT

bn=`basename $nex .nex`
echo "$bn $nex.run1.t $nex.run2.t" >> in.giraf

done

