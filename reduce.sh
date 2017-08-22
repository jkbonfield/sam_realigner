#!/bin/sh

REDUCTION=${REDUCTION:-.7}

# Repeatedly run command, reducing file input until it no longer
# crashes.  File input is assumed to be the second argument.
# file no longer crashes.

c=$1;shift
in=$1;shift
opts=${@+"$@"}

seed=0
for j in `seq 1 100`
do
    echo "Loop $j; `egrep -cv '^@' $in` lines"
    for i in `seq 1 100`
    do
	seed=`expr $seed + 1`
	echo "    trying $i / $seed"
	samtools view -h -s $seed${REDUCTION} -o _tmp $in
	eval $c _tmp $opts > /dev/null 2>&1
	if [ $? != 0 ];then break; fi
    done
    if [ $i == 100 ]; then break; fi
    cp _tmp _last.sam
    in=_last.sam
done


