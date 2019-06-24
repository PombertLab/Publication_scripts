#!/bin/bash

n=$(awk '/DiscRep_ALL:DISC_PARTIAL_PROBLEMS/{print NR}' discrep)

line=$(head -$n discrep | tail -1)

regex="::([0-9]+)"
if [[ $line =~ $regex ]]; then
	echo "Partial proteins remaining: ${BASH_REMATCH[1]}"
else
	echo "couldn't read discrep, exiting..."
	exit
fi

# Begin processing loop
while :
do

# Confirm processing
read -p "continue (y/n)? " -n 1 cont
echo

if [ "$cont" == "y" ]; then
	echo "Continuing!"
else
	echo "exiting..."
	exit
fi

# Move to next contig
((n=n+1))

# Open contig for editing
previous=""
line=$(head -$n discrep | tail -1)
regex="contig-([0-9]+):"
if [[ $line =~ $regex ]]; then
	if [ "$line" == "$previous" ]; then
		continue
	fi
	previous=${BASH_REMATCH[1]}
	echo "Opening contig-${BASH_REMATCH[1]}"
	$(art GFF3/contig-${BASH_REMATCH[1]}.embl)
else
	echo "couldn't read contig from discrep, exiting..."
	exit
fi

done
# End processing loop
