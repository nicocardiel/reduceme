#!/bin/bash
if [ $1 == "" ]
then
  echo "#> ERROR: insuficient parameters in prompt line running exheads.sh"
  exit
fi

# installation prefix of the REDUCEME programs (--program-prefix), taken
# from the name of this script (e.g. R6-exheads.sh -> R6-)
progprefix=$(basename "$0" | sed 's/exheads\.sh$//')

echo "> exheads $1"
"${progprefix}exheads"<<endexheads
$1
endexheads
