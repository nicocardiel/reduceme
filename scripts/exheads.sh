#!/bin/bash
if [ $1 == "" ]
then
  echo "#> ERROR: insuficient parameters in prompt line running exheads.sh"
  exit
fi

echo "> exheads $1"
R5-exheads<<endexheads
$1
endexheads
