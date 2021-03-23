#!/bin/bash

while read gt
do 
#cp Newicks/$gt* FixedNewicks/
cp LinkerLengthsKCO/$gt* LinkerLengthsKCORerun/
cp DomainLengthsKCO/$gt* DomainLengthsKCORerun/
done<"Newicks.with.100000.0000001.txt"


