#!bin/bash

grep "^>" Data/alignH3N2NA.fa | sed 's@ @_@g ; s@/@ @' > alignH3N2NA.fa.header
grep "^>" Data/alignH3N2.fa | sed 's@ @_@g ; s@/@ @'> alignH3N2.fa.header

awk -F" " '(NR==FNR){ h[$2]=$1;next; }( $2 in h){print h[$2],$1}' alignH3N2NA.fa.header alignH3N2.fa.header > HA_NA_map.txt.txt
