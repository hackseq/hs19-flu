#!/usr/bin/bash

python3 ./identifyClades.py ../alignments/ALIGNED-B_yam_clean.fa ../clade_defs/clades_yam_ha.tsv --reference ../references/reference_yam_ha.gb --output ../assignments/B_yam.tsv

python3 ./identifyClades.py ../alignments/ALIGNED-B_vic_clean.fa ../clade_defs/clades_vic_ha.tsv --reference ../references/reference_vic_ha.gb --output ../assignments/B_vic.tsv

python3 ./identifyClades.py ../alignments/ALIGNED-H1N1_clean.fa ../clade_defs/clades_h1n1pdm_ha.tsv --reference ../references/reference_h1n1pdm_ha.gb --output ../assignments/H1N1_h1n1pdm.tsv

python3 ./identifyClades.py ../alignments/ALIGNED-H3N2_clean.fa ../clade_defs/clades_h3n2_ha.tsv --reference ../references/reference_h3n2_ha.gb --output ../assignments/H3N2_h3n2.tsv