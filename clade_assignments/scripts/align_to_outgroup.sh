#!/bin/bash

python3 -m augur align --sequences ../../Data/B.fa --fill-gaps --nthreads auto --reference-sequence ../references/reference_yam_ha.gb --output ../alignments/ALIGNED-B_yam.fa
python3 -m augur align --sequences ../../Data/B.fa --fill-gaps --nthreads auto --reference-sequence ../references/reference_vic_ha.gb --output ../alignments/ALIGNED-B_vic.fa
python3 -m augur align --sequences ../../Data/H1N1.fa --fill-gaps --nthreads auto --reference-sequence ../references/reference_h1n1pdm_ha.gb --output ../alignments/ALIGNED-H1N1.fa
python3 -m augur align --sequences ../../Data/H3N2.fa --fill-gaps --nthreads auto --reference-sequence ../references/reference_h3n2_ha.gb --output ../alignments/ALIGNED-H3N2.fa
