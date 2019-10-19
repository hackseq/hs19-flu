#!/usr/bin/bash

python3 ./identifyClades.py ../alignments/ALIGNED-B.fa ../clade_defs/clades_yam_ha.tsv --reference ../references/reference_yam_ha.gb --output ../assignments/B_yam.tsv

python3 ./identifyClades.py ../alignments/ALIGNED-BNA.fa ../clade_defs/clades_h1n1pdm_ha.tsv --reference ../references/reference_h1n1pdm_ha.gb --output ../assignments/BNA_h1n1pdm.tsv

python3 ./identifyClades.py ../alignments/ALIGNED-H1N1.fa ../clade_defs/clades_h1n1pdm_ha.tsv --reference ../references/reference_h1n1pdm_ha.gb --output ../assignments/H1N1_h1n1pdm.tsv

python3 ./identifyClades.py ../alignments/ALIGNED-H1N1NA.fa ../clade_defs/clades_h1n1pdm_ha.tsv --reference ../references/reference_h1n1pdm_ha.gb --output ../assignments/H1N1NA_h1n1pdm.tsv

python3 ./identifyClades.py ../alignments/ALIGNED-H3N2.fa ../clade_defs/clades_h3n2_ha.tsv --reference ../references/reference_h3n2_ha.gb --output ../assignments/H3N2_h3n2.tsv

python3 ./identifyClades.py ../alignments/ALIGNED-H3N2NA.fa ../clade_defs/clades_h3n2_ha.tsv --reference ../references/reference_h3n2_ha.gb --output ../assignments/H3N2NA_h3n2.tsv

python3 ./identifyClades.py ../alignments/ALIGNED-PH3N2.fa ../clade_defs/clades_h3n2_ha.tsv --reference ../references/reference_h3n2_ha.gb --output ../assignments/PH3N2_h3n2.tsv
