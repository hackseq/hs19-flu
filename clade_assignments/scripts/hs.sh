#!/bin/bash

python3 -m augur align --sequences hackseq/fas/alignB.fa --fill-gaps --nthreads auto --reference-sequence hackseq/ref/reference_yam_ha.gb --output hackseq/ali/ALIGNED-B.fa
python3 -m augur align --sequences hackseq/fas/alignBNA.fa --fill-gaps --nthreads auto --reference-sequence hackseq/ref/reference_h1n1pdm_ha.gb --output hackseq/ali/ALIGNED-BNA.fa
python3 -m augur align --sequences hackseq/fas/alignH1N1.fa --fill-gaps --nthreads auto --reference-sequence hackseq/ref/reference_h1n1pdm_ha.gb --output hackseq/ali/ALIGNED-H1N1.fa
python3 -m augur align --sequences hackseq/fas/alignH1N1NA.fa --fill-gaps --nthreads auto --reference-sequence hackseq/ref/reference_h1n1pdm_ha.gb --output hackseq/ali/ALIGNED-H1N1NA.fa
python3 -m augur align --sequences hackseq/fas/alignH3N2.fa --fill-gaps --nthreads auto --reference-sequence hackseq/ref/reference_h3n2_ha.gb --output hackseq/ali/ALIGNED-H3N2.fa
python3 -m augur align --sequences hackseq/fas/alignH3N2NA.fa --fill-gaps --nthreads auto --reference-sequence hackseq/ref/reference_h3n2_ha.gb --output hackseq/ali/ALIGNED-H3N2NA.fa
python3 -m augur align --sequences hackseq/fas/PH3N2.fa --fill-gaps --nthreads auto --reference-sequence hackseq/ref/reference_h3n2_ha.gb --output hackseq/ali/ALIGNED-PH3N2.fa
