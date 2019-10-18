#!/bin/bash

python3 -m augur align --sequences hackseq/fas/alignB.fa --nthreads 20 --reference-sequence hackseq/ref/reference_yam_ha.gb --output hackseq/ali/ALIGNED-B.fa
python3 -m augur align --sequences hackseq/fas/alignBNA.fa --nthreads 20 --reference-sequence hackseq/ref/reference_h1n1pdm_ha.gb --output hackseq/ali/ALIGNED-BNA.fa
python3 -m augur align --sequences hackseq/fas/alignH1N1.fa --nthreads 20 --reference-sequence hackseq/ref/reference_h1n1pdm_ha.gb --output hackseq/ali/ALIGNED-H1N1.fa
python3 -m augur align --sequences hackseq/fas/alignH1N1NA.fa --nthreads 20 --reference-sequence hackseq/ref/reference_h1n1pdm_ha.gb --output hackseq/ali/ALIGNED-H1N1NA.fa
python3 -m augur align --sequences hackseq/fas/alignH3N2.fa --nthreads 20 --reference-sequence hackseq/ref/reference_h3n2_ha.gb --output hackseq/ali/ALIGNED-H3N2.fa
python3 -m augur align --sequences hackseq/fas/alignH3N2NA.fa --nthreads 20 --reference-sequence hackseq/ref/reference_h3n2_ha.gb --output hackseq/ali/ALIGNED-H3N2NA.fa
python3 -m augur align --sequences hackseq/fas/PH3N2.fa --nthreads 20 --reference-sequence hackseq/ref/reference_h3n2_ha.gb --output hackseq/ali/ALIGNED-PH3N2.fa
