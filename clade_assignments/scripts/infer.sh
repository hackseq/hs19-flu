#need to run these lines beforehand
#python3 -m augur refine --tree ../trees/flutreeB_HA.nwk --output-tree ../trees/flutreeB_HA_2.nwk
#python3 -m augur refine --tree ../trees/flutreeH1N1_HA.nwk --output-tree ../trees/flutreeH1N1_HA_2.nwk
#python3 -m augur refine --tree ../trees/flutreeH3N2_HA.nwk --output-tree ../trees/flutreeH3N2_HA_2.nwk


python3 -m augur ancestral --tree ../trees/flutreeB_HA_2.nwk --alignment ../alignments/ALIGNED-B_yam_clean_labels.fa --output ../../infer/INFERRED-B_yam.json --inference joint 

python3 -m augur ancestral --tree ../trees/flutreeB_HA_2.nwk --alignment ../alignments/ALIGNED-B_vic_clean_labels.fa --output ../../infer/INFERRED-B_vic.json --inference joint 

python3 -m augur ancestral --tree ../trees/flutreeH1N1_HA_2.nwk --alignment ../alignments/ALIGNED-H1N1_clean_labels.fa --output ../../infer/INFERRED-H1N1.json --inference joint 

python3 -m augur ancestral --tree ../trees/flutreeH3N2_HA_2.nwk --alignment ../alignments/ALIGNED-H3N2_clean_labels.fa --output ../../infer/INFERRED-H3N2.json --inference joint


