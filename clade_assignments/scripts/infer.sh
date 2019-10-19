#need to run these lines beforehand
#python3 -m augur refine --tree ../trees/flutreeB_HA.nwk --output-tree ../trees/flutreeB_HA_2.nwk
#python3 -m augur refine --tree ../trees/flutreeH1N1_HA.nwk --output-tree ../trees/flutreeH1N1_HA_2.nwk
#python3 -m augur refine --tree ../trees/flutreeH3N2_HA.nwk --output-tree ../trees/flutreeH3N2_HA_2.nwk


python3 -m augur ancestral --tree ../trees/flutreeB_HA_2.nwk --alignment ../hackseq/ali/ALIGNED-B_clean_labels.fa --output ../../infer/INFERRED-B.json --inference joint 


#python3 -m augur ancestral --tree ../trees/flutreeB_HA_2.nwk --alignment ../hackseq/ali/ALIGNED-BNA_clean_labels.fa --output ../../infer/INFERRED-BNA.json --inference joint 


#python3 -m augur ancestral --tree ../trees/flutreeH1N1_HA_2.nwk --alignment ../hackseq/ali/ALIGNED-H1N1_clean_labels.fa --output ../../infer/INFERRED-H1N1.json --inference joint 


#python3 -m augur ancestral --tree ../trees/flutreeH1N1_HA_2.nwk --alignment ../hackseq/ali/ALIGNED-H1N1A_clean_labels.fa --output ../../infer/INFERRED-H1N1A.json --inference joint



#python3 -m augur ancestral --tree ../trees/flutreeH3N2_HA_2.nwk --alignment ../hackseq/ali/ALIGNED-H3N2_clean_labels.fa --output ../../infer/INFERRED-H3N2.json --inference joint

#python3 -m augur ancestral --tree ../trees/flutreeH3N2_HA_2.nwk --alignment ../hackseq/ali/ALIGNED-H3N2NA_clean_labels.fa --output ../../infer/INFERRED-H3N2NA.json --inference joint 


#python3 -m augur ancestral --tree ../trees/flutreeH3N2_HA_2.nwk --alignment ../hackseq/ali/ALIGNED-PH3N2_clean_labels.fa --output ../../infer/INFERRED-PH3N2.json --inference joint 

