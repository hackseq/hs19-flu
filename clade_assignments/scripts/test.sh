python3 -m augur refine --tree ../trees/flutreeB_HA.nwk --output-tree ../trees/flutreeB_HA_2.nwk
python3 -m augur ancestral --tree ../trees/flutreeB_HA_2.nwk --alignment ../../ali/ALIGNED-B_clean_labels.fa --output ../../infer/INFERRED-B.json --inference joint 
