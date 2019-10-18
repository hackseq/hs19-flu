python3 -m augur refine --tree ../trees/flutreeH3N2_HA.tree --output-tree ../../infer/flutreeH3N2_HA.nwk

python3 -m augur ancestral --tree ../trees/flutreeH3N2_HA.nwk --alignment ../../ali/ALIGNED-PH3N2.fa --output ../../infer/INFERRED-PH3N2.json --inference joint 
