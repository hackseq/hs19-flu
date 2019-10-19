#!/bin/bash

python3 align_to_outgroup.py
#Rscript makeNodeLabel.R #for adding node labels
python3 labelRenamer.py
sh infer.sh
sh assign_clades.sh

