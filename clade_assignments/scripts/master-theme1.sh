#!/bin/bash

sh align_to_outgroup.sh
#Rscript makeNodeLabel.R #for adding node labels
python3 labelRenamer.py
sh infer.sh
sh assign_clades.sh

