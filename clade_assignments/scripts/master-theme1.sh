#!/bin/bash

python3 automagic.py
#Rscript makeNodeLabel.R #for adding node labels
py labelRenamer.py
sh infer.
sh assign_clades.sh

