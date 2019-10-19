import os

inputPath = os.fsencode(".")
outputPath = os.fsencode("../clade_assignments/alignments")
refPath = os.fsencode("../clade_assignments/references")

for file in os.listdir(inputPath):
    filename = os.fsdecode(file)
    if filename.endswith(".fa"):

    else:
        continue
