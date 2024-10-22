import os

directory = os.fsencode("../alignments")

for file in os.listdir(directory):
     filename = os.fsdecode(file)
     if filename.endswith(".fa"):
         in_file = filename
         prefix, postfix = in_file.split(".")[0], in_file.split(".")[1]
         newfname = f"{prefix}_clean.{postfix}"

         with open("../alignments/" + in_file, 'r') as f:
             with open("../alignments/" + newfname, 'w') as fout:
                 for line in f:
                     sline = line.strip()
                     if line.startswith('>'):
                         dat = sline.split('_')
                         fout.write(dat[0] + '\n')
                     else:
                         fout.write(sline + '\n')
     else:
         continue

