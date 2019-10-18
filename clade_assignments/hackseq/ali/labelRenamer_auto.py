import os

directory = os.fsencode(".")

for file in os.listdir(directory):
     filename = os.fsdecode(file)
     if filename.endswith(".fa"):
         in_file = filename;
         prefix, postfix = in_file.split(".")[0], in_file.split(".")
         newfname = f"{prefix}_clean_labels{postfix}"

         with open(in_file, 'r') as f:
             with open(newfname, 'w') as fout:
                 for line in f:
                     sline = line.strip()
                     if line.startswith('>'):
                         dat = sline.split(',')
                         fout.write(dat[0] + '\n')
                     else:
                         fout.write(sline + '\n')
     else:
         continue

