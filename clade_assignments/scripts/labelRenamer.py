import sys
try:    
    in_file = sys.argv[1]
except IndexError:
    print(f"Clean fasta labels. Usage: {sys.argv[0]} alignment.fa")
    sys.exit(0)

prefix, postfix = in_file.split(".")[0], in_file.split(".")
newfname = f"{prefix}_clean_labels{postfix}"

with open(in_file, 'r') as f:
	with open(newfname, 'w') as fout:
		for line in f:
			sline = line.strip()
			if line.startswith('>'):
				dat = sline.split(',')
				fout.write(dat[0]+'\n')
			else:
				fout.write(sline + '\n')
