fname = 'aligned.fasta'
newfname = 'aligned_correct.fasta'

with open(fname, 'r') as f:
	with open(newfname, 'w') as fout:
		for line in f:
			sline = line.strip()
			if line.startswith('>'):
				dat = sline.split(',')
				fout.write(dat[0]+'\n')
			else:
				fout.write(sline + '\n')
