import dendropy
import re
import functools
import itertools
import sys
import argparse
import os
from Bio import SeqIO


###########################################
# Command line interface
this_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
default_out = os.path.join(this_dir, "clade_assignments.out")

parser = argparse.ArgumentParser(description='Influenza sequence Clade assignment')
parser.add_argument('alignment', 
        type=str, 
        help='FASTA file containing aligned sequences to assign')
parser.add_argument('clade_criteria', 
        type=str, 
        help='tsv file containing criteria for a sequence to belong to a clade. '
        'Note: clade file is influenza strain specific (eg H3N2)'
        )
parser.add_argument('--reference', 
        type=str, 
        help='Genbank file containing the reference gene')
parser.add_argument('--output',
        type = str,
        help = f'path to desired output. Default: {default_out}. '
        'The output file is a tsv matching sequences names from FASTA to clades',
        default = default_out
        )

args = parser.parse_args()
###########################################

###########################################
# Input/Output files: example

ffasta   = 'aligned.fasta'      # Input: FASTA file correctly aligned (1701 nucleotides) 
fclades  = 'clades_h3n2_ha.tsv' # Input: Criteria to belong to a clade
fnameout = 'FASTA2018-5_clades.tsv' # Output: TSV file with sequence names and its clades 

###########################################
# A/Beijing/32/1992 --> h3h2_ha_outgroup.gb
# CDS             1..1701
#                 /db_xref="GI:857408"
#                 /product="hemagglutinin"
# CDS             1..48
#                 /product="Signal peptide"
#                 /gene="SigPep"
# CDS             49..1035
#                 /product="HA1 protein"
#                 /gene="HA1"
# CDS             1036..1698
#                 /product="HA2 protein"
#                 /gene="HA2"
with open(args.reference, 'r')as gbfile:
	genbank = SeqIO.read(gbfile, 'gb')
	for feature in genbank.features:
		if('host' in feature.qualifiers):
			sequenceSize = int(feature.location.end)
		if('gene' in feature.qualifiers):
			if(feature.qualifiers['gene'][0] == 'HA1'):
				HA1_beg = int(feature.location.start) + 1
				HA1_end = int(feature.location.end)
			if(feature.qualifiers['gene'][0] == 'HA2'):
				HA2_beg = int(feature.location.start)  + 1
				HA2_end = int(feature.location.end)

###########################################
# Read clade constrains 
# Example of clade constrain used to identify clades:
#clade	gene	site	alt
#3b	HA2	158	N
cladesConstrains = {}
with open(args.clade_criteria, 'r') as f:
	for line in f:
		info = (line.strip()).split("\t")

		clade = info[0] 
		gene  = info[1]
		site  = info[2]
		alt   = info[3] # Nucleotide (ACGT) or Amino-Acid 

		if clade != "clade":
			if clade not in cladesConstrains.keys():
				cladesConstrains[clade] = []
			cladesConstrains[clade].append([gene, site, alt])

#print(cladesConstrains)

###########################################
# Parameters used in the amino-acid translation: 
# 1- Table to translate codons to AA, 
# 2- Table to translate an uncertain codon to one or more certain codons.
bases = "tcag".upper()
codons = [a + b + c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

basesComb = functools.reduce(lambda a,b:a+b, [[''.join(a) for a in list(itertools.combinations(bases, i))] for i in range(1,5)]) + ['tcag'.upper()]
basesUncertain = 'tcagywkmsrhbdvxn'.upper()
bases_table = dict(zip(basesUncertain, basesComb))

codonSize = 3

basepos_aa = {}
for i in range(codonSize):
	basepos_aa[i] = {}
	for b in bases:
		basepos_aa[i][b] = set()
		for c in codons:
			if c[i] == b:
				basepos_aa[i][b].add(c)

# Translation works in 3 steps: (1) 1 uncertain codon --> (2) 1 or more certain codons --> (3) Each certain codon translates to 1 amino-acid  
def translateToAA(seq):
	seqAA = []
	for i in range(0, len(seq), codonSize):
		uncertainCodon = seq[i:i+codonSize] # Codon is made of 3 (un)certain nucleotides
				
		# Identify amino-acid for a given codon
		certainCodons = None
		for j, nuc_j in enumerate(uncertainCodon):
			# Check all possible codons with nucleotide c_j at position j
			certainCodons_j = set()

			if nuc_j in bases_table.keys(): 
				for b in bases_table[nuc_j]:
					certainCodons_j = certainCodons_j.union(basepos_aa[j][b])
				# Intersection with existent codons
				if certainCodons:
					certainCodons = certainCodons.intersection(certainCodons_j)
				else:
					certainCodons = certainCodons_j
			else:
				print("ERROR Sequence " + name + " has unidentified nucleotide " + nuc_j + " / Codon: " + uncertainCodon)
				print(str(bases_table.keys()))

		possibleAA    = set()
		for cod in certainCodons:
			possibleAA.add(codon_table[cod])
		seqAA.append([uncertainCodon, certainCodons, possibleAA])
	return seqAA

###########################################
# Read sequences and assign clades
cladesInfo  = {}

genomesOneClade = []
genomesMoreThanOneCladeButOneSpecific = []
genomesMoreThanOneCladeMoreThanOneSpecific = []
genomesNoClade = []

seqsWithError = []

with open(args.alignment, 'r') as f:

	p_name = re.compile(">(.+)")
	p_seq  = re.compile("([AaCcGgTtUuMmRrWwSsYyKkVvHhDdBbXxNn]+)")

	name = ""
	seq  = ""
	
	for line in f:
		m = p_name.search(line)
		if m:
			if name:

				# Check sequence size
				if len(seq) == sequenceSize:

					# Translate to amino-acid
					seqAA_HA1 = translateToAA(seq[(HA1_beg-1):HA1_end]) # HA1_beg-1: Python starts at 0 / HA1_end it should include the end position
					seqAA_HA2 = translateToAA(seq[(HA2_beg-1):HA2_end])

					# Check if sequence belongs to one or more clades
					clades = [] 
					for cladeName in cladesConstrains.keys():
						cinf = cladesConstrains[cladeName] # list of constrains

						# Check if sequence satisfies criteria of a specific clade
						satisfiedConstrains   = []
						unsatisfiedConstrains = []
						for constrain in cinf:

							gene  = constrain[0]
							site  = int(constrain[1])
							alt   = constrain[2] # Nucleotide (ACGT) or Amino-Acid 

							# Nucleotide
							if gene == "nuc":
								nuc = seq[site-1] # Since Python starts at 0
								if(nuc == alt):
									satisfiedConstrains.append([gene,site,alt,nuc])
								else:
									unsatisfiedConstrains.append([gene,site,alt,nuc])
							# Amino-Acide
							else:
								aa = None
								if gene == "HA1":
									aa = seqAA_HA1[site-1] # Since Python starts at 0
								elif gene == "HA2":
									aa = seqAA_HA2[site-1] # Since Python starts at 0

								if aa:
									uncertainCodon = aa[0] 
									certainCodons  = aa[1]
									possibleAA     = aa[2]
									if(alt in possibleAA):
										satisfiedConstrains.append([gene,site,alt,aa])
									else:
										unsatisfiedConstrains.append([gene,site,alt,aa])
								else:
									print("ERROR! unidentified gene " + gene)

							#A4	HA1	31	S

						# Check if clade is OK
						if (len(satisfiedConstrains) + len(unsatisfiedConstrains)) == len(cinf):							
							
							if len(satisfiedConstrains) == len(cinf):
								clades.append(cladeName)
								#print(name + " = Clade " + cladeName + " OK ")
							
							# Add info about genome regardless if it belongs to the clade or not
							if cladeName not in cladesInfo.keys():
								cladesInfo[cladeName] = []
							cladesInfo[cladeName].append([name,satisfiedConstrains,unsatisfiedConstrains])


					# Check how many clades the sequence satisfied
					if len(clades) == 1:
						genomesOneClade.append([name, clades])
					elif len(clades) > 1:
						# Identify specific clades
						specificClades = []
						for c1 in clades:
							addClade = True
							removeClade = ""
							for c2 in specificClades:
								if c1.startswith(c2):
									removeClade = c2
									break
								elif c2.startswith(c1):
									addClade = False
									break 
							if addClade:
								specificClades.append(c1)
							if removeClade:
								specificClades.remove(c2)
						
						if len(specificClades) == 1:
							genomesMoreThanOneCladeButOneSpecific.append([name, clades, specificClades])
						elif len(specificClades) > 1:
							genomesMoreThanOneCladeMoreThanOneSpecific.append([name, clades, specificClades])
						else:
							print("ERROR! Problem in the specific clades identification : " + name + " " + str(clades))

					else:
						genomesNoClade.append([name, clades])

				# Sequence not ok in size; ignore it.
				else:
					seqsWithError.append([name, len(seq)])

			name = m.group(1)
			seq  = "" 
		else:
			m = p_seq.search(line)
			if m:
				seq = seq + m.group(1) 


#######################################
# Summary and write files with clades
with open(args.output, 'w') as fout:

	fout.write("name\tnbCladesAll\tnbCladesSpec\tcladesAll\tcladesSpec\n")

	allClades = {}
	print("Sequences that satisfied 1 clade = " + str(len(genomesOneClade)))
	for seqInfo in genomesOneClade:
		name   = seqInfo[0]
		clades = seqInfo[1]

		key = ",".join(sorted(clades))
		if key not in allClades.keys():
			allClades[key] = 0 
		allClades[key] += 1

		fout.write(name + "\t" + str(len(clades)) + "\t" + str(len(clades)) + "\t" + ",".join(clades) + "\t" + ",".join(clades) + "\n")
		#print("\t" + name + " = " + str(clades))

	print("Sequences that satisfied >1 clade BUT only 1 specific clade  = " + str(len(genomesMoreThanOneCladeButOneSpecific)))
	for seqInfo in genomesMoreThanOneCladeButOneSpecific:
		name   = seqInfo[0]
		clades = seqInfo[1]
		specificClades = seqInfo[2]

		key = ",".join(sorted(specificClades))
		if key not in allClades.keys():
			allClades[key] = 0  
		allClades[key] += 1	

		fout.write(name + "\t" + str(len(clades)) + "\t" + str(len(specificClades)) + "\t" + ",".join(clades) + "\t" + ",".join(specificClades) + "\n")
		#print("\t" + name + " = " + str(clades) + " / " + str(specificClades))


	print("Sequences that satisfied >1 clade AND >1 specific clade  = " + str(len(genomesMoreThanOneCladeMoreThanOneSpecific)))
	for seqInfo in genomesMoreThanOneCladeMoreThanOneSpecific:
		name   = seqInfo[0]
		clades = seqInfo[1]
		specificClades = seqInfo[2]
		#print("\t" + name + " = " + str(clades) + " / " + str(specificClades))

		key = ",".join(sorted(specificClades))
		if key not in allClades.keys():
			allClades[key] = 0
		allClades[key] += 1	

		fout.write(name + "\t" + str(len(clades)) + "\t" + str(len(specificClades)) + "\t" + ",".join(clades) + "\t" + ",".join(specificClades) + "\n")
		#print("\t" + name + " = " + str(specificClades))

	print("Sequences that satisfied NO clade = " + str(len(genomesNoClade)))
	allClades["Unassigned"] = len(genomesNoClade)
	for seqInfo in genomesNoClade:
		name   = seqInfo[0]
		clades = seqInfo[1]
		fout.write(name + "\t" + "0" + "\t" + "0" + "\t" + "Unassigned" + "\t" + "Unassigned" + "\n")

	print("Sequences with Errors: " + str(seqsWithError))
	for seqInfo in seqsWithError:
		name    = seqInfo[0]
		seqSize = seqInfo[1]		
		fout.write(name + "\t" + "0" + "\t" + "0" + "\t" + "Error" + "\t" + "Error" + "\n")

	print("Clades identified")
	allKeys = sorted(allClades.keys(), key=lambda x: (x.count(","),x))
	for cladeName in allKeys:
		print(cladeName + " = " + str(allClades[cladeName]))

