"""
append outgroup
Take a sequence from raw alignment, make a query, run BLASTn against a 
given database and append sequence to the raw alignment
"""

evalue = "1e-10"
blasttask = "blastn"
outfmt = "6"

from os import system, remove
from Bio import SeqIO

def extract_seq(multal_file):
	name = ""
	seq = ""
	for q in open(multal_file):
		q = q.strip()
		if q[0] == '>' and name != "":
			break
		elif q[0] == '>':
			name = q[1:]
		else:
			seq += q
	temp_query = open("temp_query.fa", "w")
	temp_query.write('>'+name+'\n')
	temp_query.write(seq+'\n')
	temp_query.close()

def blast_search(db_name, multal_file):
	system("blastn -task "+blasttask+" -query temp_query.fa -out temp.blastout -outfmt "+outfmt+" -evalue "+evalue+" -db "+db_name)
	for q in open("temp.blastout"):
		try:
			q = q.split()
			_id = q[1]
		except:
			return 0
		handle = open(db_name)
		for record in SeqIO.parse(handle, "fasta"):
			if record.id == _id or record.name == _id or record.description == _id:
				alignment = open(multal_file, 'a')
				SeqIO.write(record, alignment, "fasta")
				return 1
	return 0

def append_outgroup(db_name, multal_file):
	extract_seq(multal_file)
	bs = blast_search(db_name, multal_file)
	system("rm temp_query.fa temp.blastout")
	return bs

#append_outgroup("outgroups/prochlorococcus-synechococcus/GCF_001632105.1_ASM163210v1_cds_from_genomic.fna","Prochlorococcus_marinus_sampled/1_results/aligns/1.fa.raw")

