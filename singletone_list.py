"""
singletone_list_construct
input:
1. Proteinortho_file
2. Fasta_dir
"""

from optparse import OptionParser

def get_fasta_names(proteinortho_file):
	for q in open(proteinortho_file):
		q = q.strip()
		if q[0] == '#':
			q = q.split('\t')
			return q[3:]


def get_ogg_genes(proteinortho_file):
	genedict = dict()
	for q in open(proteinortho_file):
		q = q.strip().split('\t')
		for p1 in q[3:]:
			for p2 in p1.split(','):
				genedict[p2] = 0
	return genedict


def getsingletones(proteinortho_file, fasta_dir, outfile_name):
	fasta_names = get_fasta_names(proteinortho_file)
	genedict_ogg = get_ogg_genes(proteinortho_file)
	outfile = open(outfile_name, "w")
	for fasta in fasta_names:
		for string in open(fasta_dir+fasta):
			if string[0] == '>':
				_id = string.strip()[1:]
				try:
					k = genedict_ogg[_id]
				except:
					outfile.write(_id + '\n')
	outfile.close()


parser = OptionParser()
parser.add_option("-p", "--proteinortho_file", default="./myproject.proteinortho")
parser.add_option("-d", "--fasta_dir", help="directory with fastas")
parser.add_option("-o", "--outfile_name", help="Name of the output file", default="singletones.myproject.proteinortho")
opt, args = parser.parse_args()

getsingletones(opt.proteinortho_file, opt.fasta_dir, opt.outfile_name)
