"""
singletone_list_construct
input:
1. Proteinortho_file
2. Fasta_dir
(3. phyletic pattern)
"""

from optparse import OptionParser

def get_fasta_names(proteinortho_file, phyletic_pattern):
	spec_phylpat_list = []
	for q in open(proteinortho_file):
		q = q.strip()
		if q[0] == '#':
			q = q.strip().split('\t')
			fastalist = q[3:]
			for i in range(len(fastalist)):
				spec_id = fastalist[i].split('.')[0]
				if phyletic_pattern[i] == '+':
					spec_phylpat_list.append(spec_id)
			return fastalist, spec_phylpat_list


def get_ogg_genes(proteinortho_file):
	genedict = dict()
	for q in open(proteinortho_file):
		q = q.strip().split('\t')
		for p1 in q[3:]:
			for p2 in p1.split(','):
				genedict[p2] = 0
	return genedict


def getsingletones(proteinortho_file, fasta_dir, outfile_name, phyletic_pattern):
	fasta_names, spec_phylpat_list = get_fasta_names(proteinortho_file, phyletic_pattern)
	genedict_ogg = get_ogg_genes(proteinortho_file)
	outfile = open(outfile_name, "w")
	for fasta in fasta_names:
		for string in open(fasta_dir+fasta):
			if string[0] == '>':
				_id = string.strip()[1:]
				try:
					k = genedict_ogg[_id]
				except:
					spec_id = _id.split('|')[0]
					if spec_id in spec_phylpat_list:
						pat = "pat1"
					else:
						pat = "pat2"
					outfile.write(_id + '\t'+pat+'\n')
					print '1\t'+pat
	outfile.close()


parser = OptionParser()
parser.add_option("-p", "--proteinortho_file", default="./myproject.proteinortho")
parser.add_option("-d", "--fasta_dir", help="directory with fastas")
parser.add_option("-o", "--outfile_name", help="Name of the output file", default="singletones.myproject.proteinortho")
parser.add_option("-s", "--phyletic_pattern")
opt, args = parser.parse_args()

getsingletones(opt.proteinortho_file, opt.fasta_dir, opt.outfile_name, opt.phyletic_pattern)
