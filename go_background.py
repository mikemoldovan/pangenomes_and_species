"""
Pangenome gene sampler for GO-analysis background
Input:
1. proteinortho_file
2. Directory with protein fastas
3. Output file name
"""

from optparse import OptionParser

def spclist(proteinortho_file):
	for q in open(proteinortho_file):
		q = q.strip().split('\t')
		return [p for p in q[3:]]


def fastawrite(outfile_name, spc_list, fasdir):
	outfile = open(outfile_name, "w")
	for spc in spc_list:
		for q in open(fasdir + spc):
			outfile.write(q)

parser = OptionParser()
parser.add_option("-p", "--proteinortho_file", default="./myproject.proteinortho")
parser.add_option("-f", "--fasdir",help = "directory with protein fastas")
parser.add_option("-o", "--outfile_name",help = "Name of the output file")
opt, args = parser.parse_args()

spc_list = spclist(opt.proteinortho_file)
fastawrite(opt.outfile_name, spc_list, opt.fasdir)

