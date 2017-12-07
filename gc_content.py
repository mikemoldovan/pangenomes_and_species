"""
GC content count
"""

from os import listdir
from optparse import OptionParser

def gc_content_dict(filelist, alphabet):
	gc_contents = dict()
	for fafile in filelist:
		gc_content = dict()
		for c in alphabet:
			gc_content[c] = 0
		for s in open(fafile):
			s = s.strip()
			if s[0] == '>':
				continue
			for c in s:
				try:
					gc_content[c] += 1
				except:
					pass
		gc_contents[fafile.split('/')[-1]] = gc_content
#		print fafile
	return gc_contents

def main_pr(fadir, alphabet):
	filelist = [fadir + q for q in listdir(fadir) if "_n.fa" in q]
	gc_contents = gc_content_dict(filelist, alphabet)
	for k1 in gc_contents.keys():
		outstr = []
		for c in alphabet:
			outstr.append(gc_contents[k1][c])
		s = sum(outstr)
		for i in range(len(outstr)):
			outstr[i] = str(float(outstr[i])/s)
		print k1.split('_')[0] + '\t' + '\t'.join(outstr)


parser = OptionParser()
parser.add_option("-f", "--fasta_dir", help = "Directory nucleotide fastas")
parser.add_option("-a", "--alphabet", help = "Alphabet (def. 'atgc')", default = "atgc")
opt, args = parser.parse_args()

main_pr(opt.fasta_dir, opt.alphabet)