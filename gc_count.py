from os import listdir
import optparse


def gc_count(seqfile):
	gc_dict = {'A':0,'T':0,'G':0,'C':0}
	for q in open(seqfile):
		q = q.strip()
		if q[0] == '>':
			continue
		for c in q:
			if c in gc_dict.keys():
				gc_dict[c] += 1
	return float(gc_dict['G'] + gc_dict['C'])/(gc_dict['A']+gc_dict['C']+gc_dict['G']+gc_dict['C'])


def gc_count_all(seqdir):
	for fasta in listdir(seqdir):
		print fasta, gc_count(seqdir + fasta)


parser = optparse.OptionParser()
parser.add_option("-d", "--seqdir", help="Directory with fastas")
opt, args = parser.parse_args()

gc_count_all(opt.seqdir)