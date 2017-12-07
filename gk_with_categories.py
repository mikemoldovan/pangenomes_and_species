"""
Pangenome partial build
Returns the list of OGGs with the phyletic pattern it belongs to
"""

from optparse import OptionParser


def phylpat_to_numlist(phylpat):
	numlist = []
	for i in range(len(phylpat)):
		if phylpat[i] == '+':
			numlist.append(i + 3)
	return numlist

def protortho_string_work(string, numlist1, numlist2):
	only_1 = True
	only_2 = True
	count = 0
	for i in range(3, len(string)):
		if string[i] != '*':
			count += 1
			if i in numlist1:
				only_2 = False
			elif i in numlist2:
				only_1 = False
	if (not only_1) and (not only_2):
		print str(count) + '\t' + "univ"
	elif only_1:
		print str(count) + '\t' + "pat1"
	elif only_2:
		print str(count) + '\t' + "pat2"
	else:
		print str(count) + '\t' + "mixd"

#phylpat -> numpat
def proteinortho_file_work(proteinortho_file, phylpat1, phylpat2):
	numlist1 = phylpat_to_numlist(phylpat1)
	numlist2 = phylpat_to_numlist(phylpat2)
	print "k\ttype"
	for q in open(proteinortho_file):
		if q[0] == '#':
			continue
		q = q.strip().split('\t')
		protortho_string_work(q, numlist1, numlist2)

parser = OptionParser()
parser.add_option("-a", "--phylpat1", help="Phyletic pattern 1")
parser.add_option("-b", "--phylpat2", help="Phyletic pattern 2")
parser.add_option("-p", "--proteinortho_file", help="File with proteinortho matrix")
opt, args = parser.parse_args()


proteinortho_file_work(opt.proteinortho_file, opt.phylpat1, opt.phylpat2)

# -------+---+--------++---+--+-----++----+-
# +++++++-+++-++++++++--+++-++-+++++--++++-+