"""
Pangenome matrix construction according to a pattern
Input:
1. Pattern (+/-) -s
2. Proteinortho file -p
3. Singletone list -n
Output:
mat_x	mat_y	ogg_n	genes = |OGG1|OGG2|...	
"""

from optparse import OptionParser

class Mat_string():
	def __init__(self, ogg_n, genes_str):
		self.ogg_n = ogg_n
		self.genes_str = genes_str


def retrieve_spcs(proteinortho_file, pattern):
	spcs = dict()
	for q in open(proteinortho_file):
		q = q.strip().split('\t')
		if q[0][0] == '#':
			for i in range(len(q[3:])):
				spcs[q[i].split('.')[0]] = pattern[i]
			return spcs

def countpos(pattern):
	pos_plus = 0
	pos_minus = 0
	for c in pattern:
		if c == '+':
			pos_plus += 1
		else:
			pos_minus += 1
	return pos_plus, pos_minus

def ogg_genes(proteinortho_file_str, pattern):
	genes_plus = []
	genes_minus = []
	for i in range(len(pattern)):
		if proteinortho_file_str[i] == '*':
			continue
		if pattern[i] == '+':
			genes_plus.append(proteinortho_file_str[i])
		else:
			genes_minus.append(proteinortho_file_str[i])
#	print genes_plus, genes_minus
	return genes_plus, genes_minus


def buildmatdict_proteinortho_file(proteinortho_file, pattern):
	matdict = dict()
	for q in open(proteinortho_file):
		if q[0] == '#':
			continue
		q = q.strip().split('\t')
		genes_plus, genes_minus = ogg_genes(q[3:], pattern)
		genes_plus_n = len(genes_plus)
		genes_minus_n = len(genes_minus)
		try:
			matdict[(genes_plus_n,genes_minus_n)].ogg_n += 1
			genstr = '^' + ','.join(genes_plus + genes_minus)
			matdict[(genes_plus_n,genes_minus_n)].genes_str += genstr
		except:
			genstr = ','.join(genes_plus + genes_minus)
			matdict_val_init = Mat_string(ogg_n = 1, genes_str = genstr)
			matdict[(genes_plus_n,genes_minus_n)] = matdict_val_init
	return matdict


def append_3_cases(singletone_list, proteinortho_file, matdict, pattern):
	spcs = retrieve_spcs(proteinortho_file, pattern)
	singletones = [q.strip() for q in open(singletone_list)]
	plus_single_genes = Mat_string(ogg_n=0, genes_str = "*")
	minus_single_genes = Mat_string(ogg_n=0, genes_str = "*")
	for spc in singletones:
		try:
			sign = spcs[spc.split('|')[0]]
		except:
			sign = 0
		if sign == 0:
			continue
		genstr = '^' + spc
		if sign == '+':
			plus_single_genes.ogg_n += 1
			plus_single_genes.genes_str += genstr
		else:
			minus_single_genes.ogg_n += 1
			minus_single_genes.genes_str += genstr
	
	try:
		matdict[(1,0)].ogg_n += plus_single_genes.ogg_n
		plus_single_genes.genes_str = plus_single_genes.genes_str[1:]
		matdict[(1,0)].genes_str += plus_single_genes.genes_str
	except:
		plus_single_genes.genes_str = plus_single_genes.genes_str[2:]
		matdict[(1,0)] = plus_single_genes

	try:
		matdict[(0,1)].ogg_n += minus_single_genes.ogg_n
		minus_single_genes.genes_str = minus_single_genes.genes_str[1:]
		matdict[(0,1)].genes_str += minus_single_genes.genes_str
	except:
		minus_single_genes.genes_str = minus_single_genes.genes_str[2:]
		matdict[(0,1)] = minus_single_genes

	zero_genes = Mat_string(ogg_n=0,genes_str = '*')
	matdict[(0,0)] = zero_genes
	return matdict


def print_matdict(matdict, pattern, print_genes):
	pos_plus, pos_minus = countpos(pattern)
	if print_genes=='y':
		print "pos_val\tneg_val\togg_number\tgenes"
	else:
		print "pos_val\tneg_val\togg_number"
	for i in range(pos_plus + 1):
		for j in range(pos_minus + 1):
			try:
				string = matdict[(i,j)]
			except:
				string = Mat_string(ogg_n=0, genes_str= '*')
			if print_genes=='y':
				print str(i)+'\t'+str(j)+'\t'+str(string.ogg_n) +'\t'+string.genes_str
			else:
				print str(i)+'\t'+str(j)+'\t'+str(string.ogg_n)


parser = OptionParser()
parser.add_option("-s", "--phylpat", help="Phyletic pattern string (+-)")
parser.add_option("-p", "--proteinortho_file", help="File with proteinortho matrix")
parser.add_option("-n", "--singletone_list", help="List with singletones")
parser.add_option("-g", "--print_genes", help="Print the list of all genes in a peak? (y/n)", default='n')
opt, args = parser.parse_args()

matdict = buildmatdict_proteinortho_file(opt.proteinortho_file, opt.phylpat)
matdict = append_3_cases(opt.singletone_list, opt.proteinortho_file, matdict, opt.phylpat)
print_matdict(matdict, opt.phylpat, opt.print_genes)


"""
pattern = "-------+---+--------++---+--+-----++----+-"
proteinortho_file = "2/myproject.proteinortho"
buildmatdict_proteinortho_file(proteinortho_file, pattern)
"""
