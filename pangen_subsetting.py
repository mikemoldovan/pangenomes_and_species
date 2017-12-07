"""
Pangenome sampler
Rebuild pangenome from a subset
"""

from optparse import OptionParser

def protortho_file_sample(proteinortho_file, selected_nums, outfile_name):
	outfile = open(outfile_name, 'w')
	selected_nums = [0,1,2] + [i + 3 for i in selected_nums]
	for q in open(proteinortho_file):
		count = 0
		str_arr = []
		q = q.strip()
		if q[0] == '#':
			q = q.split('\t')
			for i in range(len(q)):
				if i in selected_nums:
					str_arr.append(q[i])
			outfile.write('\t'.join(str_arr) + '\n')
			continue
		for p in q.split('\t'):
			if count in selected_nums:
				str_arr.append(p)
			count += 1
		spc_num = 0
		gen_num = 0
		for gene_group in str_arr[3:]:
			if gene_group != '*':
				spc_num += 1
				for gene in gene_group.split(','):
					gen_num += 1
		str_arr[0] = str(spc_num)
		str_arr[1] = str(gen_num)
		outfile.write('\t'.join(str_arr) + '\n')
	outfile.close()


def select_from_phylpat(phylpat):
	selected_nums = []
	for i in range(len(phylpat)):
		if phylpat[i] == '+':
			selected_nums.append(i)
	return selected_nums


def select_from_spclist(spclist, proteinortho_file):
	spcs = [q.strip() for q in open(spclist)]
	for q in open(proteinortho_file):
		if q[0] == '#':
			q = q.split('\t')
			selected_nums = []
			count = 0
			for p in q[3:]:
				if p.split('.')[0] in spcs:
					selected_nums.append(count)
				count += 1
			return selected_nums

parser = OptionParser()
parser.add_option("-p", "--proteinortho_file", default="./myproject.proteinortho")
parser.add_option("-a", "--phylpat", help="phylpat [optional]", default=None)
parser.add_option("-s", "--spclist", help="list with spc column [optional]", default=None)
parser.add_option("-o", "--outfile_name", help="Name of the output file", default="myproject_2.proteinortho")
opt, args = parser.parse_args()

if opt.phylpat:
	selected_nums = select_from_phylpat(opt.phylpat)
else:
	selected_nums = select_from_spclist(opt.spclist, opt.proteinortho_file)


protortho_file_sample(opt.proteinortho_file, selected_nums, opt.outfile_name)

