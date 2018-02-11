"""
permut_partition_adlib
additional library for permut_partition.py
"""

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


def find_phyletic_pattern(peak_pos, proteinortho_file):
	phylpat_dict = dict()
	spc_dict = dict()
	for q in open(proteinortho_file):
		q = q.strip().split('\t')
		if q[0][0] == '#':
			count = 0
			for p in q[3:]:
				spc_dict[count] = p.split('.')[0]
				count += 1
			continue
		if eval(q[0]) != peak_pos:
			continue
		pat_str = ""
		for p in q[3:]:
			if p == '*':
				pat_str += '-'
			else:
				pat_str += '+'
		try:
			phylpat_dict[pat_str] += 1
		except:
			phylpat_dict[pat_str] = 1
	maxnum = 0
	for pat_str in phylpat_dict.keys():
		if phylpat_dict[pat_str] > maxnum:
			max_pat_str = pat_str
			maxnum = phylpat_dict[pat_str]
	pos_phyletic_pattern = []
	neg_phyletic_pattern = []

	return {"spc_dict":spc_dict, "phylpat_dict":phylpat_dict, "max_pat_str":max_pat_str}