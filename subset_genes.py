"""
Specific genes barplot
input:
1. File in format:
SPC<n>	k
-without header
-k is the group number

2. Proteinortho file

Output:
-file with all combinations and numbers
"""

from optparse import OptionParser

def infile_parse(infile, proteinortho_file):
	group_dict = {q.strip().split()[0]:q.strip().split()[1] for q in open(infile)}
	num_dict = dict()
	for q in open(proteinortho_file):
		q = q.strip().split('\t')
		count = 0
		for p in q[3:]:
			num_dict[count + 3] = group_dict[p.split('.')[0]]
			count += 1
		break
	return num_dict


def makesubsets(maxnum):
	subsets = dict()
	refnum = bin(2**(maxnum + 1))[3:]
	l = len(refnum)
	for i in range(1, 2**(maxnum + 1)):
		bini = bin(i)
		bini_num = str(bini[2:])
		bini = (l - len(bini_num))*'0' + bini_num
		subset = []
		for j in range(len(bini)):
			if bini[j] == '1':
				subset.append(str(j))
		print subset
		subsets[tuple(subset)] = 0
	return subsets


def step_str(proteinortho_file_str, subsets_dict, num_dict):
	str_subset = []
	for i in range(3, len(proteinortho_file_str)):
		if proteinortho_file_str[i] != '*':
			if num_dict[i] not in str_subset:
				str_subset.append(num_dict[i])
		else:
			if num_dict[i] in str_subset:
				return None
	ref = tuple(sorted(str_subset))
	if ref:
		subsets_dict[ref] += 1

def main_process(proteinortho_file, groupfile):
	num_dict = infile_parse(groupfile, proteinortho_file)
	vals = []
	for i in num_dict.values():
		vals.append(eval(i))
	subsets_dict = makesubsets(max(vals))
	print subsets_dict
	for q in open(proteinortho_file):
		if q[0] == '#':
			continue
		q = q.strip().split()
		step_str(q, subsets_dict, num_dict)
	return subsets_dict

def print_res(subsets_dict):
	for k in subsets_dict.keys():
		print '&'.join(list(k)) + '\t' + str(subsets_dict[k])


parser = OptionParser()
parser.add_option("-i", "--groupfile", help="File with SPCs with assigned groups")
parser.add_option("-p", "--proteinortho_file", help="Proteinortho myproject file")
opt, args = parser.parse_args()

subsets_dict = main_process(opt.proteinortho_file, opt.groupfile)
print_res(subsets_dict)


print makesubsets(3)
