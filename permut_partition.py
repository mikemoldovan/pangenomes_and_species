"""
permut_partition
Obtain the best partition by permutations

Input:
1. Proteinortho file [proteinortho_file] (myproject.proteinortho)
2. Phyletic pattern [phylpat]
3. Maximal number of permutations for one sample size [n] (10)
4. Precision -- parameter for the random sampling: it will stop after p iterations that
   yielded strings that were yielded before

Output:
file in format:
<number of strains in partition>\t<peak height>\t<category: corresponds to phyletic pattern [y] or not [n]>


prchl
-------+---+--------++---+--+-----++----+-
"""

import numpy
import os
from optparse import OptionParser
from permut_partition_adlib import *
#selected_nums = select_from_phylpat(phylpat)
#protortho_file_sample(proteinortho_file, selected_nums, outfile_name)
#find_phyletic_pattern(peak_pos, proteinortho_file)
#return {"spc_dict":spc_dict, "phylpat_dict":phylpat_dict, "max_pat_str":max_pat_str}


def peak_position(phylpat_init, phylpat_subsamp):
	l = len(phylpat_init)
	phylpat_new = ""
	for i in range(l):
		if phylpat_subsamp[i] == '+' and phylpat_init[i] == '+':
			phylpat_new += '+'
		elif phylpat_subsamp[i] == '+' and phylpat_init[i] == '-':
			phylpat_new += '-'
	return phylpat_new


def make_subsamples(phylpat, samp_size, n, precision):
	subsamps = dict()
	subsamp_num = 0
	failures = 0
	stop = False
	l = len(phylpat)
	pos_arr = list(range(l))
	while not stop:
		pa = pos_arr[:]
		random_positions = numpy.random.choice(pa, size=samp_size, replace=False)
		phylpat_subsamp = ""
		count = 0
		for i in range(l):
			if i in random_positions:
				phylpat_subsamp += '+'
			else:
				phylpat_subsamp += '-'
		if subsamps.get(phylpat_subsamp):
			failures += 1
			if failures == precision:
				stop = True
		else:
			failures = 0
			phylpat_new = peak_position(phylpat, phylpat_subsamp)
			if phylpat_new.count('+') < 2:
				continue
			if phylpat_new.count('-') < 2:
				continue
			subsamps[phylpat_subsamp] = phylpat_new
			subsamp_num += 1
		if subsamp_num >= n:
			stop = True
#phylpat_subsamp : phyletic pattern corresponding to initial partition
	return subsamps


def reverse_phylpat(phylpat):
	res = phylpat.replace('+', '1')
	res = res.replace('-','+')
	res = res.replace('1','-')
	return res


def main(proteinortho_file, phylpat, n, precision):
	l = len(phylpat)
	for i in range(4, l + 1):
		subsamps = make_subsamples(phylpat, i, n, precision)
		count_subsamps = 0
		for subsamp_phylpat in subsamps.keys():
			outfile_name = proteinortho_file + ".{}_{}".format(i, count_subsamps)
			selected_nums = select_from_phylpat(subsamp_phylpat)
			protortho_file_sample(proteinortho_file, selected_nums, outfile_name)

			peak_pos_forward = subsamps[subsamp_phylpat].count('+')
			peak_pos_reverse = subsamps[subsamp_phylpat].count('-')

			phylpat_res_forward = find_phyletic_pattern(peak_pos_forward, outfile_name)
			phylpat_res_reverse = find_phyletic_pattern(peak_pos_reverse, outfile_name)

			phylpat_res_phylpat_dict_concat = dict()
			for k in phylpat_res_forward["phylpat_dict"].keys():
				rev_phylpat_freq = phylpat_res_reverse["phylpat_dict"].get(reverse_phylpat(k))
				if rev_phylpat_freq:
					phylpat_freq = phylpat_res_forward["phylpat_dict"][k] + phylpat_res_reverse["phylpat_dict"][reverse_phylpat(k)]
				else:
					phylpat_freq = phylpat_res_forward["phylpat_dict"][k]
				phylpat_res_phylpat_dict_concat[k] = phylpat_freq

			for k in phylpat_res_phylpat_dict_concat.keys():
				if k == subsamps[subsamp_phylpat] or k == reverse_phylpat(subsamps[subsamp_phylpat]):
					flag = 'y'
				else:
					flag = 'n'
				print str(i) + '\t' + str(phylpat_res_phylpat_dict_concat[k]) + '\t' + flag + '\t' + k
			os.remove(outfile_name)
			count_subsamps += 1


parser = OptionParser()
parser.add_option("-p", "--proteinortho_file", default="./myproject.proteinortho")
parser.add_option("-a", "--phylpat", help="phylpat [optional]")
parser.add_option("-n", "--maxperm", help="maximal number of permutations", default="10")
parser.add_option("-r", "--precision", help="Maximal number of tries to construct a new phyletic pattern", default="10")
opt, args = parser.parse_args()


main(opt.proteinortho_file, opt.phylpat, eval(opt.maxperm), eval(opt.precision))
#print make_subsamples(opt.phylpat, 20, eval(opt.maxperm), eval(opt.precision))
