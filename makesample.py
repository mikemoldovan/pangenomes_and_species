"""
select strains
build the genome size distribution and select strains by the following criteria:

1. No repeating genome sizes
2. Deviation no more than 3SDs
3. Sample size not above the given value

Input:
1. Directory with multiple gzipped fastas
2. SD threshold number | def=3
3. Maximal sample volume | def=100
4. Sample number | def=1
5. Bin: delete initial directory? y/n | def=n
6. Bin: exclude overlapping samples? y/n | def=n
"""

import gzip
import os
import random
import optparse
from Bio import SeqIO


def mean(arr):
	return sum(arr)/len(arr)

def sd(arr):
	arrmean = mean(arr)
	sd = 0
	for i in arr:
		sd += (i - arrmean)*(i-arrmean)
	sd = (sd/(len(arr) - 1))**0.5
	return sd

def fastacount_z(fafile_z):
	prots = 0
	with gzip.open(fafile_z, "rb") as fafile:
		for q in fafile:
			if q[0] == '>': prots += 1
	return prots


def makeprotlendict(fasta_dir): #fasta_dir with /
	protlendict = dict()
	for q in os.listdir(fasta_dir):
		if ".faa.gz" not in q:
			continue
		count = fastacount_z(fasta_dir + q)
		protlendict[count] = q
	return protlendict


#--------Select fasta functions--------

def select_by_sd(protlenarr, sd_num):
	arrmean = mean(protlenarr)
	arrsd = sd(protlenarr)
	sel_arr = []
	for i in protlenarr:
		if (i > arrmean - sd_num*arrsd) and (i < arrmean + sd_num*arrsd):
			sel_arr.append(i)
	return sel_arr


def select_1_sample(protlenarr, max_volume, exclude_overlapping):
	l = len(protlenarr)
	if l == 0:
		return 0
	elif l <= max_volume:
		return 1
	else:
		pass

	outarr = []
	if exclude_overlapping == 'n':
		protlenarr_copy = protlenarr[:]
	else:
		protlenarr_copy = protlenarr

	for i in range(max_volume):
		r = random.choice(range(len(protlenarr_copy)))
		outarr.append(protlenarr_copy.pop(r))
	return outarr


def makesamples(protlenarr, max_volume, n_samples, exclude_overlapping):
	samples = []
	for i in range(n_samples):
		sample = select_1_sample(protlenarr, max_volume, exclude_overlapping)
		if sample == 0:
			return samples
		elif sample == 1:
			samples.append(protlenarr)
			return samples
		else:
			samples.append(sample)
	return samples

def copyfasta(fafile_init, fafile_new, dest_dir):
	handle_in = gzip.open(fafile_init)
	handle_out = open(dest_dir + fafile_new + '.fa', "w")
	genetable = open(dest_dir + fafile_new + '_genetable.txt', "w")

	genetable.write("init_ID\tnew_ID\n")

	count = 0
	for record in SeqIO.parse(handle_in, "fasta"):
		count += 1
		genetable.write(record.id + '\t' + fafile_new + '|' + 'g' + str(count) + '\n')
		record.id = fafile_new + '|' + 'g' + str(count)
		record.name = ''
		record.description = ''
		SeqIO.write(record, handle_out, "fasta")
	handle_out.close()
	genetable.close()


def select_fastas(fasta_dir, sd_num, max_volume, n_samples, exclude_overlapping, spec_id):
	fd = fasta_dir.split('/')[-2]
	try:
		os.mkdir(fd + '_sampled')
	except:
		pass
	protlendict = makeprotlendict(fasta_dir)
	protlenarr = list(protlendict.keys())
	protlenarr = select_by_sd(protlenarr, sd_num)
	samples = makesamples(protlenarr, max_volume, n_samples, exclude_overlapping)

	for i in range(len(samples)):
		dirname = fd + '_sampled/' + str(i+1) + '/'
		try:
			os.mkdir(dirname)
		except:
			pass
		spectable = open(dirname + "spectable_"+spec_id+".txt", "w")
		proteinortho_command = open(dirname + "proteinortho.command_" + spec_id, "w")
		spectable.write("old_name\tnew_name\n")
		count = 0
		for genlen in samples[i]:
			count += 1
			spec = protlendict[genlen]
			spectable.write(spec + '\t' + spec_id + str(count) + '\n')
			proteinortho_command.write(spec_id + str(count) + '.fa\t')
			fafile_init = fasta_dir + spec
			fafile_new = spec_id + str(count)
			copyfasta(fafile_init, fafile_new, dirname)
		spectable.close()
		proteinortho_command.close()


def remove_initial(fasta_dir):
	os.system("rm -R " + fasta_dir)


parser = optparse.OptionParser()
parser.add_option("-f", "--fasta_dir", help="Directory with multiple gzipped fastas", default="a")
parser.add_option("-s", "--sd_num", help="SD threshold number", default=3)
parser.add_option("-v", "--max_volume", help="Maximal sample volume", default=100)
parser.add_option("-n", "--n_samples", help="Sample number", default=1)
parser.add_option("-o", "--exclude_overlapping", help="Bin: delete initial directory? y/n", default='n')
parser.add_option("-r", "--remove_initial", help="Bin: exclude overlapping samples? y/n", default='n')
parser.add_option("-i", "--spec_id", help="short species id", default="SPC")
opt, args = parser.parse_args()

#print opt.max_volume
select_fastas(opt.fasta_dir, float(opt.sd_num), int(opt.max_volume), int(opt.n_samples), opt.exclude_overlapping, opt.spec_id)

if opt.remove_initial == 'y':
	remove_initial(opt.fasta_dir)


#protlendict = makeprotlendict(args.fasta_dir)
#sel_arr = select_lengths(protlendict, args.sd_num, args.maxgroup)
#remove_files(args.fasta_dir, sel_arr, protlendict)


"""
def select_lengths(protlendict, sd_num, num_threshold):
	sel_arr = []
	protlenarr = list(protlendict.keys())
	arrmean = mean(protlenarr)
	arrsd = sd(protlenarr)
	for i in protlenarr:
		if (i > arrmean - sd_num*arrsd) and (i < arrmean + sd_num*arrsd):
			sel_arr.append(i)

	lsel = len(sel_arr)
	if lsel <= num_threshold:
		return sel_arr

	fin_arr = []
	for i in range(num_threshold):
		fin_arr.append(random.choice(sel_arr))
	return fin_arr
"""
