"""
ProteinOrtho G(k) builder
Builds G(k) table and finds the most common phyletic pattern corresponding to the highest peak
"""

from sys import stdout, stderr
from os import listdir, mkdir, system
from Bio import SeqIO
from optparse import OptionParser
import append_outgroup

def fastacount(fafile):
	count = 0
	for q in fafile:
		if q[0] == '>':
			count += 1
	return count

#~~~~~~~ main: gk_build(*) Returns the G(k) spectrum as a dictionary ~~~~~~~~ 
def gk_prebuild(proteinortho_file):
	gk_dict = dict()
	for q in open(proteinortho_file):
		if q[0] == '#':
			continue
		q = q.strip().split()
		nspec = eval(q[0])
		try:
			gk_dict[nspec] += 1
		except:
			gk_dict[nspec] = 1
	return gk_dict

def append_singletones(singletone_table, proteinortho_file):
	for q in open(proteinortho_file):
		if q[0] == '#':
			spclist = []
			q = q.split('\t')
			for p in q[3:]:
				spclist.append(p.split('.')[0])
			break

	singletones = 0
	for q in open(singletone_table):
		q = q.strip().split('|')
		if q[0] in spclist:
			singletones += 1
	return singletones

def gk_build(proteinortho_file, singletone_table):
	gk_dict = gk_prebuild(proteinortho_file)
	singletones = append_singletones(singletone_table, proteinortho_file)
	try:
		gk_dict[1] += singletones
	except:
		gk_dict[1] = singletones
	gk_dict2 = dict()
	for k in gk_dict.keys():
		if k != 0:
			gk_dict2[k] = gk_dict[k]
	return gk_dict2

def gk_print(gk_dict, outfile):
	outfile.write("strain_num\tgk_value\n")
	for k in sorted(gk_dict.keys()):
		outfile.write(str(k)+'\t'+str(gk_dict[k])+'\n')
	outfile.close()

#~~~~~~~ Find the highest inner peak and build the corresponding phyletic pattern ~~~~~~
def findpeak(gk_dict):
	highest_peak = 0   
	highest_peak_pos = 0
	for i in range(2, max(gk_dict.keys())-1):
		if (gk_dict[i] > gk_dict[i-1]) and (gk_dict[i] > gk_dict[i+1]):
			peak = gk_dict[i] - (gk_dict[i-1] + gk_dict[i+1])/2
			if peak > highest_peak:
				highest_peak_pos = i
				highest_peak = peak
	return highest_peak, highest_peak_pos

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


#~~~~~~~~~ Test if there actually is a second peak ~~~~~~~

def print_fpp(fpp):
	print "spc_dict"
	for k in fpp["spc_dict"].keys():
		print k, fpp["spc_dict"][k]
	print "\nphylpat_dict"
	for k in fpp["phylpat_dict"].keys():
		print k, fpp["phylpat_dict"][k]
	print "\nmax_pat_str"
	print fpp["max_pat_str"]


def build_reverse_pattern(pattern):
	revpat = ""
	for c in pattern:
		if c == '+':
			revpat += '-'
		else:
			revpat += '+'
	return revpat

def reverse_pattern_test(peak_pos, proteinortho_file):
	fpp = find_phyletic_pattern(peak_pos, proteinortho_file)
	revpat = build_reverse_pattern(fpp["max_pat_str"])
	rev_peak_pos = len(revpat) - peak_pos
	rev_fpp_pat = find_phyletic_pattern(rev_peak_pos, proteinortho_file)["max_pat_str"]
	return rev_fpp_pat == revpat

#~~~~~~~~~ Retrieve genes from the peak ~~~~~~~
def ispat(proteinortho_file_str, pat_str):
	count = 0
	for q in proteinortho_file_str[3:]:
		if q == '*' and pat_str[count] == '+':
			return False
		elif q != '*' and pat_str[count] == '-':
			return False
		count += 1
	return True

def retrieve_gene_ids(peak_pos, proteinortho_file):
	fpp = find_phyletic_pattern(peak_pos, proteinortho_file)
	gene_ids = []
	for q in open(proteinortho_file):
		if q[0] == '#':
			continue
		q = q.strip().split()
		if not ispat(q, fpp["max_pat_str"]):
			continue
		for p in q[3:]:
			if p != '*':
				gene_ids += p.split(',')
	return gene_ids, fpp

def getseqs(peak_pos, proteinortho_file, seq_dir, outfile_name="peak_seqs.fa"):
	gene_ids, fpp = retrieve_gene_ids(peak_pos, proteinortho_file)
#	print gene_ids
	outfile = open(outfile_name, "w")
	for q in listdir(seq_dir):
		if (q.split('.')[-1] != "fa") or ("_n.fa" in q):
			print q, q.split('.')[-1]
			continue
		handle = open(seq_dir + q)
		for record in SeqIO.parse(handle, "fasta"):
			if record.id in gene_ids:
				SeqIO.write(record, outfile, "fasta")
		handle.close()
	outfile.close()
	print "Hello!"

#~~~~~~ Make alignments for tree ~~~~~~~
def retrieve_gene_ids2(proteinortho_file_str):
	gene_ids = []
	if eval(proteinortho_file_str[0]) != eval(proteinortho_file_str[1]):
		return None
	for q in proteinortho_file_str[3:]:
		gene_ids.append(q)
		if q == '*':
			return None
	return gene_ids


def makeal(gene_ids, seq_dir, al_dir_name, alname, phylpat_str, outgroup_db):
	multal = open(al_dir_name + alname + ".raw", 'w')
	for i in range(len(gene_ids)):
		gene_id = gene_ids[i]
		spc_nfa = seq_dir + gene_id.split('|')[0] + "_n.fa"
		handle = open(spc_nfa)
		for record in SeqIO.parse(handle, "fasta"):
			if record.id == gene_id:
				record.id = phylpat_str[i] + gene_id.split('|')[0]
				record.name = None
				record.descripption = None
				SeqIO.write(record, multal, "fasta")
				break
	multal.close()
	if outgroup_db:
		res = append_outgroup.append_outgroup(outgroup_db, al_dir_name + alname + ".raw")
#		if res:
#			system("muscle -in "+al_dir_name+alname+".raw"+" -out "+al_dir_name+alname)
		return res
#	system("muscle -in "+al_dir_name+alname+".raw"+" -out "+al_dir_name+alname)
	return 1


def makeals(proteinortho_file, seq_dir, al_dir_name, phylpat_str, fuse_aligns=True, use_aligns="all", outfile="catalign.fa", blastdb=None):
	try:
		mkdir(al_dir_name)
	except:
		pass
	count = 0
	for q in open(proteinortho_file):
		q = q.strip()
		if q[0] == '#':
			continue
		gene_ids = retrieve_gene_ids2(q.split())
		if not gene_ids:
			continue
		print gene_ids
		app = makeal(gene_ids, seq_dir, al_dir_name, str(count)+".fa", phylpat_str, blastdb)
		if app == 0:
			continue
		count += 1
		if use_aligns != "all":
			if use_aligns == count:
				break
#	system("rm "+al_dir_name+"*.raw")
#	if fuse_aligns:
#		system("python alignment_fusion.py -d "+al_dir_name+" -o "+outfile)


#~~~~~~ Make G(k)s for pattern ~~~~~~~
def make_pat_gk(proteinortho_file, seq_dir, peak_pos):
	fpp = find_phyletic_pattern(peak_pos, proteinortho_file)
	max_pat_str = fpp["max_pat_str"]
	spc_dict = fpp["spc_dict"]
	gk = dict()
	spc_fa = [spc_dict[i]+'.fa' for i in sorted(spc_dict.keys())[1:] if max_pat_str[i-1] == '+']
	pangen_arr = []
	for q in open(proteinortho_file):
		if q[0] == '#':
			continue
		num_pos = 0
		q = q.split()
		for i in range(len(q[3:])):
			if max_pat_str[i-3] == '+' and q[i] != '*':
				num_pos += 1
				pangen_arr += q[i].split(',')
		if num_pos != 0:
			try:
				gk[num_pos] += 1
			except:
				gk[num_pos] = 1
	count_all = 0
	for q in spc_fa:
		fafile = seq_dir + q
		count_all += fastacount(open(fafile))
		print fastacount(open(fafile))
	pangen_sum = len(set(pangen_arr))
	print count_all, len(pangen_arr), len(set(pangen_arr))
	print gk
	try:
		gk[1] += count_all - pangen_sum
	except:
		gk[1] = count_all - pangen_sum
	print gk
	return gk

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~ Task functions ~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Print G(k) spectrum to a given file
def gk_spectrum(proteinortho_file, singletone_table, outfile):
        print "hello!"
	if (outfile != stdout) and (outfile != stderr):
		outfile = open(outfile, 'w')
	gk_dict = gk_build(proteinortho_file, singletone_table)
	gk_print(gk_dict, outfile)
	if (outfile != stdout) and (outfile != stderr):
		outfile.close()

# Print G(k) spectrum of a subset
def gk_spectrum2(proteinortho_file, seq_dir, peak_pos, outfile):
	if (outfile != stdout) and (outfile != stderr):
		outfile = open(outfile, 'w')
	gk_dict = make_pat_gk(proteinortho_file, seq_dir, peak_pos)
	gk_print(gk_dict, outfile)
	if (outfile != stdout) and (outfile != stderr):
		outfile.close()

# Print peak info: 
# 1. Phyletic pattern distribution
# 2. The most common phyletic pattern
# 3. The double peak test
def peak_info(peak_pos, proteinortho_file):
	fpp = find_phyletic_pattern(peak_pos, proteinortho_file)
	print_fpp(fpp)
	print "\nREVERSE PATTERN TEST"
	print reverse_pattern_test(peak_pos, proteinortho_file)

# Get genes prom the most common phyletic pattern in position peak_pos
def peakgenes(peak_pos, proteinortho_file, seq_dir, outfile_name):
	getseqs(peak_pos, proteinortho_file, seq_dir, outfile_name)

# Make alignments of universal genes [for further tree building]
def tree_alignment(proteinortho_file, seq_dir, al_dir_name, fuse_aligns, use_aligns, outfile, peak_pos, blastdb):
	phylpat_str = find_phyletic_pattern(peak_pos, proteinortho_file)["max_pat_str"]
	makeals(proteinortho_file, seq_dir, al_dir_name, phylpat_str, fuse_aligns, use_aligns, outfile, blastdb)


parser = OptionParser()
parser.add_option("-t", "--task", help="""task code:\n\t\t\t\t\t
=gk_spectrum for G(k) spectrum building.\n\r
\n\t\t\t\t\tRequires options:
\n\t\t\t\t\t-p : proteinortho file [myproject.proteinortho]
\n\t\t\t\t\t-t : table with singletones
\n\t\t\t\t\t-o : outfile name [in txt or tsv format]
\t\t\t\t\t=gk_spectrum2 for G(k) spectrum building of a subset corresponding to a peak.\n
\n\t\t\t\t\tRequires options:
\n\t\t\t\t\t-p : proteinortho file [myproject.proteinortho]
\n\t\t\t\t\t-s : directory with sequences in fasta format
\n\t\t\t\t\t-o : outfile name [in txt or tsv format]
\n\t\t\t\t\t-k : peak position
\t\t\t\t\t=peak_info   for information about phyletic patterns distribution and double peak test\n
\n\t\t\t\t\tRequires options:
\n\t\t\t\t\t-p : proteinortho file [myproject.proteinortho]
\n\t\t\t\t\t-k : peak position
\t\t\t\t\t=peakgenes   for obtaining genes from the most popular phyletic pattern in a given peak\n
\n\t\t\t\t\tRequires options:
\n\t\t\t\t\t-p : proteinortho file [myproject.proteinortho]
\n\t\t\t\t\t-k : peak position
\n\t\t\t\t\t-s : directory with sequences in fasta format
\n\t\t\t\t\t-o : outfile name [in fasta format]
\t\t\t\t\t=tree_align  for building alignments of universal genes [for further tree building]\n
\t\t\t\t\t\nRequires options:
\t\t\t\t\t\n-p : proteinortho file [myproject.proteinortho]
\t\t\t\t\t\n-s : directory with sequences in fasta format
\t\t\t\t\t\n-a : directory with alignments
\t\t\t\t\t\n-f : fuse alignments or leave them be? (y/n)
\t\t\t\t\t\n-u : number of genes used in the fusion building or "all" if all
\t\t\t\t\t\n-o : outfile name [in fasta format]
\t\t\t\t\t\n-k : peak position [for tree marking]
\t\t\t\t\t\n-b : name of the blast database for the outgroup""", default="gk_spectrum")
parser.add_option("-p", "--proteinortho_file", default="./myproject.proteinortho")
parser.add_option("-g", "--genetables_direct", default="./")
parser.add_option("-o", "--outfile", default = "outfile")
parser.add_option("-k", "--peak_pos", default = 0)
parser.add_option("-s", "--seq_dir", default="")
parser.add_option("-a", "--al_dir_name", default="tree_alignment/")
parser.add_option("-f", "--fuse_aligns", default='y')
parser.add_option("-u", "--use_aligns", default="all")
parser.add_option("-b", "--blastdb", default=None)
parser.add_option("-x", "--singletone_table", default=None)
opt, args = parser.parse_args()
print opt.task
if opt.task == "gk_spectrum":
        print "hello!"
	gk_spectrum(opt.proteinortho_file, opt.singletone_table, opt.outfile)
elif opt.task == "gk_spectrum2":
	gk_spectrum2(opt.proteinortho_file, opt.seq_dir, eval(opt.peak_pos), opt.outfile)
elif opt.task == "peak_info":
	peak_info(eval(opt.peak_pos), opt.proteinortho_file)
elif opt.task == "peakgenes":
	peakgenes(eval(opt.peak_pos), opt.proteinortho_file, opt.seq_dir, opt.outfile)
elif opt.task == "tree_align":
	if opt.use_aligns != "all":
		use_aligns = eval(opt.use_aligns)
	else:
		use_aligns = opt.use_aligns
	if opt.fuse_aligns == 'y':
		fuse_aligns = True
	else:
		fuse_aligns = False
	tree_alignment(opt.proteinortho_file, opt.seq_dir, opt.al_dir_name, fuse_aligns, use_aligns, opt.outfile, eval(opt.peak_pos), opt.blastdb)


#test: 

"""
def find_phyletic_pattern(peak_pos, proteinortho_file):
	phylpat_dict = dict()
	for q in open(proteinortho_file):
		if q[0] == '#':
			continue
		q = q.strip().split()
		if eval(q[0]) == peak_pos:
			for p in q[3:]:
				spec = p.split('|')[0]
				try:
					phylpat_dict[spec] += 1
				except:
					phylpat_dict[spec] = 1
	peak_vals = sorted(phylpat_dict.values())[-peak_pos:]
	peak_spcs = []
	for k in phylpat_dict.keys():
		if phylpat_dict[k] in peak_vals:
			peak_spcs.append(k)
	return peak_spcs, list(phylpat_dict.keys())

def reverse_pattern_test(peak_pos, proteinortho_file):
	phylpat, all_spcs = find_phyletic_pattern(peak_pos, proteinortho_file)
	reverse_phylpat = [q for q in all_spcs if q not in phylpat]
	corresp_phylpat, all_spcs = find_phyletic_pattern(len(all_spcs) - peak_pos, proteinortho_file)
	reverse_phylpat = sorted(reverse_phylpat)
	corresp_phylpat = sorted(corresp_phylpat)
	for i in range(len(reverse_phylpat)):
		if reverse_phylpat[i] != corresp_phylpat[i]:
			return False
	return True
"""
#-------+---+--------++---+--+-----++----+-
#+++-++++--++-+-+++-++++-+++++++-++++++++++


