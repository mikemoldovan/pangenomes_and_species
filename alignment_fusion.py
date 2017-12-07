"""
alignment fusion
"""

from os import listdir, remove
from optparse import OptionParser
from Bio import AlignIO

def fastaparse(fafile):
	fafile_dict = dict()
	for q in fafile:
		q = q.strip().split()[0]
		if q[0] == '>':
			gene_name = q[1:].split('|')[0]
			fafile_dict[gene_name] = ""
		else:
			fafile_dict[gene_name] += q
	return fafile_dict

def cat_align_dicts(aldict1, aldict2):
	if not aldict1:
		return aldict2
	if not aldict2:
		return aldict1
	cat_dict = dict()
	first = True
	for k in aldict1.keys():
		l1 = len(aldict1[k])
		break
	for k in aldict2.keys():
		l2 = len(aldict2[k])
		break

	for k in aldict1.keys():
		if k in aldict2.keys():
			cat_dict[k] = aldict1[k] + aldict2[k]
		else:
			cat_dict[k] = aldict1[k] + '-'*l2
	for k in aldict2.keys():
		if k not in aldict1.keys():
			cat_dict[k] = '-'*l1 + aldict2[k]
	return cat_dict


def cataligns(align_dir, outfile):
	cat_dict = dict()
	for al in listdir(align_dir):
		if ".fa" not in al or ".raw" in al:
			continue
		fafile_dict = fastaparse(open(align_dir + al))
		cat_dict = cat_align_dicts(cat_dict, fafile_dict)
		remove(align_dir + al)

	for k in cat_dict.keys():
		outfile.write('>' + k + '\n')
		outfile.write(cat_dict[k] + '\n')


def make_phylip(align_name):
	newname = '.'.join(align_name.split('.')[:-1]) + ".phy"
	input_handle = open(align_name)
	output_handle = open(newname, "w")
	alignment = AlignIO.read(input_handle, "fasta")
	AlignIO.write(alignment, output_handle, "phylip")
	output_handle.close()
	input_handle.close()
	remove(align_name)


parser = OptionParser()
parser.add_option("-d", "--align_dir", help="directory with fasta alignments", default="./")
parser.add_option("-o", "--outfile", help="Name of the output file", default="catalign.fa")
parser.add_option("-p", "--make_phylip", help="output in the phylip format? (y/n)", default='y')
opt, args = parser.parse_args()

cataligns(opt.align_dir, open(opt.outfile, "w"))
if opt.make_phylip == 'y':
	make_phylip(opt.outfile)


