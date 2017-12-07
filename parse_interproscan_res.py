"""
Interproscan result parser
Input:
interproscan file
Output:
standard GO-stat input
"""

import re
from optparse import OptionParser

def format1_print(go_list, outfile):
	for q in go_list:
		outfile.write(q+'\n')

def format2_print(gene_id, go_list, outfile):
	go_list_new = [q[3:] for q in go_list]
	outfile.write(gene_id + ' ' + ';'.join(go_list_new) + '\n')

def interpro_parse(interpro_file_name, output_file_name, outfmt):
	gene_id = ""
	go_list = []
	outfile = open(output_file_name, 'w')
	count = 1
	for q in open(interpro_file_name):
		q_list = q.split()
		if not gene_id:
			gene_id = q_list[0]
		elif q_list[0] != gene_id:
			go_list = list(set(go_list))
			if outfmt == "1" and go_list:
				format1_print(go_list, outfile)
			elif outfmt == "2" and go_list:
				format2_print(str(count), go_list, outfile)
			go_list = []
			gene_id = q_list[0]
			count += 1
		go_list += re.findall(r'GO:\d+', q)
	outfile.close()

parser = OptionParser()
parser.add_option("-i", "--input_file", help="Interproscan output file", default="myproject.proteinortho-graph")
parser.add_option("-o", "--output_file", help="Output file name", default="myproject.proteinortho-graph")
parser.add_option("-f", "--outfmt", help="Output format: 1 for GO list; 2 for GO database format (default 2)", default='2')
opt, args = parser.parse_args()

interpro_parse(opt.input_file, opt.output_file, opt.outfmt)