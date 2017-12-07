"""
Parse interPROscan results and obtain information about maximal number of OGGs
input:
proteinortho_file
interpro_file
go_database
phylpat

output:
file with goterms
file with goterms descriptions from the database
"""

import re
from optparse import OptionParser

def interpro_file_parse(interpro_file):
	go_dict = dict()
	for q in open(interpro_file):
		key = q.split()[0]
		gos = re.findall(r"GO:[0-9]+", q)
		gos = list(set(gos))
		try:
			goarr = go_dict[key]
			for go in gos:
				if go not in goarr:
					go_dict[key].append(go)
		except:
			go_dict[key] = gos
#	print go_dict
	return go_dict

def proteinortho_file_parse(proteinortho_file, phylpat, go_dict):
	ogg_dict = dict()
	ogg_id_dict = dict()
	count = 0
	for q in open(proteinortho_file):
		if q[0] == '#':
			continue
		q = q.strip().split('\t')
		genarr = []
		for i in range(len(q[3:])):
			coef = i+3
			if q[coef] != '*' and phylpat[i] == '+':
				genarr.append(q[coef])
			elif q[coef] == '*' and phylpat[i] == '-':
				pass
			else:
				genarr = []
				break
		if not genarr:
			continue
		ogg_dict[count] = []
		for gen in genarr:
			try:
				ogg_id_dict[count].append(gen)
			except:
				ogg_id_dict[count] = [gen]
			try:
				gen_gos = go_dict[gen]
				for go in gen_gos:
					if go not in ogg_dict[count]:
						ogg_dict[count].append(go)
			except:
				pass
		count += 1
#	print ogg_dict
	return ogg_dict, ogg_id_dict


def ogg_dict_transpose(ogg_dict):
	newdict = dict()
	for k in ogg_dict.keys():
		go_list = ogg_dict[k]
		for go in go_list:
			try:
				newdict[go].append(k)
			except:
				newdict[go] = [k]
	return newdict


def godata_parse(godata_file, ogg_dict_transposed, ogg_id_dict):
	goterm_info = dict()
	term_arr = []
	for q in open(godata_file):
		q = q.strip()
		if q == '':
			go = term_arr[1].split()[1]
			try:
				oggs = ogg_dict_transposed[go]
				for ogg in oggs:
					try:
						goterm_info[ogg].append(term_arr)
					except:
						goterm_info[ogg] = [term_arr]
			except:
				pass
			term_arr = []
		else:
			term_arr.append(q)
	for ogg in goterm_info.keys():
		print "\n------------ ogg:", ogg, ",".join(ogg_id_dict[ogg]),"--------------"
		for term_info in goterm_info[ogg]:
			print ''
			for line in term_info:
				print line


parser = OptionParser()
parser.add_option("-i", "--interpro_file", help="Interproscan filtered file")
parser.add_option("-g", "--godata_file", help="GO database")
parser.add_option("-s", "--phylpat", help="phyletic pattern")
parser.add_option("-p", "--proteinortho_file", help="proteinortho output table")
opt, args = parser.parse_args()

go_dict = interpro_file_parse(opt.interpro_file)
ogg_dict, ogg_id_dict = proteinortho_file_parse(opt.proteinortho_file, opt.phylpat, go_dict)
ogg_dict_t = ogg_dict_transpose(ogg_dict)
godata_parse(opt.godata_file, ogg_dict_t, ogg_id_dict)

"""
def godata_parse(godata_file, ogg_dict_transposed):
	goterm_info = dict()
	wasterm = False
	term_arr = []
	for q in open(godata_file):
		q = q.strip()
		if q == "[Term]":
			if term_arr:
				for ogg in go_ogg_list:
					try:
						goterm_info[]
			wasterm = True
			term_arr.append(q)
		if not wasterm:
			continue
		if q[0:3] == "id:":
			go = q.strip().split()[1]
			try:
				go_ogg_list = ogg_dict_transposed[go]
			except:
				wasterm = False
				term_arr = []
		term_arr.append(q)
"""