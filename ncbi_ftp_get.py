"""
get from NCBI FTP server
Input:
1. Bacterial NCBI directory 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/'
2. Bacteria list
3. Add_string ('latest_assembly_versions/')
"""

#init_string = "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Ralstonia_solanacearum/latest_assembly_versions/"

import urllib2
import os
import optparse
import random

def getstrainlist(init_string):
	hfo = urllib2.urlopen(init_string)
	html = hfo.read().split()
	strainlist = []
	for i in range(len(html)):
		if (i-8)%11 == 0:
			strainlist.append(html[i])
	return strainlist

def sample_reduce(strainlist, sample_volume):
	l = len(strainlist)
	if l <= sample_volume:
		return strainlist
	index_list = list(range(l))
	new_strainlist = []
	for i in range(sample_volume):
		strain = strainlist.pop(random.randint(0, l-i))
		new_strainlist.append(strain)
	return new_strainlist


def faa_get(strainlist, init_string, sample_volume, file_identifier):
	strainlist = sample_reduce(strainlist[:], sample_volume)
	for strain in strainlist:
		url_dir = init_string + strain + '/'
		hfo = urllib2.urlopen(url_dir)
		html = hfo.read().split('\r\n')
		for file in html:
			file = file.split()
			try:
				if file_identifier in file[8] and ("_from_" not in file[8] or "_from_" in file_identifier):
                                        print "wget\t"+url_dir + file[8]
#					os.system("wget\t"+url_dir + file[8])
					break
			except:
				pass

def download_bacs(ncbi_dir, baclist_file, add_string, sample_volume, file_identifier):
	baclist = [q.strip() for q in open(baclist_file)]
	for bac in baclist:
		init_string = ncbi_dir + bac + '/' + add_string
		try:
			os.mkdir(bac)
		except:
			pass
		os.chdir(bac)
		strainlist = getstrainlist(init_string)
		faa_get(strainlist, ncbi_dir + bac + '/' + add_string, sample_volume, file_identifier)
		os.chdir('..')

parser = optparse.OptionParser()
parser.add_option("-n", "--ncbi_dir", help="Directory on the NCBI FTP server ", default='ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/')
parser.add_option("-b", "--baclist", help="Bacterial species list", default="baclist.txt")
parser.add_option("-a", "--add_string", help="Subdirectory name", default='latest_assembly_versions/')
parser.add_option("-k", "--sample_volume", help="Number of downloaded genomes", default="all")
parser.add_option("-i", "--file_identifier", help="String that marks the file with needed data, def=.faa.gz", default = ".faa.gz")
opt, args = parser.parse_args()

download_bacs(opt.ncbi_dir, opt.baclist, opt.add_string, opt.sample_volume, opt.file_identifier)
#strainlist = getstrainlist(init_string)
#faa_get(strainlist)
