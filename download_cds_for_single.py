"""
get nucleotide data from NCBI FTP server
Input:
1. Bacterial NCBI directory 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/'
2. Bacteria name on the ncbi server
2.1. spectable
3. Add_string ('latest_assembly_versions/')
4. Target directory

gets files _cds_from_genomic.fna.gz -- 
FASTA format of the nucleotide sequences corresponding to all CDS 
features annotated on the assembly, based on the genome sequence. See 
the "Description of files" section below for details of the file format.
"""

import urllib2
import optparse
from os import listdir, chdir, system, getcwd

def getstrainlist(init_string):
	hfo = urllib2.urlopen(init_string)
	html = hfo.read().split()
	strainlist = []
	for i in range(len(html)):
		if (i-8)%11 == 0:
			strainlist.append(html[i])
	return strainlist

def getspecnames(spectable_name):
	specnames = []
	for q in open(spectable_name):
		if "old_name\tnew_name" in q:
			continue
		q = q.strip().split()
		name = q[0][:-15]
		specnames.append(name)
	return specnames

def cds_get(strainlist, specnames, init_string, file_identifier):
	for strain in strainlist:
		if strain not in specnames:
			continue
		url_dir = init_string + strain + '/'
		hfo = urllib2.urlopen(url_dir)
		html = hfo.read().split('\r\n')
		for file in html:
			file = file.split()
			try:
				if (file_identifier in file[8]) and ("_from_" not in file[8]):
					print file[8]
#					print "help", url_dir + file[8]
					system("wget\t"+url_dir + file[8])
					break
			except:
				return strain

def download_bacs(ncbi_dir, bacname_ncbi, add_string, target_dir, spectable, file_identifier):
	init_string = ncbi_dir + bacname_ncbi + '/' + add_string
	work_dir1 = getcwd()
	chdir(target_dir)
	strainlist = getstrainlist(init_string)
	specnames = getspecnames(spectable)
	a = cds_get(strainlist, specnames, init_string, file_identifier)
	chdir(work_dir1)


parser = optparse.OptionParser()
parser.add_option("-n", "--ncbi_dir", help="Directory on the NCBI FTP server ", default='ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/')
parser.add_option("-b", "--bacname_ncbi", help="Name of the bacterium on the ncbi server")
parser.add_option("-a", "--add_string", help="Subdirectory name", default='latest_assembly_versions/')
parser.add_option("-d", "--target_dir", help="Directory to which the sequences must be placed")
parser.add_option("-s", "--spectable", help="standard spectable")
parser.add_option("-i", "--file_identifier", help="String that marks the file with needed data, def=.faa.gz", default = ".faa.gz")
opt, args = parser.parse_args()

print "hello!!"
download_bacs(opt.ncbi_dir, opt.bacname_ncbi, opt.add_string, opt.target_dir, opt.spectable, opt.file_identifier)

