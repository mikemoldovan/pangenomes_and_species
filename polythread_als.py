"""
polythread alignment building
Input:
1. Alignment_dir
2. N-treads
"""

from os import listdir, system
from optparse import OptionParser

def makefilelists(n_threads, aldir):
	filelists = []
	count = 0
	for q in listdir(aldir):
		if ".raw" in q:
			filename = aldir + q
			try:
				filelists[count%n_threads].append(filename)
			except:
				filelists.append([filename])
			count += 1
	print len(filelists)
	return filelists

def maketempsh(files, num):
	tempsh = open("tempsh_"+str(num)+".sh", "w")
	tempsh.write("#!/bin/bash\n")
	for rawfile in files:
		nonrawfile = '.'.join(rawfile.split('.')[:-1])
		tempsh.write("muscle -in "+rawfile+" -out "+nonrawfile+'\n')
	tempsh.close()

def main(n_threads, aldir):
	filelists = makefilelists(n_threads, aldir)
	albuild_sh = open("albuild_sh.sh", "w")
	albuild_sh.write("#!/bin/bash\n")
	for i in range(len(filelists)):
		maketempsh(filelists[i], i)
		albuild_sh.write("nohup	./"+"tempsh_"+str(i)+".sh	&\n")
	system("chmod +x tempsh_*.sh")
	system("chmod +x albuild_sh.sh")
	albuild_sh.close()


parser = OptionParser()
parser.add_option("-a", "--alignment_dir", help = "Directory with *.raw alignments")
parser.add_option("-n", "--n_scrs", help = "number of parallel scrs. Outfile: albuild_sh.sh")
opt, args = parser.parse_args()

main(eval(opt.n_scrs), opt.alignment_dir)
