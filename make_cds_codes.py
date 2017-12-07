"""
Raname fastas
"""

import gzip
from os import listdir, chdir
from optparse import OptionParser

def make_id_dict(genetable_file):
	id_dict = dict()
	for q in open(genetable_file):
		if "init_ID" in q:
			continue
		q = q.strip().split()
		id_dict[q[0]] = q[1]
	return id_dict

def make_spec_dict(spectable):
	spec_dict = dict()
	for q in open(spectable):
		if "old_name" in q:
			continue
		q = q.strip().split()
		spec_dict[q[0][:-15]] = q[1]
	return spec_dict

def gzfasta_rewrite(gz_fafile, genetable_file, spec_dict):
	id_dict = make_id_dict(genetable_file)
	wd = '/'.join(gz_fafile.split('/')[:-1]) + '/'
	spc = wd + spec_dict[gz_fafile[:-24].split('/')[-1]]
	spc = open(spc + '_n.fa', "w")
	seq = ""
	_id = ""
	for q in gzip.open(gz_fafile):
		if q[0] == '>':
			if seq != "":
				for k in id_dict.keys():
					if k in _id:
						spc.write('>'+id_dict[k]+'\n')
						spc.write(seq+'\n')
						seq = ""
						break
			_id = q
		else:
			seq += q.strip()
	spc.close()

def rewrite_all(species_dir):
        specdicts = []
        for q in listdir(species_dir):
                if "spectable" in q:
                        specdicts.append(make_spec_dict(species_dir+q))
        spec_dict = dict()
        for d in specdicts:
                for k in d.keys():
                        spec_dict[k] = d[k]
	for q in listdir(species_dir):
		_file = species_dir + q
		if "_cds_from_genomic.fna.gz" not in _file:
			continue
		spc = q[:-24]
		genetable_file = species_dir + spec_dict[spc] + "_genetable.txt"
		gzfasta_rewrite(_file, genetable_file, spec_dict)

parser = OptionParser()
parser.add_option("-d", "--species_dir", help="Directory with everything about a species", default="./")
#parser.add_option("-s", "--sspectables", help="Directory with everything about a species", default="./")
opt, args = parser.parse_args()

print opt.species_dir
rewrite_all(opt.species_dir)
