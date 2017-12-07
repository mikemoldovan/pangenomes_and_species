"""
Refine biult pangenome:
Delete strains highly similar to some other strain
"""

import matplotlib.pyplot as plt
from optparse import OptionParser
from os import listdir


class Vertex:
	def __init__(self, name, cont_dict, weight = 0):
		self.name = name
		self.cont_dict = cont_dict #Vertex : weight
		self.weight = weight
#		self.weight = sum([cont_dict[k] for k in cont_dict.keys()])


class Graph:
	def __init__(self, vertex_list):
		self.vertex_list= vertex_list

	def delete_vertex(self, del_vertex):
		self.vertex_list.remove(del_vertex)
		for name in self.vertex_dict.keys():
			v = self.vertex_dict[name]
			v.weight -= v.cont_dict[del_vertex]
			v.cont_dict.pop(del_vertex, None)


def build_weight_dist(graph, plot=True):
	l = len(graph.vertex_list)
	weights = []
	for i in range(l-1):
		vertex1 = graph.vertex_list[i]
		for j in range(i+1,l):
			vertex2 = graph.vertex_list[j]
			weights.append(vertex1.cont_dict[vertex2])
	if plot:
		plt.plot(list(range(l(l-1)/2)), sorted(weights), 'b.')
		plt.show()
	return sorted(weights)

#---Works with graphs obtained by rebuild_graph function
def refine(graph, vertex_weight_cutoff):
	if len(graph.vertex_list) <= 2:
		return graph
	min_weight = graph.vertex_list[0].cont_dict[graph.vertex_list[1]]
	min_pair = [graph.vertex_list[0],graph.vertex_list[1]]
	l = len(graph.vertex_list)
	for i in range(l-1):
		vertex1 = graph.vertex_list[i]
		for j in range(i+1,l):
			vertex2 = graph.vertex_list[j]
			if vertex1.cont_dict[vertex2] < min_weight:
				min_weight = vertex1.cont_dict[vertex2]
				min_pair = [vertex1, vertex2]

	if min_weight > vertex_weight_cutoff:
		return graph

	if min_pair[0] < min_pair[1]:
		graph.delete_vertex(min_pair[0])
	else:
		graph.delete_vertex(min_pair[1])
	refine(graph, vertex_weight_cutoff)

# ---Builds graph with weight=number of common OGGs
def build_graph(proteinortho_graph_file):
	gfile = open(proteinortho_graph_file)
	vertex_dict = dict()
	vertex_list = []
	for q in gfile:
		if q[0] == '#':
			continue
		q = q.split()
		name1 = q[0].split('|')[0]
		name2 = q[1].split('|')[0]
		if name1 not in vertex_dict.keys():
			vertex_dict[name1] = Vertex(name=name1, cont_dict=dict(), weight=0)
		if name2 not in vertex_dict.keys():
			vertex_dict[name2] = Vertex(name=name2, cont_dict=dict(), weight=0)
		try:
			vertex_dict[name1].cont_dict[vertex_dict[name2]] += 1
			vertex_dict[name2].cont_dict[vertex_dict[name1]] += 1
			vertex_dict[name1].weight += 1
			vertex_dict[name2].weight += 1
		except:
			vertex_dict[name1].cont_dict[vertex_dict[name2]] = 1
			vertex_dict[name2].cont_dict[vertex_dict[name1]] = 1
			vertex_dict[name1].weight += 1
			vertex_dict[name2].weight += 1
	for k in vertex_dict.keys():
		vertex_list.append(vertex_dict[k])
	return Graph(vertex_list)

# ---Rebuilds graph with weight=number of different OGGs
def make_genenumbers_dict(genetables_direct):
	genenumbers = dict()
	for file in listdir(genetables_direct):
		filename = file.split('.')
		if filename[-1] != "txt":
			continue
		if "genetable" not in filename[0]:
			continue
		spc_name = file[0].split('_')[0]
		genenum = -1
		for strin in open(genetables_direct + filename):
			genenum += 1
		genenumbers[spc_name] = genenum
	return genenumbers

def rebuild_graph(graph, genetables_direct):
	genenumbers = make_genenumbers_dict(genetables_direct)
	for vertex in graph.vertex_list:
		for k in vertex.cont_dict.keys():
			similar = vertex.cont_dict[k]
			different = (genenumbers[vertex.name] - similar) + (genenumbers[k.name] - similar)
			vertex.cont_dict[k] = different
	for vertex in graph.vertex_list:
		weight = sum(list(vertex.cont_dict.values()))
		vertex.weight = weight
	return graph
# ----

def write_weight_dist(weights_sorted, filename):
	outfile = open(filename, "w")
	outfile.write("number\tweight\n")
	for i in range(len(weights_sorted)):
		outfile.write(str(i)+'\t'+str(weights_sorted[i])+'\n')


parser = OptionParser()
parser.add_option("-g", "--proteinortho_graph", help="ProteinORTHO defaulf output", default="myproject.proteinortho-graph")
parser.add_option("-p", "--plot", help="Bin: build plot? y/n", default="n")
parser.add_option("-c", "--cutoff", help="Minimal acceptible edge weight", default=0)
parser.add_option("-w", "--weight_file", help="Resulting weight distribution file name", default="edge_weights.txt")
parser.add_option("-r", "--rebuild_graph", help="Bin: make weights as differences? y/n", default="y")
parser.add_option("-t", "--gene_tables", help="Directory with gene tables: if you have -r y", default='./')
opt, args = parser.parse_args()

graph = build_graph(opt.proteinortho_graph)

if opt.rebuild_graph == 'y':
	graph = rebuild_graph(graph, opt.gene_tables)

if opt.cutoff != 0 and opt.rebuild_graph == 'y':
	graph = refine(graph, opt.cutoff)

weights = build_weight_dist(graph, opt.plot == 'y')
write_weight_dist(weights, opt.weight_file)





