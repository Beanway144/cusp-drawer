# Takes in a regina iso sig of geometric triangulation and 
# outputs a drawing of the cusp
from sage.all import arg, pi, real, imag, I, polygon, point, line, show
import regina
import snappy
import matplotlib.pyplot as plt

###TODO
# use python classes instead of nested lists
# test on noncusped, nonhyperbolic, consider number of cusps?
# check + error message for weird edge case (figure8)
# recongifure to be able to output many pics w/ names

pallet1 = ['coral', 'green', 'lightblue', 'gray', 'red', 'lightblue', 'pink', 'blueviolet', 'seashell', 'sienna', 'teal', 'thistle', 'peru']
pallet2 = ["#EF476F", "#FFD166", "#06D6A0", "#118AB2", "#073B4C"]
purply = ["#dabfff", "#907ad6", "#4f518c", "#2c2a4a", "#7fdeff", "#E06BC7", "#C329A2", "#621551", "#E06C8B", "#E0886C"]
pinky = ["#f72585","#b5179e","#7209b7","#560bad","#480ca8","#3a0ca3"]

test1 = 'cPcbbbiht' 
test2 = 'eLAkaccddngbqw'
test3 = 'fLAMcbccdeehhnaqw' # figure 8 sister geometric
test4 = 'fLLQcacdedejbqqww' # isolated figure 8 1
test5 = 'fLLQccecddehqrwwn' # isolated figure 8 2
test6 = 'dLQbcccdegj' #failing!


##### PARAMETERS  #####
ISO_SIG = 'fLLQccecddehqrwwn'
pallet = purply
DRAW_VERTICES = False
DRAW_EDGES = False

#######################


def edge_of(d, shapes):
	"""
	Given a Regina edge and a list of the triangulation's shapes,
	parses the edge into a list 
	(tetrahedron index, tail vertex, head vertex, edge parameter).
	"""
	s_tet, s_edges = d.detail().split("(")
	tet = int(s_tet)
	v1, v2 = map(int, s_edges.split(")\n")[0])
	z = shapes[tet]
	param = 'fail'
	match (v1, v2):
		case (0, 1) | (1, 0) | (2 , 3) | (3, 2):
			param = z
		case (0, 2) | (2, 0) | (1, 3) | (3, 1):
			param = 1 / (1 - z)
		case (0, 3) | (3, 0) | (1, 2) | (2, 1):
			param = (z - 1) / z
	return (tet, v1, v2, param)


def organize(sig):
	"""
	Organizes edge data into a table given isosig of a triangulation.

	Output table is a doubly-nested list; each element of the output 
	is a list representing an edge class. Each element of an edge class 
	is an edge (ordered properly) taking the form 
	[tetrahedron index, tail vertex, head vertex, shape parameter].
	"""
	triangulation_edges = []
	T = regina.Triangulation3.fromIsoSig(sig)
	T.orient()
	M = snappy.Manifold(sig)
	shapes = M.tetrahedra_shapes(part='rect')
	for e in T.edges():
		edge = [] #equivalence class of edges
		for d in e:
			edge.append(edge_of(d, shapes))
		triangulation_edges.append(edge)
	return triangulation_edges

def plot_complex_triangle(z1, z2, z3, color='blue', vertex_color1='red', vertex_color2='red', vertex_color3='red'):
	vertices = [(z1.real(), z1.imag()), (z2.real(), z2.imag()), (z3.real(), z3.imag())]
	vertices.append(vertices[0])

	triangle = polygon(vertices, fill=True, color=color, alpha=1)
	if DRAW_VERTICES:
		# triangle  += sum(point((z.real(), z.imag()), size=30, color='red') for z in [z1, z2, z3])
		triangle += point((z1.real(), z1.imag()), size=30, color=vertex_color1)
		triangle += point((z2.real(), z2.imag()), size=30, color=vertex_color2)
		triangle += point((z3.real(), z3.imag()), size=30, color=vertex_color3)

	if DRAW_EDGES:
		edge_lines = line(vertices, color='white', thickness=0.5)
		triangle += edge_lines

	return triangle

# Current implementation fails and gives bad output on triangulations that have the pattern
# around an edge v1, v2, v1. So we just give up in these cases...(how often does this happen?)
# returns true if failure reached, else false
def checkFailure(triangulation_edges, sig):
	for edge_class_i in range(len(triangulation_edges)):
		edge_class = triangulation_edges[edge_class_i]
		for i in range(len(edge_class)):
			k = (i+2) % len(edge_class)
			if edge_class[i][0] == edge_class[k][0]:
				if (edge_class[i][1] == edge_class[k][1]) or (edge_class[i][2] == edge_class[k][2]):
					print(f"The triangulation {sig} has a particular structure that this method can't work on." +
						f" (v1,v2,v1 sequence at Edge Class {edge_class_i}, Edge {i})")
					return True
	return False


def checkEdgeCase(t1, v1, t2, v2, last_edge_class_i, lht, last_edge_i, last_orientation, triangulation_edges):
	edge_class = triangulation_edges[last_edge_class_i]
	last_edge = edge_class[last_edge_i]
	last_edge_minus_one = edge_class[(last_edge_i-1) % len(edge_class)]
	last_edge_plus_two = edge_class[(last_edge_i+2) % len(edge_class)]

	if last_orientation == -1:
		if (last_edge[0], last_edge[lht+1]) == (t1,v1): # ..., [v1], v2, ...
			if (last_edge_minus_one[0], last_edge_minus_one[lht+1]) == (t2, v2):
				print("-A")
				return (True, (last_edge_class_i, last_edge_i, 1, lht)) 
			if (last_edge_plus_two[0], last_edge_plus_two[lht+1]) == (t1, v1):
				print("-B")
				return (True, (last_edge_class_i, (last_edge_i + 2) % len(edge_class), 1, lht))

		if (last_edge[0], last_edge[lht+1]) == (t2,v2): # ..., [v2], v1, ...
			if (last_edge_minus_one[0], last_edge_minus_one[lht+1]) == (t1, v1):
				print("-C")
				return (True, (last_edge_class_i, last_edge_i, 1, lht)) 
			if (last_edge_plus_two[0], last_edge_plus_two[lht+1]) == (t2, v2):
				print("-D")
				return (True, (last_edge_class_i, (last_edge_i + 2) % len(edge_class), 1, lht))

	last_edge_minus_two = edge_class[(last_edge_i-2) % len(edge_class)]
	last_edge_plus_one = edge_class[(last_edge_i+1) % len(edge_class)]
	if last_orientation == 1:
		if (last_edge[0], last_edge[lht+1]) == (t1,v1): # ..., v2, [v1], ...
			if (last_edge_plus_one[0], last_edge_plus_one[lht+1]) == (t2, v2):
				print("A")
				return (True, (last_edge_class_i, last_edge_i, -1, lht)) 
			if (last_edge_minus_two[0], last_edge_minus_two[lht+1]) == (t1, v1):
				print("B")
				return (True, (last_edge_class_i, (last_edge_i - 1) % len(edge_class), -1, lht))

		if (last_edge[0], last_edge[lht+1]) == (t2,v2): # ..., v1, [v2], ...
			if (last_edge_plus_one[0], last_edge_plus_one[lht+1]) == (t1, v1):
				print("C")
				return (True, (last_edge_class_i, last_edge_i, -1, lht)) 
			if (last_edge_minus_two[0], last_edge_minus_two[lht+1]) == (t2, v2):
				print("D")
				return (True, (last_edge_class_i, (last_edge_i - 1) % len(edge_class), -1, lht))

	return (False, 0)


	# for edge_class_i in range(len(triangulation_edges)):
	# 	edge_class = triangulation_edges[edge_class_i]
	# 	for i in range(len(edge_class)):
	# 		j = (i+1) % len(edge_class)
	# 		k = (i+2) % len(edge_class)
	# 		if edge_class[i][0] == edge_class[k][0]:
	# 			if (edge_class[i][1] == edge_class[k][1]): #tails aba
	# 				if (edge_class[i][0] == t1) and (edge_class[j][0] == t2): # t1,t2,t1
	# 					if (edge_class[i][1] == v1) and (edge_class[j][1] == v2): #v1,v2,v1
	# 						if last_orientation == 1:
	# 							return (True, (edge_class_i, k, -1, 0))
	# 						if last_orientation == -1:
	# 							return (True, (edge_class_i, i, 1, 0))
	# 				if (edge_class[i][0] == t2) and (edge_class[j][0] == t1): # t2,t1,t2
	# 					if (edge_class[i][1] == v2) and (edge_class[j][1] == v1): #v2,v1,v2
	# 						if last_orientation == 1:
	# 							return (True, (edge_class_i, k, -1, 0))
	# 						if last_orientation == -1:
	# 							return (True, (edge_class_i, i, 1, 0))

	# 			if (edge_class[i][2] == edge_class[k][2]): #heads aba
	# 				if (edge_class[i][0] == t1) and (edge_class[j][0] == t2): # t1,t2,t1
	# 					if (edge_class[i][2] == v1) and (edge_class[j][2] == v2): #v1,v2,v1
	# 						if last_orientation == 1:
	# 							return (True, (edge_class_i, k, -1, 1))
	# 						if last_orientation == -1:
	# 							return (True, (edge_class_i, i, 1, 1))
	# 				if (edge_class[i][0] == t2) and (edge_class[j][0] == t1): # t2,t1,t2
	# 					if (edge_class[i][2] == v2) and (edge_class[j][2] == v1): #v2,v1,v2
	# 						if last_orientation == 1:
	# 							return (True, (edge_class_i, k, -1, 1))
	# 						if last_orientation == -1:
	# 							return (True, (edge_class_i, i, 1, 1))
	# return (False, 0)

def findSequence(t1, v1, t2, v2, last_edge_class_i, lht, last_edge_i, last_orientation, triangulation_edges):
	"""
	Looks through the edge table for the adjacency occurence of the input vertices, but not
	the given occurence. 

	Input: 
		t1, v1: tetrahedron and vertex indices for first vertex
		t2, v2: tetrahedron and vertex indices for second vertex
		last_edge_class_i: edge class index of previous occurence
		lht: whether the last occurence was in heads or tails
		last_edge_i: edge index of first vertex (wrt orientation) of last occurence
		last_orientation: previous orientation, 1 or -1
		triangulation_edges: edge table of triangulation

	Returns the tuple (edge class index, edge index, orientation, head or tail),
	where 1 is ascending index orientation and -1 is descending index orientation,
	and 0 is tails and 1 is heads.
	"""
	# handle edge case v,w,v
	(tf, ret) = checkEdgeCase(t1, v1, t2, v2, last_edge_class_i, lht, last_edge_i, last_orientation, triangulation_edges)
	if tf:
		oli, ol, oo, oht = ret
		print(f'Input:\n' +
			f'    > T{t1}V{v1} and T{t2}V{v2} at {last_edge_i} going {last_orientation}.\n' +
			f'Ret:\n' + f'	- Continue from {oli}:{ol}, going {oo}.')
		return ret

	for edge_class_i in range(len(triangulation_edges)):
		edge_class = triangulation_edges[edge_class_i]
		for i in range(len(edge_class)):
			j = (i+1) % len(edge_class)
			k = (i+2) % len(edge_class)
			h = (i-1) % len(edge_class)
			(tt1, tv1, tt2, tv2) = (edge_class[i][0], edge_class[i][1], edge_class[j][0], edge_class[j][1]) #tail
			(ht1, hv1, ht2, hv2) = (edge_class[i][0], edge_class[i][2], edge_class[j][0], edge_class[j][2]) #head
			if (last_edge_i not in [i,j,k,h]) or (lht == 1) or (last_edge_class_i != edge_class_i): # tails
				if (tt1, tv1, tt2, tv2) == (t1, v1, t2, v2):
					return (edge_class_i, i, 1, 0)
				elif (tt1, tv1, tt2, tv2) == (t2, v2, t1, v1):
					return (edge_class_i, (i + 1) % len(edge_class) , -1, 0)
			if (last_edge_i not in [i,j,k,h]) or (lht == 0) or (last_edge_class_i != edge_class_i): # heads
				if (ht1, hv1, ht2, hv2) == (t1, v1, t2, v2):
					return (edge_class_i, i, 1, 1)
				elif (ht1, hv1, ht2, hv2) == (t2, v2, t1, v1):
					return (edge_class_i, (i + 1) % len(edge_class), -1, 1)

	print(f'Failed to find sequence: T{t1}V{v1}, T{t2}V{v2}; Last Edge: {last_edge_class_i} {'heads' if lht else 'tails'}, {last_edge_i}. Exiting...')
	exit()

# task = (origin offset : complex num, p1, (t1,v1), (t2,v2), edge just visited : (edgeclass index, head vs tail, edgei)))
# head = 1, tail = 0; head is edge[2], tail is edge[1]
### CLEANME
def draw(triangulation_edges, depth, output='output.png'):
	"""
	DOCME
	"""
	# First pass around origin
	TODO = []
	edge0 = triangulation_edges[0]
	pl = plot_complex_triangle(0+0*I, 1+0*I, edge0[0][3], pallet[edge0[0][0]])
	TODO.append((1+0*I, edge0[0][3], (edge0[0][0], edge0[0][1]), (edge0[-1][0], edge0[-1][1]), (0, 0, 0, 1)))
	prod = edge0[0][3]
	for e in range(1, len(edge0)):
		pl += plot_complex_triangle(0+0*I, prod, edge0[e][3]*prod, pallet[edge0[e][0]])
		if e < 3:
			TODO.append((prod, edge0[e][3]*prod, (edge0[e][0], edge0[e][1]), (edge0[e-1][0], edge0[e-1][1]), (0, 0, e, 1)))
		prod *= edge0[e][3]


	# Points around first pass
	# NOTE: prod always based at 0; adjust accordingly
	while len(TODO) > 0:
		offset, p1, (t1,v1), (t2,v2), (lastedgeclass, lht, lastedgei, lastorientation) = TODO[0]
		TODO = TODO[1:]
		edge_class_i, edge, orientation, ht = findSequence(t1, v1, t2, v2, lastedgeclass, lht, lastedgei, lastorientation, triangulation_edges)
		edge_class = triangulation_edges[edge_class_i]
		current_edge_degree = len(edge_class)
		pl += plot_complex_triangle(offset, p1, (p1-offset)*edge_class[edge][3] + offset, pallet[edge_class[edge][0]])
		prod = (p1 - offset)*edge_class[edge][3]
		for i in range(1, current_edge_degree):
			ce = (edge + (i * orientation)) % current_edge_degree        #current edge
			le = (edge + (i - 1) * orientation) % current_edge_degree    #last edge
			pl += plot_complex_triangle(offset, prod+offset, edge_class[ce][3]*prod + offset, pallet[edge_class[ce][0]])
			if depth > 0: #don't add anymore 
				(TODO.append((prod + offset, edge_class[ce][3]*prod + offset, 
					(edge_class[ce][0], edge_class[ce][ht+1]),
			 		(edge_class[le][0], edge_class[le][ht+1]), 
			 		(edge_class_i, ht, ce, orientation))))
				depth -= 1
			prod *= edge_class[ce][3]
		

	pl.save(output, dpi=600, axes=False)

def pp(edges):
	print(f"Edges of Triangulation ({ISO_SIG}):")
	for e in range(len(edges)):
		print(">>> Edge " + str(e) + ": ")
		for i in edges[e]:
			print("      - Tetrahedron " + str(i[0]) + " with vertices " + str(i[1]) + " and " + str(i[2]) + " || Dihedral angle: " + str(round(float(arg(i[3]) * 180 / pi), 2)))

def main():
	print("Getting shapes...")
	triangulation_edges = organize(ISO_SIG)
	pp(triangulation_edges)
	# if checkFailure(triangulation_edges, ISO_SIG): return
	
	print("Drawing...")
	draw(triangulation_edges, 100, ISO_SIG+'.png')




main()