# Takes in a regina iso sig of geometric triangulation and 
# outputs a drawing of the cusp
from sage.all import arg, pi, real, imag, I, polygon, point, line, show
import regina
import snappy
import matplotlib.pyplot as plt

###TODO
# pretty up code + comments
# test on noncusped, nonhyperbolic, consider number of cusps?
# check + error message for weird edge case (figure8)
# recongifure to be able to output many pics w/ names
# !actually parse the edge data in a sane way

pallet1 = ['coral', 'green', 'lightblue', 'gray', 'red', 'lightblue', 'pink', 'blueviolet', 'seashell', 'sienna', 'teal', 'thistle', 'peru']
pallet2 = ["#EF476F", "#FFD166", "#06D6A0", "#118AB2", "#073B4C"]
purply = ["#dabfff", "#907ad6", "#4f518c", "#2c2a4a", "#7fdeff", "#E06BC7", "#C329A2", "#621551", "#E06C8B", "#E0886C"]
pinky = ["#f72585","#b5179e","#7209b7","#560bad","#480ca8","#3a0ca3"]

pallet = purply

test1 = 'cPcbbbiht' 
test2 = 'eLAkaccddngbqw'
test3 = 'fLAMcbccdeehhnaqw' # figure 8 sister geometric
test4 = 'fLLQcacdedejbqqww' # isolated figure 8 1
test5 = 'fLLQccecddehqrwwn' # isolated figure 8 2
test6 = 'dLQbcccdegj' #failing!
ISO_SIG = 'dLQacccjgnb'

DRAW_VERTICES = False
DRAW_EDGES = False


def edge_of(d, shapes):
	tet = int(d.detail()[0]) #very dumb!! parse this properly
	v1 = int(d.detail()[3])
	v2 = int(d.detail()[4])
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

# edge: (tet, v1, v2, shape)
# triangulation_edges: [[edges]]
def organize(sig):
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

# given tet# and vertex#, returns list
# e.g., [((tet, v1, v, shape), 1), ((tet, v, v2, shape), 3), 
#          ((tet, v, v3, shape, 4))]
# unused?
def getVertexTriangle(tet, v):
	vertexData = []
	for c in range(len(triangulation_edges)):
		for e in triangulation_edges[c]:
			if e[0] == tet:
				if e[1] == v or e[2] == v:
					vertexData.append((e, c))
	return vertexData

def plot_complex_triangle(z1, z2, z3, color='blue', vertex_color='red'):
	vertices = [(z1.real(), z1.imag()), (z2.real(), z2.imag()), (z3.real(), z3.imag())]
	vertices.append(vertices[0])

	triangle = polygon(vertices, fill=True, color=color, alpha=1)
	if DRAW_VERTICES:
		triangle  += sum(point((z.real(), z.imag()), size=30, color='red') for z in [z1, z2, z3])
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
						" (v1,v2,v1 sequence) ")
					return True
	return False


# given two vertices, find the edge (which was not just visted = class+ht)
# in which the two vertices are adjacent
# return (edgeclass index, edge index, 1/-1 for orientation, ht)
# 1 ascending index, -1 descending index
def findSequence(t1, v1, t2, v2, last_edge_class, lht, last_edge_i, triangulation_edges):
	for edge_class_i in range(len(triangulation_edges)):
		edge_class = triangulation_edges[edge_class_i]
		for i in range(len(edge_class)):
			j = (i+1) % len(edge_class)
			k = (i+2) % len(edge_class)
			h = (i-1) % len(edge_class)
			(tt1, tv1, tt2, tv2) = (edge_class[i][0], edge_class[i][1], edge_class[j][0], edge_class[j][1]) #tail
			(ht1, hv1, ht2, hv2) = (edge_class[i][0], edge_class[i][2], edge_class[j][0], edge_class[j][2]) #head
			if (last_edge_i not in [i,j,k,h]) or (lht == 1) or (last_edge_class != edge_class_i): # tails
				if (tt1, tv1, tt2, tv2) == (t1, v1, t2, v2):
					return (edge_class_i, i, 1, 0)
				elif (tt1, tv1, tt2, tv2) == (t2, v2, t1, v1):
					return (edge_class_i, (i + 1) % len(edge_class) , -1, 0)
			if (last_edge_i not in [i,j,k,h]) or (lht == 0) or (last_edge_class != edge_class_i): # heads
				if (ht1, hv1, ht2, hv2) == (t1, v1, t2, v2):
					return (edge_class_i, i, 1, 1)
				elif (ht1, hv1, ht2, hv2) == (t2, v2, t1, v1):
					return (edge_class_i, (i + 1) % len(edge_class), -1, 1)

	print(f'Failed to find sequence: T{t1}V{v1}, T{t2}V{v2}; Last Edge: {last_edge_class} {'heads' if lht else 'tails'}, {last_edge_i}. Exiting...')
	exit()

# task = (origin offset : complex num, p1, (t1,v1), (t2,v2), edge just visited : (edgeclass index, head vs tail, edgei)))
# head = 1, tail = 0; head is edge[2], tail is edge[1]
### CLEANME
def draw(triangulation_edges, depth, output='output.png'):
	# First pass around origin
	TODO = []
	edge0 = triangulation_edges[0]
	pl = plot_complex_triangle(0+0*I, 1+0*I, edge0[0][3], pallet[edge0[0][0]])
	TODO.append((1+0*I, edge0[0][3], (edge0[0][0], edge0[0][1]), (edge0[-1][0], edge0[-1][1]), (0, 0, 0)))
	prod = edge0[0][3]
	for e in range(1, len(edge0)):
		pl += plot_complex_triangle(0+0*I, prod, edge0[e][3]*prod, pallet[edge0[e][0]])
		if e < 3:
			TODO.append((prod, edge0[e][3]*prod, (edge0[e][0], edge0[e][1]), (edge0[e-1][0], edge0[e-1][1]), (0, 0, e)))
		prod *= edge0[e][3]


	# Points around first pass
	# NOTE: prod always based at 0; adjust accordingly
	while len(TODO) > 0:
		offset, p1, (t1,v1), (t2,v2), (lastedgeclass, lht, lastedgei) = TODO[0]
		TODO = TODO[1:]
		edge_class_i, edge, orientation, ht = findSequence(t1, v1, t2, v2, lastedgeclass, lht, lastedgei, triangulation_edges)
		edge_class = triangulation_edges[edge_class_i]
		current_edge_degree = len(edge_class)
		pl += plot_complex_triangle(offset, p1, (p1-offset)*edge_class[edge][3] + offset, pallet[edge_class[edge][0]])
		prod = (p1 - offset)*edge_class[edge][3]
		for i in range(1, current_edge_degree):
			ce = (edge + (i * orientation)) % current_edge_degree        #current edge
			le = (edge + (i - 1) * orientation) % current_edge_degree    #last edge
			pl += plot_complex_triangle(offset, prod+offset, edge_class[ce][3]*prod + offset, pallet[edge_class[ce][0]])
			if depth > 0: #don't add anymore 
				TODO.append((prod + offset, edge_class[ce][3]*prod + offset, 
					(edge_class[ce][0], edge_class[ce][ht+1]),
			 		(edge_class[le][0], edge_class[le][ht+1]), 
			 		(edge_class_i, ht, ce)))
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
	if checkFailure(triangulation_edges, ISO_SIG): return
	pp(triangulation_edges)
	print("Drawing...")
	draw(triangulation_edges, 100, ISO_SIG+'.png')




main()