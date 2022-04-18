# vertex_nodes is list of 3 lists of integers containing the numbers
# of elements of ls that are vertex_nodes
# similar for edge_nodes interior_nodes is a list of integers for
# interior_nodes 
#

class DualBasis:
	def __init__( self, ls , vertex_nodes , edge_nodes , interior_nodes ):
		self.ls = ls
		self.vertex_nodes = vertex_nodes
		self.edge_nodes = edge_nodes
		self.interior_nodes = interior_nodes
	def getBasis( self ):
		return self.ls
	def getVertexNodeIDs( self ):
		return self.vertex_nodes
	def getEdgeNodeIDs( self ):
		return self.edge_nodes
	def getInteriorNodeIDs( self ):
		return self.interior_nodes


