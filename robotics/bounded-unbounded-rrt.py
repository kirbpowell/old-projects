# Code Sample
#
# This code will take a hypothetical robot arm
# from it's "neutral" state to a given destination.
# The algorithm used in this case is
# Rapidly exploring Random Trees - in both the bounded
# and the unbounded case.

import matplotlib.pyplot as plt
from random import uniform
from math import sin, sqrt, pi

"""
Part A:

theta_t+1 = theta_t + d*theta_dot_t
theta_dot_t+1 = theta_dot + d(u-sin(theta))
Xt = (theta_t; theta_dot_t)
Xt+1 = (theta_t+1; theta_dot_t+1)
"""

class Node():
	"""Class for the Nodes of the RRT. They basically form a
	reverse-linked-list, wherein each Node knows only about
	its own data, and its immediate parent. The Root
	(start point) of the tree has a parent of None.
	"""

	def __init__(self, data, parent=None):
		self.data = data
		self.parent = parent

"""GLOBAL VARIABLES"""

timestep = .05
maxTorque = .5
start = Node(((pi/2), 0), None)
goal = Node((0,0), None)

class RRT:
	"""Class for the RRT data structure. Its essentially
	a list made up of nodes of the Node class shown
	above. Does not enforce de-duplication natively,
	and for proper behavior this must be implemented
	separately.
	"""

	def __init__(self):
		self.nodeList = []

	def insertNode(self, node):
		self.nodeList.append(node)

	def addToTree(self, node, parent):
		toAdd = Node(node, parent)
		self.insertNode(toAdd)

def getXandY(nodeList):
    """Takes a list of tuples of the form (x, y) and
    returns all X values in an array, and all
    values in an array.

    Used to format output for graphing.
    """

    xArray = []
    yArray = []

    for node in nodeList:
        if isinstance(node, Node):
            node = node.data

        xArray.append(node[0])
        yArray.append(node[1])

    return xArray, yArray

def existsInTree(tree, sample):
	"""Compares the given node with each element
	in the tree, and returns TRUE iff a match is
	found and FALSE otherwise.
	"""

	for node in tree.nodeList:
		if sample.data == node.data:
			return True

	return False

def randomSample(tree):
	"""Sample from the configuration space where theta must be -PI >= theta >= PI
	and where theta2 must be -2 >= theta >= 2. If the sample already exists in the
	given tree, sample again until we find a samplewe haven't already explored.
	"""

	# Take an initial sample
	sample = Node((uniform(-pi, pi), uniform(-2, 2)))

	while existsInTree(tree, sample): # sample again until we haven't see said sample
		sample = Node((uniform(-pi, pi), uniform(-2, 2)))

	return sample

def distance(p1, p2):
	"""The distance metric calculates the euclidean distance
	between the points (theta, theta_dot) in configuration space
	"""

	return sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

def nearestToSample(tree, sample):
	"""This method uses the above distance metric to find
	the node in the tree that is nearest to the random sample
	and then returns that node.
	"""

	distances = {}

	for node in tree.nodeList:
		distances[node] = distance(node.data, sample.data)

	return min(distances, key=distances.get)

def getNewNode(nearestNode, sampleNode, bounded):
	"""Calculates a new node so that the new node is
	feasible, within torque limits (if applicable) and
	possible within the gien timestep allotment.

	Returns the new node's configuration (theta, thetaDot)
	"""

	thetaPrime = nearestNode.data[0] + nearestNode.data[1]*timestep
	u = ((sampleNode.data[1] - nearestNode.data[1])/timestep) - sin(nearestNode.data[0])

	if bounded:
		if u > .5:
			u = .5
		elif u < -.5:
			u = -.5

	thetaDotPrime = nearestNode.data[1] + timestep*(u - sin(nearestNode.data[0]))

	if thetaPrime > 2*pi:
		thetaPrime -= 2*pi

	return (thetaPrime, thetaDotPrime)

def linkTree(sample, tree, bounded):
	"""This function is the meat of the algorithm. It moves
	the robot arm  towards the random sample iff that new
	configuration is valid and unique. It will also add this
	new node to the tree and link it to the node nearest to
	the random sample.

	Additionally, this method will draw a line between the
	new node and its parent on the graph, as well as keeping
	track of those rejected new nodes that were rejected
	due to them colliding with the obstacle.
	"""

	found = False # tracks whether or not the goal has been reached. RETURNED
	nearestNode = nearestToSample(tree, sample) # get the node nearest to the sample

	newNode = getNewNode(nearestNode, sample, bounded)

	if not existsInTree(tree, Node(newNode)): #if the new node doesn't already exist in the tree, continue
		tree.addToTree(newNode, nearestNode) # add the node to the tree

		# plot it on the graph
		plt.plot([nearestNode.data[0], newNode[0]], [nearestNode.data[1], newNode[1]], 'b.-', animated=True)

		nearestGoalDist = distance(goal.data, newNode)

		if nearestGoalDist <= .1: # If we're close enough, call it good

			found = True
			plt.plot(newNode[0],newNode[1], 'go')

		return found
	else: # failed to add due to new node already existing in tree
		return found

def doRRT(bounded):
	"""This method contains the actaul RRT loop. It Kicks
	off and continues the RRT process until the goal state has
	been found.

	Returns the tree that contains the goal as its final node
	"""

	# Initialize the tree, the counter, and the stop condition
	tree = RRT()
	tree.insertNode(start)
	i = 0
	found = False

	# Plot the start point on the graph
	plt.plot(start.data[0], start.data[1], 'go')

	# The UNBOUNDED RRT loop
	while not found:
		sample = randomSample(tree)
		found = linkTree(sample, tree, bounded)

		if i % 500 == 0: # progress messages
			print "iteration: ",i

		i += 1

	if bounded:
		typeOfRRT = " Bounded"
	else:
		typeOfRRT = " Unbounded"

	print "*** RRT",typeOfRRT," Complete: Path from start to goal achieved. ***"
	print 'RRT took ',i,' iterations to produce a tree containing',len(tree.nodeList),'nodes.'

	return tree

def plotPathShowGraph(pathNode, bounded):
	"""Plot the path from start to finish (shown with larger, green colored dots)
	and also title the graph, and show it."""

	# Highlight the path from start to finish using bigger dots on the graph.
	while pathNode.parent is not None:
		plt.plot(pathNode.data[0], pathNode.data[1], pathNode.parent.data[0], pathNode.parent.data[1], 'go-', lw=10)
		pathNode = pathNode.parent

	# Add titles to the graph, and show it
	if bounded:
		plt.title("RRT for Inverted Pendulum\nBounded Torques")
	else:
		plt.title("RRT for Inverted Pendulum\nUnbounded Torques")

	plt.ylabel('Theta Dot')
	plt.xlabel('Theta')
	plt.show()

def main():
	"""The main method. Kicks off the RRT process,
	calls the linkTree method, and graphs final results.
	"""

	unbounded = False
	bounded = True

	unboundedTree = doRRT(unbounded)
	pathNode = unboundedTree.nodeList[-1]
	plotPathShowGraph(pathNode, unbounded)

	boundedTree = doRRT(bounded)
	bpathNode = boundedTree.nodeList[-1]
	plotPathShowGraph(bpathNode, bounded)

if __name__ == '__main__':
	main()
