import random
#Instead of using the strategy framework, we can define a strategy manually with the apply function.
#The advantage of using apply is that we can make any strategy we want.
#In this dummy example we randomly expand the reaction network for a number of iterations.
#We import the "random" library and use it for making random decisions.

#String containing a GML rule representing keto-enol isomerization
ketoEnolGML = """rule [
   ruleID "Keto-enol isomerization" 
   left [
      edge [ source 1 target 4 label "-" ]
      edge [ source 1 target 2 label "-" ]
      edge [ source 2 target 3 label "=" ]
   ]   
   context [
      node [ id 1 label "C" ]
      node [ id 2 label "C" ]
      node [ id 3 label "O" ]
      node [ id 4 label "H" ]
   ]   
   right [
      edge [ source 1 target 2 label "=" ]
      edge [ source 2 target 3 label "-" ]
      edge [ source 3 target 4 label "-" ]
   ]   
]"""

#String containing a GML rule representing aldol Addition 
aldolAddGML = """rule [
   ruleID "Aldol Addition"
   left [
      edge [ source 1 target 2 label "=" ]
      edge [ source 2 target 3 label "-" ]
      edge [ source 3 target 4 label "-" ]
      edge [ source 5 target 6 label "=" ]
   ]
   context [
      node [ id 1 label "C" ]
      node [ id 2 label "C" ]
      node [ id 3 label "O" ]
      node [ id 4 label "H" ]
      node [ id 5 label "O" ]
      node [ id 6 label "C" ]
   ]
   right [
      edge [ source 1 target 2 label "-" ]
      edge [ source 2 target 3 label "=" ]
      edge [ source 5 target 6 label "-" ]
      edge [ source 4 target 5 label "-" ]
      edge [ source 6 target 1 label "-" ]
   ]
]"""

#We load the rules we will use for the random expansion
ketoEnol_F = Rule.fromGMLString(ketoEnolGML)
ketoEnol_B = Rule.fromGMLString(ketoEnolGML, invert=True)
aldolAdd_F = Rule.fromGMLString(aldolAddGML)
aldolAdd_B = Rule.fromGMLString(aldolAddGML, invert=True)

#We create a list of rules which we will apply for every iteration.
rules = [aldolAdd_B, aldolAdd_F, ketoEnol_F, ketoEnol_B]

#We load the initial graphs for our reaction network
formaldehyde = Graph.fromSMILES("C=O", name="Formaldehyde")
glycolaldehyde = Graph.fromSMILES("OCC=O", name="Glycolaldehyde")

#We define our reaction network as usual
dg = DG(graphDatabase=[formaldehyde, glycolaldehyde])

#We build the reaction network so it is active
with dg.build() as b:

	#This is the list which will contain all of the graphs we will randomly select from when making a rule application
	graphs = [formaldehyde, glycolaldehyde]

	i = 0
	#The number of iterations. This controls how many random rule applications we will make
	num_iterations = 100

	while i < num_iterations:

	    i += 1

	    #We select two graphs at random from our list of graphs
	    graph_1 = random.sample(graphs, 1)[0]
	    graph_2 = random.sample(graphs, 1)[0]

	    #Now we have two graphs that we will try to apply all rules to
	    reactants = [graph_1, graph_2]

	    #We iterate through the rules and try to apply them to the reactant graphs
	    for rule in rules:

	    	#We collect the hyper edges the rule application created
	       	reaction_edges = b.apply(reactants, rule, onlyProper=False)

	       	#We iterate through the new hyper edges and collect all of the target graphs
	        for edge in reaction_edges:
	           	for vertex in edge.targets:
	           		graphs.append(vertex.graph)

#We print the resulting reaction network, which can be seen in "summary.pdf"
dg.print()





