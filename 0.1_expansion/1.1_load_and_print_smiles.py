#Creating a graph from a SMILES string and printing it

#Here we load the graph and name it. We can omit the name and it will automatically be named
glycolaldehyde = Graph.fromSMILES("OCC=O", name="Glycolaldehyde")

#Here we print the loaded graph. More advanced printing options can be found in the mød documentation.
#The illustration of the graph can be found in the created directories named summary and out. 
#Usually we are only interested in the summary directory, where "summary.pdf" can be found, but sometimes the out directory is also relevant.
glycolaldehyde.print() 