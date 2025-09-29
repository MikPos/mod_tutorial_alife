#Creating a graph from graph DFS and printing it

#Here we load the graph and name it. We can omit the name and it will automatically be named
glycolaldehyde = Graph.fromDFS("[O]([C]([C](=[O])[H])([H])[H])[H]", name="Glycolaldehyde")

#Here we print the loaded graph. More advanced printing options can be found in the mød documentation.
#The illustration of the graph can be found in the created directories named summary. 
#Usually we are only interested in the summary directory where "summary.pdf" can be found, but sometimes the out directory is also relevant.
glycolaldehyde.print() 