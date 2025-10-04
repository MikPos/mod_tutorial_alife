#Here we create a derivation graph by using the "apply" function. Then we print it.

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


#We load the rule we want to use
ketoEnol_F = Rule.fromGMLString(ketoEnolGML, name="keto_enol_forward")

#We load the initial graphs we want in our reaction network
formaldehyde = Graph.fromSMILES("C=O", name="Formaldehyde")
glycolaldehyde = Graph.fromSMILES("OCC=O", name="Glycolaldehyde")

#We create a list with our input graphs
input_molecules = [formaldehyde, glycolaldehyde]

#We define our DG object and the initial graphs present in the reaction network
dg = DG(graphDatabase=input_molecules)

#Here we build the dg making it active with the handle b.
with dg.build() as b:

   #Now we apply the forward version of keto-enol isomerization to glycoaldehyde.
   b.apply([glycolaldehyde], ketoEnol_F)

#We print the resulting reaction network, which can be seen in "summary.pdf".
dg.print()