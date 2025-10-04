#We can restrict our search with a left predicate.
#This restricts the a strategy to only apply rules to graphs that satisfy the predicate.
#Intuitively, this prevents our strategy from adding reactions to our reaction network where the reactants do not satisy some constraint

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

#We load the rules we want to use to expand the chemical space.
aldolAdd_F = Rule.fromGMLString(aldolAddGML)
aldolAdd_B = Rule.fromGMLString(aldolAddGML, invert=True)
ketoEnol_F = Rule.fromGMLString(ketoEnolGML)
ketoEnol_B = Rule.fromGMLString(ketoEnolGML, invert=True)
inputRules = [aldolAdd_B, aldolAdd_F, ketoEnol_F, ketoEnol_B]

#We define our input graphs.
formaldehyde = Graph.fromSMILES("C=O", name="Formaldehyde")
glycolaldehyde = Graph.fromSMILES("OCC=O", name="Glycolaldehyde")
input_molecules = [formaldehyde, glycolaldehyde]

#We define our DG object and the initial graphs present in the reaction network
dg = DG(graphDatabase=input_molecules)

#We make the dg active so we can modify it
with dg.build() as b:

   #We define a strategy as follows
   strategy = (

      addSubset(input_molecules) 

      #We restrict our strategy to only apply rules to graphs with at most 4 oxygen atoms
      >> leftPredicate[lambda derivation: sum([g.vLabelCount("O") for g in derivation.left])<=4](
         repeat[10](inputRules))
      

   )

   #We execute the strategy on the active dg
   b.execute(strategy)

#We print the dg and the result can be seen in "summary.pdf"
dg.print()