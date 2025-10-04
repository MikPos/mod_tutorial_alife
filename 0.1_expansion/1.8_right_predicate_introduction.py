#We can restrict our search with a right predicate.
#This restricts the strategy to only add reactions to the reaction network where the product satisfies some predicate.
#Intuitively, this prevents our strategy from adding reactions to our reaction network where the products do not satisy some constraint

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

#We might want to constrain some property of the molecules we produce.
#This pattern represents a path graph with 8 vertices with any vertex labels and edge labels.
#If we constrain our search to exclude any molecule with this pattern, we will{*}
#never produce a molecule with a chain of length 8
pattern = Graph.fromDFS("[*]{*}[*]{*}[*]{*}[*]{*}[*]{*}[*]{*}[*]{*}[*]")
pattern.print()

#We define our DG object and the initial graphs present in the reaction network.
#Additionally, we set the label setting to the "term" setting since we are using wild cards.
ls = LabelSettings(LabelType.Term, LabelRelation.Specialisation)
dg = DG(graphDatabase=input_molecules, labelSettings=ls)

#We make the dg active so we can modify it
with dg.build() as b:

   #We define a strategy as follows
   strategy = (

      addSubset(input_molecules) 

      #We restrict our strategy to only add reactions to the reaction network where the product do not contain some pattern
      >> rightPredicate[lambda derivation: all(pattern.monomorphism(g, labelSettings=ls)<1 for g in derivation.right)](
         repeat[10](inputRules))
      

   )

   #We execute the strategy on the active dg
   b.execute(strategy)

#We print the dg and the result can be seen in "summary.pdf"
dg.print()