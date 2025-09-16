
if False:
   #creating a molecule from a SMILES string
   glycolaldehyde = Graph.fromSMILES("OCC=O", name="Glycolaldehyde")
   glycolaldehyde.print()

if False:
   #Creating a molecule from graphDFS
   glycolaldehyde = Graph.fromDFS("[O]([C]([C](=[O])[H])([H])[H])[H]", name="Glycolaldehyde")
   glycolaldehyde.print()

if False:
   #Creating a molecule of your choice from a SMILES string
   my_own_smiles_string = "CCC=CO"#Change this string to a SMILES string representing a molecule of your choice
   my_molecule = Graph.fromSMILES(my_own_smiles_string, name="My own molecule")
   my_molecule.print()



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

if False:
   #Load the rule for keto-enol isomerization (rules can also be loaded from a gml file)
   ketoEnol_F = Rule.fromGMLString(ketoEnolGML)

   #We can also print rules, similar to how we print molecules
   ketoEnol_F.print()


if False:
   #We can define a rule that is the inverse of a given rule
   ketoEnol_B = Rule.fromGMLString(ketoEnolGML, invert=True)
   ketoEnol_B.print()


if False:
   #Here we load all of the four rules we use for formose chemistry and put all four rules in a list
   aldolAdd_F = Rule.fromGMLString(aldolAddGML)
   aldolAdd_B = Rule.fromGMLString(aldolAddGML, invert=True)
   ketoEnol_F = Rule.fromGMLString(ketoEnolGML)
   ketoEnol_B = Rule.fromGMLString(ketoEnolGML, invert=True)

   formose_rules = [ketoEnol_F, ketoEnol_B, aldolAdd_F, aldolAdd_B]

   #We can loop through the rules and print them 
   for rule in formose_rules:
      rule.print()



if False:
   #Now we are ready to apply a rule to a molecule. First we have to set up an environment that keeps track of all of the molecules we produce.
   #Formally we call it a reaction networ or a derivation graph, DG for short.

   #We first load the rules
   aldolAdd_F = Rule.fromGMLString(aldolAddGML)
   aldolAdd_B = Rule.fromGMLString(aldolAddGML, invert=True)
   ketoEnol_F = Rule.fromGMLString(ketoEnolGML)
   ketoEnol_B = Rule.fromGMLString(ketoEnolGML, invert=True)

   #First we define our input molecules
   formaldehyde = Graph.fromSMILES("C=O", name="Formaldehyde")

   glycolaldehyde = Graph.fromSMILES("OCC=O", name="Glycolaldehyde")

   input_molecules = [formaldehyde, glycolaldehyde]

   #We define our DG object and the initial molecules present in the reaction network
   dg = DG(graphDatabase=input_molecules)
   reaction_network = dg.build()

   #Now we apply the forward version of keto-enol isomerization to glycoaldehyde
   reaction_network.apply([glycolaldehyde], ketoEnol_F)


   #After we are done we lock the DG so we can print it
   del reaction_network


   #We print the resulting reaction network
   dg.print()


if False:
   #Instead of using the apply function, we can also a strategy framework that applys the rules for us.
   #Now we define a new DG, with the same input molecules as before.
   aldolAdd_F = Rule.fromGMLString(aldolAddGML)
   aldolAdd_B = Rule.fromGMLString(aldolAddGML, invert=True)
   ketoEnol_F = Rule.fromGMLString(ketoEnolGML)
   ketoEnol_B = Rule.fromGMLString(ketoEnolGML, invert=True)

   #First we define our input molecules
   formaldehyde = Graph.fromSMILES("C=O", name="Formaldehyde")

   glycolaldehyde = Graph.fromSMILES("OCC=O", name="Glycolaldehyde")

   input_molecules = [formaldehyde, glycolaldehyde]

   dg = DG(graphDatabase=input_molecules)
   reaction_network = dg.build()

   #We define a strategy as follows
   strategy = (


      addSubset(input_molecules) 

      >> inputRules
      #>> inputRules
      #>> inputRules

   )

   reaction_network.execute(strategy)

   del reaction_network

   dg.print()



if False:
   aldolAdd_F = Rule.fromGMLString(aldolAddGML)
   aldolAdd_B = Rule.fromGMLString(aldolAddGML, invert=True)
   ketoEnol_F = Rule.fromGMLString(ketoEnolGML)
   ketoEnol_B = Rule.fromGMLString(ketoEnolGML, invert=True)

   #First we define our input molecules
   formaldehyde = Graph.fromSMILES("C=O", name="Formaldehyde")

   glycolaldehyde = Graph.fromSMILES("OCC=O", name="Glycolaldehyde")

   input_molecules = [formaldehyde, glycolaldehyde]

   #We can make even more complicated strategies
   dg = DG(graphDatabase=input_molecules)
   reaction_network = dg.build()

   #We define a strategy as follows
   strategy = (


      addSubset(input_molecules) 

      >> repeat[4](inputRules)

   )

   reaction_network.execute(strategy)

   del reaction_network

   dg.print()

   



if True:
   aldolAdd_F = Rule.fromGMLString(aldolAddGML)
   aldolAdd_B = Rule.fromGMLString(aldolAddGML, invert=True)
   ketoEnol_F = Rule.fromGMLString(ketoEnolGML)
   ketoEnol_B = Rule.fromGMLString(ketoEnolGML, invert=True)

   #First we define our input molecules
   formaldehyde = Graph.fromSMILES("C=O", name="Formaldehyde")

   glycolaldehyde = Graph.fromSMILES( "OCC=O", name="Glycolaldehyde")

   input_molecules = [formaldehyde, glycolaldehyde]

   #We might want to constrain some property of the molecules we produce.
   #This pattern represents a path graph with 8 vertices with any vertex labels and edge labels.
   #If we constrain our search to exclude any molecule with this pattern, we will{*}
   #never produce a molecule with a chain of length 8
   pattern = graphDFS("[*]{*}[*]{*}[*]{*}[*]{*}[*]{*}[*]{*}[*]{*}[*]")
   pattern.print()


   ls = LabelSettings(LabelType.Term, LabelRelation.Unification)
   dg = DG(graphDatabase=input_molecules, labelSettings=ls)
   reaction_network = dg.build()

   #We define a strategy as follows
   strategy = (


      addSubset(input_molecules) 

      
      >> rightPredicate[lambda derivation: all(pattern.monomorphism(g, labelSettings=ls)<1 for g in derivation.right)](
         repeat[10](inputRules)
      )

   )

   reaction_network.execute(strategy)

   del reaction_network

   dg.print()




if False:
   aldolAdd_F = Rule.fromGMLString(aldolAddGML)
   aldolAdd_B = Rule.fromGMLString(aldolAddGML, invert=True)
   ketoEnol_F = Rule.fromGMLString(ketoEnolGML)
   ketoEnol_B = Rule.fromGMLString(ketoEnolGML, invert=True)

   #First we define our input molecules
   formaldehyde = Graph.fromSMILES("C=O", name="Formaldehyde")

   glycolaldehyde = Graph.fromSMILES("OCC=O", name="Glycolaldehyde")

   input_molecules = [formaldehyde, glycolaldehyde]

   #We might want to constrain some property of the molecules we consume.
   #Here we only allow reactions where the reactants together have at most 4 oxygen atoms

   dg = DG(graphDatabase=input_molecules)
   reaction_network = dg.build()

   strategy = (

      addSubset(input_molecules)
      
      >> leftPredicate[lambda derivation: sum([g.vLabelCount("O") for g in derivation.left])<=4](
         repeat[10](inputRules)
      )

   )

   reaction_network.execute(strategy)

   del reaction_network

   dg.print()





import random
if True:
   
   ketoEnol_F = Rule.fromGMLString(ketoEnolGML)
   ketoEnol_B = Rule.fromGMLString(ketoEnolGML, invert=True)
   aldolAdd_F = Rule.fromGMLString(aldolAddGML)
   aldolAdd_B = Rule.fromGMLString(aldolAddGML, invert=True)
   
   rules = [aldolAdd_B, aldolAdd_F, ketoEnol_F, ketoEnol_B]
   


   formaldehyde = Graph.fromSMILES("C=O", name="Formaldehyde")

   glycolaldehyde = Graph.fromSMILES("OCC=O", name="Glycolaldehyde")


   dg = DG(graphDatabase=[formaldehyde, glycolaldehyde])

   reaction_network = dg.build()


   molecules = [formaldehyde, glycolaldehyde]


   i = 0
   
   num_iterations = 10

   while i < num_iterations:

         i += 1

         molecule1 = random.sample(molecules, 1)[0]
         molecule2 = random.sample(molecules, 1)[0]

         reactants = [molecule1, molecule2]


         for rule in rules:

            reaction_edges = reaction_network.apply(reactants, rule, onlyProper=False)

            for edge in reaction_edges:
               for molecule in edge.targets:
                  molecules.append(molecule.graph)



   del reaction_network

   dg.print()

