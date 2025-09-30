#Here we define rules describing the formose chemistry and load one and print it.

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

#Here we load the rule for keto-enol isomerization. Rules can also be loaded from a gml file.
ketoEnol_F = Rule.fromGMLString(ketoEnolGML, name="keto_enol_forward")

#Here we print the rule, which creates the summary and out directories. 
#"summary.pdf" will contain a depiction of the rule.
ketoEnol_F.print()

#We can also define the inverse of a rule and print it.
ketoEnol_B = Rule.fromGMLString(ketoEnolGML, invert=True, name="keto_enol_backward")
ketoEnol_B.print()
