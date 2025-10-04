
# define input molecules:
formaldehyde = Graph.fromSMILES("C=O", name="Formaldehyde")
glycolaldehyde = Graph.fromSMILES( "OCC=O", name="Glycolaldehyde")

inputGraphs = [formaldehyde, glycolaldehyde]

# define molecules for flow queries
example_mol = Graph.fromSMILES("C(CO)(CO)=O", name="example molecule")


# import rules from GML files:
ketoEnol_F = Rule.fromGMLFile("rules/keto_enol_univ.gml")
ketoEnol_B = Rule.fromGMLFile("rules/keto_enol_univ.gml", invert=True)
aldolAdd_F = Rule.fromGMLFile("rules/aldol_add_univ.gml")
aldolAdd_B = Rule.fromGMLFile("rules/aldol_add_univ.gml", invert=True)


