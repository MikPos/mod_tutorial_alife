
# define input molecules:
formaldehyde = Graph.fromSMILES("C=O", name="Formaldehyde")
glycolaldehyde = Graph.fromSMILES( "OCC=O", name="Glycolaldehyde")

inputGraphs = [formaldehyde, glycolaldehyde]

# define molecules for flow queries
example_mol = Graph.fromSMILES("C(CO)(CO)=O", name="example molecule")


# import rules from GML files:
ketoEnol_F = ruleGML("rules/keto_enol_univ.gml")
ketoEnol_B = ruleGML("rules/keto_enol_univ.gml", invert=True)
aldolAdd_F = ruleGML("rules/aldol_add_univ.gml")
aldolAdd_B = ruleGML("rules/aldol_add_univ.gml", invert=True)


