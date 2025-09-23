
# define input molecules:
formaldehyde = Graph.fromSMILES("C=O", name="Formaldehyde")
glycolaldehyde = Graph.fromSMILES( "OCC=O", name="Glycolaldehyde")


# import rules from GML files:
ketoEnol_F = Rule.fromGMLFile("rules/keto_enol_univ.gml")
ketoEnol_B = Rule.fromGMLFile("rules/keto_enol_univ.gml", invert=True)
aldolAdd_F = Rule.fromGMLFile("rules/aldol_add_univ.gml")
aldolAdd_B = Rule.fromGMLFile("rules/aldol_add_univ.gml", invert=True)

