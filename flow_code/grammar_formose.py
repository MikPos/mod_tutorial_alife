
# define input molecules:
formaldehyde = smiles("C=O", name="Formaldehyde")
glycolaldehyde = smiles( "OCC=O", name="Glycolaldehyde")


# import rules from GML files:
ketoEnol_F = ruleGML("rules/keto_enol_univ.gml")
ketoEnol_B = ruleGML("rules/keto_enol_univ.gml", invert=True)
aldolAdd_F = ruleGML("rules/aldol_add_univ.gml")
aldolAdd_B = ruleGML("rules/aldol_add_univ.gml", invert=True)

