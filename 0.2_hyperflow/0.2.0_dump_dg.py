### Create a DG and dump the files.

include("formose.py")


strat = (
	addUniverse(formaldehyde)
	>> addSubset(glycolaldehyde)
	# Constrain the reactions:
	# No molecules with more than 20 atoms can be created.
	>> rightPredicate[
		lambda derivation: all(g.numVertices <= 20 for g in derivation.right)
	](
		# Iterate until nothing new is found.
		repeat(
			inputRules
		)
	)
)

dg = DG(graphDatabase=inputGraphs)
dg.build().execute(strat)
dg.print()



dg.dump("dg_autocat.dg")
