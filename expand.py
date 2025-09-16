include("formose.py")

post.summarySection("Input Graphs")

for a in inputGraphs:

   a.print()

post.summarySection("Input Rules")

for a in inputRules:

   a.print()

# Reaction networks are expaned using a strategy:

strat = (

   # A molecule can be active or passive during evaluation.

   addUniverse(formaldehyde) # passive

   >> addSubset(glycolaldehyde) # active

   # Aach reaction must have a least 1 active educt.

   >> repeat [4] (inputRules
    )

)

# We call a reaction network a 'derivation graph'.

dg = DG(graphDatabase=inputGraphs)

dg.build().execute(strat)

# They can also be visualised.

dg.print()
