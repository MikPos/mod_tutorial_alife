"""
Formose Reaction Network Definition

This module defines the molecular species and reaction rules for the formose
reaction network. The formose reaction is a complex autocatalytic reaction
that can produce various sugars from formaldehyde, making it relevant for
understanding prebiotic chemistry and the origin of life.

The network includes:
- Initial molecules: formaldehyde and glycolaldehyde
- Reaction rules: keto-enol isomerization and aldol addition (both forward and reverse)
- GML (Graph Markup Language) definitions for the reaction rules

The formose reaction is particularly interesting because it can produce
sugars that are essential for life, such as ribose, which is a component of RNA.
"""

# =============================================================================
# MOLECULAR SPECIES DEFINITION
# =============================================================================

# Define the initial molecular species for the formose reaction
# These are the starting materials that will be present at the beginning of the simulation

# Formaldehyde (H2C=O) - the simplest aldehyde and key building block
formaldehyde = Graph.fromSMILES("C=O", name="Formaldehyde")

# Glycolaldehyde (HOCH2CHO) - a simple sugar that can participate in the formose reaction
glycolaldehyde = Graph.fromSMILES("OCC=O", name="Glycolaldehyde")

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

# =============================================================================
# REACTION RULES CREATION
# =============================================================================

# Create reaction rule objects from the GML definitions
# Each rule is created in both forward and reverse directions to allow
# the reaction to proceed in both directions (reversible reactions)

# Keto-enol isomerization (forward direction)
KETO_ENOL_F = Rule.fromGMLString(ketoEnolGML)

# Keto-enol isomerization (reverse direction)
# The invert=True parameter creates the reverse reaction
KETO_ENOL_B = Rule.fromGMLString(ketoEnolGML, invert=True, name="Keto-enol isomerization reverse")

# Aldol addition (forward direction)
ALDOL_ADD_F = Rule.fromGMLString(aldolAddGML)

# Aldol addition (reverse direction)
# The invert=True parameter creates the reverse reaction
ALDOL_ADD_B = Rule.fromGMLString(aldolAddGML, invert=True, name="Aldol Addition reverse")

# =============================================================================
# REACTION RULES CONFIGURATION
# =============================================================================

# List of all reaction rules to be used in the simulation
# These rules define the possible chemical transformations in the network
reaction_rules = [KETO_ENOL_F, KETO_ENOL_B, ALDOL_ADD_F, ALDOL_ADD_B]
