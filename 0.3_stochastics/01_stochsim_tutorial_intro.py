"""
Stochastic Simulation Tutorial - Introduction

This script demonstrates the basic usage of MØD's stochastic simulation capabilities
for modeling the formose reaction network. It sets up a simple simulation with
formaldehyde and glycolaldehyde as starting materials and runs a single simulation
to show the basic workflow.

The formose reaction is a complex autocatalytic reaction that can produce various
sugars from formaldehyde, making it relevant for understanding prebiotic chemistry
and the origin of life.

Key components:
- Initial molecules: formaldehyde and glycolaldehyde
- Reaction rules: keto-enol isomerization and aldol addition (both forward and reverse)
- Constraints: prevent formation of certain forbidden subgraphs and limit molecule size
- Simulation: runs for 100 time units using mass-action kinetics
"""

import sys

# Include necessary modules for the simulation
include("formose.py")      # Contains reaction rules and initial molecules
include("callbacks.py")    # Contains callback functions for data collection
include("constraints.py")  # Contains constraint functions to limit reaction network

# Import mod packages
import mod.stochsim as stoch

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================

# Initial concentrations for starting molecules
FORMALDEHYDE_INIT = 100    # Initial concentration of formaldehyde
GLYCOLALDEHYDE_INIT = 1000 # Initial concentration of glycolaldehyde

# Reaction rate constants
# These determine how fast each reaction occurs in the simulation
ALDOL_ADDITION_RATE = 0.01        # Rate for aldol addition (forward)
KETO_ENOL_RATE = 0.1              # Rate for keto-enol isomerization (forward)
ALDOL_ADDITION_REVERSE_RATE = ALDOL_ADDITION_RATE / 2  # Reverse rate (slower)
KETO_ENOL_REVERSE_RATE = KETO_ENOL_RATE           # Reverse rate (same as forward)

# =============================================================================
# SIMULATION SETUP
# =============================================================================

def reaction_rate(hyperedge):
    """
    Calculate reaction rates for mass-action kinetics.
    
    This function is called by the stochastic simulator to determine the rate
    of each reaction. It extracts the rate constant from the rates dictionary
    based on the rule name.
    
    Args:
        hyperedge: The hyperedge representing the reaction
        
    Returns:
        tuple: (rate_constant, is_reversible) where rate_constant is the rate
               and is_reversible indicates if the reaction can go backwards
    """
    rule_rates = [rates[rule.name] for rule in hyperedge.rules]
    return rule_rates[0], False

# Label settings required for constraint checking
# This tells MØD how to match and compare molecular graphs
ls = LabelSettings(LabelType.Term, LabelRelation.Specialisation)

# Dictionary mapping reaction names to their rate constants
rates = {
    "Aldol Addition": ALDOL_ADDITION_RATE,
    "Aldol Addition reverse": ALDOL_ADDITION_REVERSE_RATE,
    "Keto-enol isomerization": KETO_ENOL_RATE,  
    "Keto-enol isomerization reverse": KETO_ENOL_REVERSE_RATE, 
}

# Initial state of the simulation
# This defines the starting concentrations of each molecule
init_state = {glycolaldehyde: GLYCOLALDEHYDE_INIT, formaldehyde: FORMALDEHYDE_INIT}
    
# Expansion strategy that applies constraints during network growth
# This ensures that only chemically reasonable molecules are generated
expansion_strategy = (rightPredicate [lambda d: all_constraints_apply(CONSTRAINT_FUNCTIONS, d)] ( reaction_rules ))

# =============================================================================
# SIMULATION EXECUTION
# =============================================================================

# Create the stochastic simulator
sim = stoch.Simulator(
        graphDatabase=[formaldehyde, glycolaldehyde],  # Starting molecules
        expandNetwork=stoch.ExpandByStrategy(expansion_strategy),  # How to grow the network
        initialState=init_state,  # Initial concentrations
        draw=stoch.DrawMassAction(reactionRate=reaction_rate)  # Kinetics model
)

# Run the simulation for 100 time units
trace = sim.simulate(time=100)

# Print the simulation trace to see what happened
trace.print()

# Clean up memory
del sim

