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
include("constraints.py")  # Contains constraint functions to limit reaction network

# Import mod packages
from mod.causality import Simulator

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================

# Initial concentrations for starting molecules
FORMALDEHYDE_INIT = 100    # Initial counts of formaldehyde
GLYCOLALDEHYDE_INIT = 1000 # Initial counts of glycolaldehyde

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
    based on the rule objects.
    
    Args:
        hyperedge: The hyperedge representing the reaction
        
    Returns:
        tuple: (rate_constant, should_be_cached) where rate_constant is the rate
               and should_be_cached indicates if the rate is cached
    """
    rule_rates = [rates[rule] for rule in hyperedge.rules]
    return rule_rates[0], False

# Dictionary mapping reaction names to their rate constants
rates = {
    ALDOL_ADD_F: ALDOL_ADDITION_RATE,
    ALDOL_ADD_B: ALDOL_ADDITION_REVERSE_RATE,
    KETO_ENOL_F: KETO_ENOL_RATE,  
    KETO_ENOL_B: KETO_ENOL_REVERSE_RATE, 
}

# Initial state of the simulation
# This defines the starting concentrations of each molecule
init_state = {glycolaldehyde: GLYCOLALDEHYDE_INIT, formaldehyde: FORMALDEHYDE_INIT}
    
# Expansion strategy that applies constraints during network growth
# This ensures that only chemically reasonable molecules are generated
expansion_strategy = rightPredicate [lambda d: all_constraints_apply(CONSTRAINT_FUNCTIONS, d)] ( reaction_rules )

# =============================================================================
# SIMULATION EXECUTION
# =============================================================================

# Create a new simulator for each run
sim = Simulator(
     graphDatabase=[formaldehyde, glycolaldehyde],  # Starting molecules
     expandNetwork=Simulator.ExpandByStrategy(expansion_strategy),  # How to grow the network
     initialState=init_state,  # Initial counts
     draw=Simulator.DrawMassAction(reactionRate=reaction_rate)  # Kinetics model
)

# Run the simulation for 100 time units
trace = sim.simulate(time=100)

# Print the simulation trace to see what happened
trace.print()
