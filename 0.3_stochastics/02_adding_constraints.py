"""
Stochastic Simulation Tutorial - Adding Constraints

This script demonstrates how to run multiple stochastic simulations with constraints
to study the formose reaction network. It runs 50 independent simulations to
generate statistical data about the reaction dynamics.

The constraints help ensure that only chemically reasonable molecules are formed
during the simulation, preventing the generation of unrealistic or unstable
intermediates that would not occur in real chemistry.

Key features:
- Multiple independent simulations (50 runs)
- Constraint application to limit reaction network growth
- Statistical sampling of reaction dynamics
- Memory management to handle large numbers of simulations
"""

import sys
# Get path to mod installation
path_to_mod = "/home/talax/xtof/local/Mod/bin/mod"
# Add mod /lib folder to sys to import packages from there
sys.path.insert(0, path_to_mod + "/lib64")

# Include necessary modules for the simulation
include("formose.py")      # Contains reaction rules and initial molecules
include("constraints.py")  # Contains constraint functions to limit reaction network
include("analysis.py")     # Contains analysis functions for processing results

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
ALDOL_ADDITION_RATE = 0.1              # Rate for aldol addition (forward)
KETO_ENOL_RATE = 0.1                  # Rate for keto-enol isomerization (forward)
ALDOL_ADDITION_REVERSE_RATE = ALDOL_ADDITION_RATE / 2  # Reverse rate (slower)
KETO_ENOL_REVERSE_RATE = KETO_ENOL_RATE                # Reverse rate (same as forward)

# Simulation parameters
SIMULATION_TIME = 100                 # Duration of each simulation
NUMBER_OF_SIMULATIONS = 5             # Number of independent simulations to run

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

# Run multiple independent simulations for statistical analysis
# Each simulation starts with the same initial conditions but follows
# a different stochastic trajectory
for index in range(NUMBER_OF_SIMULATIONS):
    print(f"Starting simulation {index+1}/{NUMBER_OF_SIMULATIONS}")
    
    # Create a new simulator for each run
    sim = stoch.Simulator(
        labelSettings=ls,  # Label settings for constraint checking
        graphDatabase=[formaldehyde, glycolaldehyde],  # Starting molecules
        expandNetwork=stoch.ExpandByStrategy(expansion_strategy),  # How to grow the network
        initialState=init_state,  # Initial concentrations
        draw=stoch.DrawMassAction(reactionRate=reaction_rate)  # Kinetics model
        )
    
    # Run the simulation for the specified time
    trace = sim.simulate(time=SIMULATION_TIME)
    
    # Print the simulation trace to monitor progress
    trace.print()
    
    # Clean up memory to prevent accumulation
    del sim
    del trace

