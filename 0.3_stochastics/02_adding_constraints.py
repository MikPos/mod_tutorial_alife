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

# Include necessary modules for the simulation
include("formose.py")      # Contains reaction rules and initial molecules
include("constraints.py")  # Contains constraint functions to limit reaction network

# Import mod packages
import mod.stochsim as stoch

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================

# Initial concentrations for starting molecules
FORMALDEHYDE_INIT = 100    # Initial counts of formaldehyde
GLYCOLALDEHYDE_INIT = 1000 # Initial counts of glycolaldehyde

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
# This defines the starting counts of each molecule
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
        graphDatabase=[formaldehyde, glycolaldehyde],  # Starting molecules
        expandNetwork=stoch.ExpandByStrategy(expansion_strategy),  # How to grow the network
        initialState=init_state,  # Initial counts
        draw=stoch.DrawMassAction(reactionRate=reaction_rate)  # Kinetics model
        )
    
    # Run the simulation for the specified time
    trace = sim.simulate(time=SIMULATION_TIME)
    
    # Print the simulation trace to monitor progress
    trace.print()
