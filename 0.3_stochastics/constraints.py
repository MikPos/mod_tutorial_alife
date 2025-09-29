"""
Constraint Functions for Stochastic Simulation

This module defines constraint functions that limit the growth of the reaction
network during stochastic simulation. These constraints ensure that only
chemically reasonable molecules are generated, preventing the formation of
unrealistic or unstable intermediates.

The constraints include:
- Forbidden subgraphs: prevent formation of certain structural patterns
- Size limits: limit the maximum size of molecules
- Radical limits: control the number of radical centers

These constraints are essential for maintaining chemical realism in the
simulation and preventing combinatorial explosion of the reaction network.
"""

# =============================================================================
# CONSTRAINT PARAMETERS
# =============================================================================

# Global constraint parameters
FORBIDDEN_SUBGRAPHS = []  # List of forbidden structural patterns
MAX_NUMBER_OF_VERTICES = 20  # Maximum number of atoms in a molecule
MAX_NUMBER_OF_RADICALS = 1    # Maximum number of radical centers

# =============================================================================
# FORBIDDEN SUBGRAPHS DEFINITION
# =============================================================================

# Define forbidden subgraphs that should not be formed during the simulation
# These represent chemically unreasonable or unstable structures

# 3-membered ring (cyclopropane-like structure)
# This is highly strained and unlikely to form in the formose reaction
CYC3 = Graph.fromDFS("[*]1{*}[*]{*}[*]{*}1")
FORBIDDEN_SUBGRAPHS.append(CYC3)

# 4-membered ring (cyclobutane-like structure)
# This is also strained and not typical for sugar chemistry
CYC4 = Graph.fromDFS("[*]1{*}[*]{*}[*]{*}[*]{*}1")
FORBIDDEN_SUBGRAPHS.append(CYC4)

# Conjugated double bonds (C=C-C=C)
# This pattern is not typical in sugar chemistry and can lead to
# unrealistic aromatic-like structures
DOUBLE_BOND = Graph.fromDFS("[*]1=[*]2=[*]3")
FORBIDDEN_SUBGRAPHS.append(DOUBLE_BOND)

# =============================================================================
# CONSTRAINT FUNCTIONS
# =============================================================================

def forbiddenSubgraphs(derivation):
    """
    Check if any forbidden subgraphs are present in the derivation products.
    
    This function prevents the formation of chemically unreasonable structures
    by checking if any of the forbidden subgraphs can be found in the products
    of the reaction derivation.
    
    Args:
        derivation: The reaction derivation to check
        
    Returns:
        bool: True if no forbidden subgraphs are found, False otherwise
    """
    return all(subgraph.monomorphism(g, labelSettings=ls) == 0
                      for subgraph in FORBIDDEN_SUBGRAPHS for g in derivation.right)

def has_at_most_x_vertices(derivation):
    """
    Check if all products have at most the maximum allowed number of vertices.
    
    This function limits the size of molecules to prevent combinatorial explosion
    and maintain computational tractability. It also reflects the fact that
    very large molecules are less likely to form in the formose reaction.
    
    Args:
        derivation: The reaction derivation to check
        
    Returns:
        bool: True if all products have at most MAX_NUMBER_OF_VERTICES, False otherwise
    """
    return all(g.numVertices <= MAX_NUMBER_OF_VERTICES for g in derivation.right)

def all_constraints_apply(constraint_functions, derivation):
    """
    Check if all constraint functions are satisfied for a given derivation.
    
    This is a helper function that applies all constraint functions to a
    derivation and returns True only if all constraints are satisfied.
    
    Args:
        constraint_functions: List of constraint functions to apply
        derivation: The reaction derivation to check
        
    Returns:
        bool: True if all constraints are satisfied, False otherwise
    """
    return all(function(derivation) for function in constraint_functions)

# =============================================================================
# CONSTRAINT CONFIGURATION
# =============================================================================

# List of all constraint functions to apply during network expansion
# These functions will be called for each potential reaction to determine
# whether it should be allowed in the network
CONSTRAINT_FUNCTIONS = [forbiddenSubgraphs, has_at_most_x_vertices]
