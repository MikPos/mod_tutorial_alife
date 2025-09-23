# Import MOD packages
import mod.stochsim as stoch

# Include all necessary rules
## Recommend listing rules using 0 indexing (Python), as the 
rule0 = ruleGML()


# Include all necessary molecules
Molecules = include()

# Define input/starting reactants
Reactants = []

# List kinetic parameters as variables for each rule, these can be changed to your liking to edit the behaviour of the system.
k_r0 = 1.0  # Placeholder
k_r1 = 0.5  # Placeholder
km_r1 = 0.2 # Placeholder

# Defining the reaction rates for each rule as callbacks within an if statement
def reactionRate(r):
    # Rate of rule 0
    if rule0 in r.rules:
        # sim._marking is simply requesting the number of counts (MOD does not use concentrations) of the specified molecule
        if sim._marking[] >=0:  # Only applying when the reactant required for the rule is present
            r = k_r0            # This constant is automatically multipled by the current count of the species already specified in the previous line
            return r, False     # False refers to binning the value
        else:
            return 0, False     # Without the presence of the reactant the rate is of course 0
    # Alternatively we can define a specific rate equation that relies on the concentration of other species (e.g. a MM equation with inhibition)
    # Rate of rule 1
    if rule1 in r.rules:
        if sim._marking[] >= 0:
            r = k_r1 * (sim._marking[])     # r = k[A][B]
            return r, False
        else:
            return 0, False
    else:
        assert False 

# Stochastic Simulation
# Defining the initial state of the system and the simulation
seed = None
initialState = {}
sim = stoch.Simulator(
    labelSettings = LabelSettings(LabelType.Term, LabelRelation.Specialisation),
    graphDatabase = inputGraphs,
    expandNetwork = stoch.ExpandByStrategy(inputRules),
    initialState = initialState,
    draw = stoch.DrawMassAction(reactionRate=reactionRate)
)

# Simulate and draw the traces of each species, defining time in seconds
trace = sim.simulate(time=1000)
trace.print()
