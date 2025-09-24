import sys
# Get path to mod installation
path_to_mod = "/home/talax/xtof/local/Mod/bin/mod"
# Add mod /lib folder to sys to import packages from there
sys.path.insert(0, path_to_mod + "/lib64")
include("formose.py")
include("callbacks.py")
include("constraints.py")

# Import mod packages
import mod.stochsim as stoch

FORMALDEHYDE_INIT = 100
GLYCOLALDEHYDE_INIT = 1000

def reaction_rate(hyperedge):
    rule_rates = [rates[rule.name] for rule in hyperedge.rules]
    return rule_rates[0], False

# We need this for the constraints to work
ls = LabelSettings(LabelType.Term, LabelRelation.Specialisation)

aldol_addition_rate = 0.01
keto_enol_rate = 0.1
aldol_addition_reverse_rate = aldol_addition_rate / 2
keto_enol_reverse_rate = keto_enol_rate / 2

rates = {
    "Aldol Addition": aldol_addition_rate,
    "Aldol Addition reverse": aldol_addition_reverse_rate,
    "Keto-enol isomerization": keto_enol_rate,  
    "Keto-enol isomerization reverse": keto_enol_reverse_rate, 
}

init_state = {glycolaldehyde: GLYCOLALDEHYDE_INIT, formaldehyde: FORMALDEHYDE_INIT}
    
expansion_strategy = (rightPredicate [lambda d: all_constraints_apply(CONSTRAINT_FUNCTIONS, d)] ( reaction_rules ))

sim = stoch.Simulator(
        graphDatabase=[formaldehyde, glycolaldehyde],
        expandNetwork=stoch.ExpandByStrategy(expansion_strategy),
        initialState=init_state,
        draw=stoch.DrawMassAction(reactionRate=reaction_rate)
)

# Run simulation
trace = sim.simulate(time=100)
trace.print()
del sim

