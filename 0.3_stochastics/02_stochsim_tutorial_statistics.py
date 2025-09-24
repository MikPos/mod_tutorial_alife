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

ls = LabelSettings(LabelType.Term, LabelRelation.Specialisation)

init_state = {glycolaldehyde: GLYCOLALDEHYDE_INIT, formaldehyde: FORMALDEHYDE_INIT}
    
expension_strategy = (rightPredicate [lambda d: all_constraints_apply(CONSTRAINT_FUNCTIONS, d)] ( reaction_rules ))

simulations = []

for index in range(50):
    print(f"Starting simulation {index+1}/50")
    sim = stoch.Simulator(
            labelSettings=ls,
            graphDatabase=[formaldehyde, glycolaldehyde],
            expandNetwork=stoch.ExpandByStrategy(expension_strategy),
            initialState=init_state,
            draw=stoch.DrawMassAction(reactionRate=reaction_rate)
    )
    # Apply callbacks with parameter information
    setCallbacks(sim, verbose=False)
        
    # Run simulation
    trace = sim.simulate(time=100)
    del sim
    simulations.append(simulation_cache)

statistics = SimulationStatistics(simulations)

# Perform statistical analysis and plot the results
statistics.statistical_analysis(output_dir='analysis_results',
    uncertainty_type='ci',
    confidence_level=0.95
)


