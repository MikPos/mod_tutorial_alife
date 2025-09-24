### flow query on a loaded dg 
## the goal is to search for a specific path through the network

# include the grammar used to create the dg 
include("grammar_formose.py")


# define the name of the dg to load
dg_name = "dg_autocat"
# this dg has been created before for this purpose, the script to make it is called make_autocat_dg.py

# load the dg in question from the previous dg dump file
dg = DG.load(inputGraphs, inputRules, f"{dg_name}.dg") 


# explore the size of the dg
def display_size_dg(dg):
    print("============ Number vertices ============")
    print(dg.numVertices)
    print("============ No edges ============")
    print(dg.numEdges)

display_size_dg(dg)

### set up a flow model 
# the flow model is the model we create that consists of linear equations to feed into an ILP solver
# give the flow a name to keep an overview
flow_name = "flow_cycle"


## this is the actual flow model ##
flow = hyperflow.Model(dg, ilpSolver = "CPLEX") # specified to use the CPLEX ILP solver
# the default solver is Gurobi, if you don't have an academic license for that use CPLEX


# to find a pathway through the network, we need inflows (sources) and outflows (sinks)
# Define the input molecules for the network
flow.addSource(formaldehyde)
flow.addSource(glycolaldehyde)

# Define the output molecules of the network
flow.addSink(glycolaldehyde) 

# Specify restrictions and amounts of used molecules
flow.addConstraint(inFlow[formaldehyde] == 2)
flow.addConstraint(inFlow[glycolaldehyde] == 1)
flow.addConstraint(outFlow[glycolaldehyde] == 2)


# Disable "strange" misleading input/output flows:
flow.allowIOReversal = False

if False:
    #### autocatalysis
    flow.overallAutocatalysis.enable()


### Set objective function
flow.objectiveFunction = isEdgeUsed

# # you can try different objective functions and their outcome:
# flow.objectiveFunction = edgeFlow
# # you can combine and formulate any linear equation as an objective function:
# flow.objectiveFunction = isEdgeUsed * 1000 + edgeFlow 


##### run flow ####
# Find a solution:
flow.findSolutions()

# # you can find more than one solution:
# flow.findSolutions(maxNumSolutions = 5) 

# Show solution information in the terminal
flow.solutions.list()

# Print solutions
flow.solutions.print()

# save the flow as a file
flow.dump(f"{dg_name}_{flow_name}.flow")


