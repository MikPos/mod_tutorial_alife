# ALIFE tutorial

### ILP ###

# include the grammar used to create the dg 
include("grammar_formose.py")


# define the name of the dg to load
dg_name = "dg_ex_3"

## dg_ex_3 is very big, potentially useful for ILP

## dg_ex_4 is also interesting, little small but could be extended more


# load the dg in question from the previous dg dump file
dg = DG.load(inputGraphs, inputRules, f"{dg_name}.dg") 


# explore the size of the dg
def display_size_dg(dg):
    print("============ Number vertices ============")
    print(dg.numVertices)
    print("============ No edges ============")
    print(dg.numEdges)

display_size_dg(dg)


##### set up the flow model #####


## simple flow from A to B
if True:

    flow_name = "flow_1"
    # set up a flow model 
    flow = hyperflow.Model(dg) 
    # the default solver is Gurobi, if you don't have an academic license for that use CPLEX
    flow = hyperflow.Model(dg, ilpSolver = "CPLEX") # this is a different ILP solver
    

    # Define the input molecules for the network
    flow.addSource(formaldehyde)
    flow.addSource(glycolaldehyde)
    # Define the output molecules of the network
    flow.addSink(glycolaldehyde) # for a flow to happen, some molecules must be defined for in and outflow

    # Specify restrictions and amounts of used molecules
    flow.addConstraint(inFlow[formaldehyde] == 2)
    flow.addConstraint(inFlow[glycolaldehyde] == 1)
    flow.addConstraint(outFlow[glycolaldehyde] == 2)

    # Disable "strange" misleading input/output flows:
    flow.allowIOReversal = False

    # Set objective function
    flow.objectiveFunction = isEdgeUsed
    # # you can try different objective functions and their outcome:
    # flow.objectiveFunction = edgeFlow
    # # you can combine and formulate any linear equation as an objective function:
    # flow.objectiveFunction = isEdgeUsed * 1000 + edgeFlow 


    

## forbidden molecules
if False:

    # by setting the flow of a specific graph to zero we forbid it from being used:
    flow.addConstraint(vertex["OCC(=O)CO"] == 0)

    # We can also forbid patterns within molecules:
    for v in dg.vertices:
        if v.graph.vLabelCount("C") > 4:
            flow.addConstraint(vertex[v] == 0)

## autocatalytic cycle
if False:

    # define the inflow and outflow
    flow.addSource(formaldehyde)
    flow.addSource(glycolaldehyde)

    # what can flow out of the network
    flow.addSink(glycolaldehyde)

    # define the autocatalysis constraint
    flow.overallAutocatalysis.enable()

    # Set objective function
    flow.objectiveFunction = isEdgeUsed




## input requirement


## enumeration?



##### run flow ####
# Find a solution:
flow.findSolutions()

# # you can find more than one solution:
# flow.findSolutions(maxNumSolutions = 5) 

# Show solution information in the terminal
flow.solutions.list()

# Print solutions
flow.solutions.print()


# flowPrinter = hyperflow.Printer()
# flowPrinter.printUnfiltered = False

# save the flow as a file
flow.dump(f"{dg_name}_{flow_name}.flow")
