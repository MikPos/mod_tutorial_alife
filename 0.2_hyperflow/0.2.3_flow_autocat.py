### flow query on a loaded dg 
include("formose.py")

dg_name = "dg_autocat"
dg = DG.load(inputGraphs, inputRules, f"{dg_name}.dg") 

# explore the size of the dg
def display_size_dg(dg):
    print("============ Number vertices ============")
    print(dg.numVertices)
    print("============ No edges ============")
    print(dg.numEdges)

display_size_dg(dg)

### set up a flow model 
flow_name = "flow_autocat"

flow = hyperflow.Model(dg) 

flow.addSource(formaldehyde)
flow.addSource(glycolaldehyde)

flow.addSink(glycolaldehyde) 

flow.addConstraint(inFlow[formaldehyde] == 2)
flow.addConstraint(inFlow[glycolaldehyde] == 1)
flow.addConstraint(outFlow[glycolaldehyde] == 2)

flow.allowIOReversal = False


#### autocatalysis
# Enables the autocatalysis search for autocatalytic cylces
flow.overallAutocatalysis.enable()
####

flow.objectiveFunction = isEdgeUsed

##### run flow ####
flow.findSolutions()
flow.solutions.list()
flow.solutions.print()
flow.dump(f"{dg_name}_{flow_name}.flow")



