### Loading a dg 

# include the grammar used to create the dg (includes the molecules and the rules used to create the dg)
include("grammar_formose.py")

# Define the filename of the DG to load.
dg_name = "dg_ex_3"

# load the dg in question from the previous dg dump file
# inputGraphs and inputRules are parsed from the grammar file that was loaded earlier
dg = DG.load(inputGraphs, inputRules, f"{dg_name}.dg") 

# explore the size of the dg
def display_size_dg(dg):
    # print the number of vertices of the dg
    print("============ Number vertices ============")
    print(dg.numVertices)
    # print the number of edges of the dg
    print("============ No edges ============")
    print(dg.numEdges)

display_size_dg(dg)
