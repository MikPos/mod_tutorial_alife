FORBIDDEN_SUBGRAPHS = []
MAX_NUMBER_OF_VERTICES = 10 
MAX_NUMBER_OF_RADICALS = 1
# This restricts the derivation to graphs not containing any of the forbidden subgraphs
cyc3 = Graph.fromDFS("[*]1{*}[*]{*}[*]{*}1")
FORBIDDEN_SUBGRAPHS.append(cyc3)
cyc4 = Graph.fromDFS("[*]1{*}[*]{*}[*]{*}[*]{*}1")
FORBIDDEN_SUBGRAPHS.append(cyc4)
double_bond = Graph.fromDFS("[*]1=[*]2=[*]3")
FORBIDDEN_SUBGRAPHS.append(double_bond)

def forbiddenSubgraphs(derivation):
    return all(subgraph.monomorphism(g, labelSettings=ls) == 0
                      for subgraph in FORBIDDEN_SUBGRAPHS for g in derivation.right)


def has_more_radicals_than(derivation):
    radical_counts = []
    index = 0
    for g in derivation.right:
        radical_counts.append(0)
        for v in g.vertices:
            radical_counts[index] += "." in v.stringLabel
        index += 1
    return not any(count > MAX_NUMBER_OF_RADICALS for count in radical_counts)

def has_at_most_x_vertices(derivation):
    return all(g.numVertices <= MAX_NUMBER_OF_VERTICES for g in derivation.right)


def all_constraints_apply(constraint_functions, derivation):
    return all(function(derivation) for function in constraint_functions)


CONSTRAINT_FUNCTIONS = [forbiddenSubgraphs, has_more_radicals_than, has_at_most_x_vertices]
