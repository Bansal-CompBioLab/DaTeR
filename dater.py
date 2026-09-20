
from dendropy import Tree
import argparse
import preprocessor
import solver
import utility_operations as utility
import statistics



parser = argparse.ArgumentParser()

parser.add_argument("-i" ,"--input" ,required=True ,help="Input tree")
parser.add_argument("-o" ,"--output" ,required=True ,help="The output file name")
parser.add_argument("-m" ,"--model" ,required=False
                    ,help="The model for constrained optimization: select one from SBD/SDD/SLRB, default model is SLRB")
parser.add_argument("-c" ,"--constraints" ,required=True ,help="List of relative constraints")

args = vars(parser.parse_args())


output_filename = args["output"]
model = args['model'] if args['model'] else 'SLRB'

constraints_file_name = args["constraints"] if args["constraints"] else None

def get_tree_size_by_traversal(tree):
    i = 0
    for node in tree.postorder_node_iter():

        i+=1
    return i


def get_aggregated_output_tree(output_filename):
    output_file = open(output_filename, 'r')
    tree_strings = output_file.readlines()
    dates = []
    all_dates = []

    for tree_string in tree_strings:
        tree = Tree.get(data=tree_string, schema='newick', preserve_underscores=True)

        preprocessor.set_traversal_numbers(tree)

        while len(all_dates) < get_tree_size_by_traversal(tree):
            all_dates.append([])

        for node in tree.postorder_node_iter():
            index = node.post_order_number
            date = node.distance_from_tip()

            all_dates[index - 1].append(date)

    for date_list in all_dates:

        dates.append(statistics.mean(date_list))

    for node in tree.postorder_node_iter():
        if node is not tree.seed_node:
            index = node.post_order_number
            parent_index = node.parent_node.post_order_number
            node.edge_length = dates[parent_index - 1] - dates[index - 1]
        else:
            node.edge_length = 0


    return tree



with open(args["input"] ,'r') as fin:
    tree_strings = fin.readlines()


i=0

for treestr in tree_strings:
    i += 1

    input_tree = Tree.get(data=treestr, schema='newick', preserve_underscores=True)

    preprocessor.set_traversal_numbers(input_tree)


    results= solver.solve(input_tree, constraints_file_name, model)


    utility.print_results(results, input_tree, output_filename)
    print("Writing output tree #"+str(i)+" to "+output_filename)

#if we have more than one input tree then output an aggregated output tree as well
if i>1:
    aggregated_tree = get_aggregated_output_tree(output_filename)
    print("Writing aggregated dated tree to "+"aggregated_"+output_filename)
    aggregated_tree.write(path="aggregated_"+output_filename,schema="newick")


