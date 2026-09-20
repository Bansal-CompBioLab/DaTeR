from dendropy import Tree
import  preprocessor
import numpy


# get LCA of two nodes manually
def get_LCA( root, first_taxon, second_taxon):
    if root is None:
        return None
    if root == first_taxon or root==second_taxon:
        return root

    elif root.is_leaf():
        return None
    left_child = root.child_nodes()[0]
    right_child = root.child_nodes()[1]

    left_LCA = get_LCA(left_child, first_taxon, second_taxon)
    right_LCA = get_LCA(right_child, first_taxon, second_taxon)

    if left_LCA is not None and right_LCA is not None:
        return root
    else:
        if left_LCA is not None:
            return left_LCA
        else:
            return right_LCA


# set sampling times from the input file sampling.txt
def set_sampling_time_from_real_data(sampling_time, tree):
    # set every node's sampling time in tree to None first
    for node in tree.postorder_node_iter():
        node.sampling_time = None
        node.sampling_time_upper = None
        node.sampling_time_lower = None

    # open input sampling time file by name
    file = open(sampling_time, "r")


    # for each line in sampling time file
    for line in file:
        # split the line by whitespaces
        words = line.split()

        if len(words) == 1:
            # get total number of constraints
            total_sampling_times  = int(words[0])+ int((len(tree.nodes())+1)/2)

            continue
        # get the taxon name or internal node name
        first_taxon_name = words[0]
        first_taxon = tree.find_node_with_taxon_label(first_taxon_name)

        # get the taxon name or internal node name
        second_taxon_name = words[1]
        second_taxon = tree.find_node_with_taxon_label(second_taxon_name)

        #print(first_taxon)
        # print(second_taxon)
        #taxon_labels = [first_taxon_name, second_taxon_name]
        #node = tree.mrca(taxon_labels=taxon_labels)
        #print(node)
        # get LCA manually from root to tip path from stackOverFlow.com
        node = get_LCA(tree.seed_node, first_taxon, second_taxon)
        #print(node)

        # if the node has both upper and lower bound then add them as an attribute to the node
        if len(words) == 4:
            #print(words)
            # save lower and upper bound as attributes in the node (converting difference with maximum upper bound value)
            upper_bound = float(words[2])
            lower_bound = float(words[3])
            # if no upper bound then it will be coverted to maximum upper bound eventually
            if upper_bound == -1:
                upper_bound = 2*tree.seed_node.distance_from_tip()
            # if no lower bound then it will be coverted to zero eventually
            if lower_bound == -1:
                lower_bound = 0


            node.sampling_time_lower = lower_bound
            node.sampling_time_upper = upper_bound



#get tree from filename
def get_tree_from_filename(file_name):
    #get input file from input file name
    input_file = open(file_name,'r' )
    # get  tree string
    tree_string = input_file.readline()
    # get the newick tree from the dendropy library
    tree = Tree.get(data=tree_string, schema='newick', preserve_underscores=True)

    return tree
def get_combined_tree_from_filename(filename):
    tree = get_tree_from_filename(filename)
    preprocessor.set_traversal_numbers(tree)
    labeled_tree = get_tree_from_filename('..\\labels.txt')
    tree = preprocessor.get_combined_tree(labeled_tree, tree)
    return tree
# get equivalent node from another tree by name
def get_equivalent_node_by_name(node, another_tree):
    # if the node is internal
    if node.is_leaf()==False:
        # fina node by internal node label
        equivalent_node = another_tree.find_node_with_label(node.label)

    else:
        # as the node is a taxon, so find equivalent node with taxon label
        equivalent_node = another_tree.find_node_with_taxon_label(node.taxon.label)
    return equivalent_node


# print detailed and compared results with true output from the final output
def print_results(results, tree, output_filename):
    for node in tree.postorder_node_iter():
        #convert distance from subtrating maximum upper bound value
        node.output_date =  results[node.post_order_number]
    for node in tree.postorder_node_iter():
        #get true node from output true tree
        #true_node = get_equivalent_node_by_name(node, output_true_tree)
        if node is not tree.seed_node:
            #save original date and edge length to write later
            node.original_date = node.distance_from_tip()
            node.original_edge_length = node.edge_length
        else:
            # node is root
            node.original_date = node.distance_from_tip()
    for node in tree.postorder_node_iter():
        if node.parent_node!=None:
            node.edge_length = node.parent_node.output_date - node.output_date
    output_treestring = tree.as_string(schema="newick")
    output_file = open(output_filename, "a")
    output_file.write(output_treestring)


