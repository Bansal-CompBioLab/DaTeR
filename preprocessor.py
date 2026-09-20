
from dendropy import Tree
#set pre, post and inorder node traversal numbers in the tree
def set_traversal_numbers(tree):
    i = 1
    for node in tree.postorder_node_iter():
        node.post_order_number = i
        i+=1
    #print("Number of nodes by traversal: "+ str(i-1))
    i = 1
    for node in tree.preorder_node_iter():
        node.pre_order_number = i
        i += 1

    i = 1
    for node in tree.inorder_node_iter():
        node.in_order_number = i
        i += 1
        # set sampling time as none
        node.sampling_time = None

#get fixed tree by combining chronogram file and labels file
def compare_labeled_node_with_chronogram_node( chronogram_node, labeled_node, chronogram_tree, labeled_tree):
    #make an empty list of leaves under labeled node
    leaf_list_of_labeled_node = []
    #get all the leaves under labeled node
    for leaf_of_labeled_node in labeled_node.leaf_iter():
        leaf_list_of_labeled_node.append(leaf_of_labeled_node.taxon.label)

    # make an empty list of leaves under labeled node
    leaf_list_of_chronogram_node = []
    # get all the leaves under chronoogram node
    for leaf_of_chronogram_node in chronogram_node.leaf_iter():
            leaf_list_of_chronogram_node.append(leaf_of_chronogram_node.taxon.label)

    if set(leaf_list_of_chronogram_node) != set(leaf_list_of_labeled_node):
        return False
    else:
        return True


def get_combined_tree(labeled_tree, input_tree):

    #for each chronogram node
    for input_node in input_tree.postorder_node_iter():
        # set sampling time of input node as None
        input_node.sampling_time = None
        input_node.sampling_time_lower = None
        input_node.sampling_time_upper = None
        #skip leaf nodes
        if input_node.is_leaf():
            continue
        # for each internal node in the labeled tree
        for labeled_node in labeled_tree.postorder_node_iter():
            # get equivalent chronogram node from the labeled node
            if  labeled_node.is_leaf() == False:
                if compare_labeled_node_with_chronogram_node(input_node, labeled_node, input_tree, labeled_tree):
                    input_node.label = labeled_node.label
                    break


    # return chronogram tree as we have now the labels imported from the labeled  tree
    return input_tree

