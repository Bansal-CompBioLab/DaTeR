import numpy
from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
from scipy.optimize import minimize
from scipy.sparse import csr_matrix
from dendropy import Tree
import networkx

from docplex.mp.model import Model

import os

HGT_constraints_file_name = None

# get number of nodes by traversal in the tree
def get_tree_size_by_traversal(tree):
    i = 0
    # traverse tree by post order
    for node in tree.postorder_node_iter():
        i += 1
    return i





# get maximum upper bound =root to tip distance
def get_max_possible_upper_bound(tree):
    # get root to tip distance and return it
    return 2 * tree.seed_node.distance_from_tip()


# create a manual array for linear constraints from two post order numbers of high one and low node
def get_manual_array(number_of_nodes, higher_node_number, lower_node_number):
    manual_array = numpy.full((number_of_nodes + 1), 0.0)
    manual_array[higher_node_number] = 1
    manual_array[lower_node_number] = -1
    return manual_array


# get a hashmap of nodes by post order using  dictionary
def get_dictionary_of_nodes_by_post_order(tree):
    # initiate empty dictionary
    post_order_dictionary = {}
    # create hashmap like dictionary,
    # where key = post order number of node
    # and value = node
    for node in tree.postorder_node_iter():
        post_order_dictionary[node.post_order_number] = node
    return post_order_dictionary


# check if there is a cycle in the manual constraints
def check_cycles(tree, manual_constraints):
    # create an empty set of edges
    edges = []
    # get a hashmap of nodes by post order using  dictionary
    post_order_dictionary = get_dictionary_of_nodes_by_post_order(tree)
    # iterate all the nodes
    for node in tree.postorder_node_iter():
        # make an empty set of descendants list
        node.descendants = []
        # add each child into descendants list
        for child in node.child_node_iter():
            node.descendants.append(child)
    # get number of manual constraints
    number_of_constraints = len(manual_constraints)
    for i in range(number_of_constraints):
        # create a descendant relationship for each pair of manual constraints
        # get higher and lower node numbers first
        higher_node_number = manual_constraints[i][0]
        lower_node_number = manual_constraints[i][1]
        # get the nodes from the dictionary
        higher_node = post_order_dictionary[higher_node_number]
        lower_node = post_order_dictionary[lower_node_number]
        # lower node is a descendant of the higher node
        higher_node.descendants.append(lower_node)
        # lower_node.descendants.append(higher_node)
    # iterate all the nodes
    for node in tree.postorder_node_iter():
        # create edge for each node to each of its descendant

        for descendant in node.descendants:
            new_edge = (str(node.post_order_number), str(descendant.post_order_number))
            edges.append(new_edge)
            # print(edges)
    # print(edges)
    # create a directed graph using networkx library using our edge list
    G = networkx.DiGraph(edges)

    # checker of cycle
    cycle_flag = False
    ''' 
    # if there is a cycle then print it
    for cycle in networkx.simple_cycles(G):
        cycle_flag = True
        print(cycle)
    if cycle_flag:
        print("Cycle exists!!")
        return True
    '''
    cycles = []
    try:
        cycles = list(networkx.find_cycle(G, orientation="original"))
    except:
        pass

    if len(cycles) > 0:
        cycle_flag = True
        for cycle in cycles:
            print("Cycle exists!")
            print(cycle)
        return True
    # no cycle, return False
    else:
        return False


# get a node object of tree using its label
def get_node_from_label(tree, label):
    # print(tree)
    # find node with its label
    node = tree.find_node_with_label(label)
    if node is None:
        # if not found then find node with taxon label
        node = tree.find_node_with_taxon_label(label)
    return node


# Add manual constraints as post order number pairs, first ones time is greater than the second one, 20 constraints are here
# read it from a file
def get_contstraints_from_file(tree, constraints_file_name):
    if constraints_file_name is None:
        return []
    # open the file named constraints_file_name to get true values for each node of the tree
    file = open(constraints_file_name, "r")
    # initiate empty set of constraints
    constraint_pairs = []
    # for each line, get the post order number of that two nodes
    i = 0
    for line in file:
        i += 1
        # split the line into parts
        words = line.split()
        # first part is the first node label
        first_node_label = words[0]
        # second part is the second node label
        second_node_label = words[1]
        # get first node from the first node label
        first_node = get_node_from_label(tree, first_node_label)
        # print(first_node_label)
        # get second node from the second node label
        second_node = get_node_from_label(tree, second_node_label)
        # obtain first_node and second node post order numbers
        first_node_post_order_number = first_node.post_order_number
        second_node_post_order_number = second_node.post_order_number
        # make a constraint pair using two node's post order numbers
        manual_constraint = [first_node_post_order_number, second_node_post_order_number]
        # append the constraint to the constraints list
        constraint_pairs.append(manual_constraint)
    # print(constraint_pairs)
    return constraint_pairs


# create all the linear constraints for optimization
def get_linear_constraints(tree, constraints_file_name):
    # get maximum upper bound = maximum input sample time of trees
    max_upper_bound = get_max_possible_upper_bound(tree)
    # get number of nodes
    number_of_nodes = get_tree_size_by_traversal(tree)

    parent_child_upper_bounds = []
    parent_child_lower_bounds = []
    all_parent_child_constraints = []
    # iterate over all the nodes and set parent child constraints
    for node in tree.levelorder_node_iter():
        if node is not tree.seed_node:
            # create an empty array of n+1 nodes
            index_array = numpy.full((number_of_nodes + 1), 0.0)
            node_index = node.post_order_number
            parent_index = node.parent_node.post_order_number
            # print(parent_index)
            # child time > parent time in linear constraints
            index_array[node_index] = -1
            index_array[parent_index] = 1
            # set minimum difference of parent and child near to zero
            parent_child_lower_bounds.append(-0.0001)
            # Maximum parent child difference can be maximum upper bound = maximum input sample time of trees
            parent_child_upper_bounds.append(max_upper_bound)
            all_parent_child_constraints.append(index_array)
            # break
    # Add manual constraints as post order number pairs, first ones time is greater than the second one, 20 constraints are here
    manual_constraints = get_contstraints_from_file(tree, constraints_file_name)

    # check cycles in manual constraints
    exists_cycles = check_cycles(tree, manual_constraints)
    # Get number of constraints
    if exists_cycles == False:
        number_of_constraints = len(manual_constraints)
        for i in range(number_of_constraints):
            # if i==11:
            #    continue
            # set minimum difference of parent and child near to zero
            parent_child_lower_bounds.append(-0.0001)
            # Maximum parent child difference can be maximum upper bound = maximum input sample time of trees
            parent_child_upper_bounds.append(max_upper_bound)
            # Create manual constraints using the manual costraints array from above
            higher_node_number = manual_constraints[i][0]
            lower_node_number = manual_constraints[i][1]
            manual_array = get_manual_array(number_of_nodes, higher_node_number, lower_node_number)
            # append the manual constraints
            all_parent_child_constraints.append(manual_array)
    else:
        print("Cycles Exist!!")
    # return all types of linear constraints
    return LinearConstraint(csr_matrix(numpy.array(all_parent_child_constraints)), parent_child_lower_bounds,
                            parent_child_upper_bounds, keep_feasible=False)




def rosen(x):
    sum = 0
    # get HGT constraints
    HGT_constraints_list = get_contstraints_from_file(global_tree, HGT_constraints_file_name)
    for node in global_tree.levelorder_node_iter():
        if node is not global_tree.seed_node:
            node_index = node.post_order_number
            parent_index = node.parent_node.post_order_number
            zero_val_constant = .000000001
            branch_length = node.edge_length
            # print(branch_length)
            if branch_length == 0:
                branch_length = branch_length + zero_val_constant

            penalty = 0
            if x[parent_index] - x[node_index] > 0:
                penalty = numpy.log((x[parent_index] - x[node_index]) / branch_length)
            # elif x[node_index]-x[parent_index] == 0:
            #    penalty = numpy.log(zero_val_constant)
            else:
                penalty = numpy.log(zero_val_constant)  # *(x[node_index]-x[parent_index] )

            # SLRB
            sum += (numpy.sqrt(branch_length)) * (penalty) ** 2

    return sum


def jacobian(x):
    """The Jacobian function"""
    sum = 0
    number_of_nodes = get_tree_size_by_traversal(global_tree)
    jacobian = numpy.full((number_of_nodes + 1), 0.0)

    for node in global_tree.levelorder_node_iter():
        if node is not global_tree.seed_node:
            node_index = node.post_order_number
            parent_index = node.parent_node.post_order_number
            branch_length = node.edge_length
            jacobian[node_index] = -2 * branch_length * numpy.sqrt(branch_length) * numpy.log(
                (x[parent_index] - x[node_index]) / branch_length) * (
                                               branch_length / (-x[node_index] + x[parent_index]))
    return jacobian




# set default bounds for each taxon where a node with a given sample time and root node of the tree
def get_default_bounds(tree):
    # get number of nodes and create two empty arrays of lower and upper bounds
    # max upper bound = max sampling time
    max_upper_bound = get_max_possible_upper_bound(tree)
    number_of_nodes = get_tree_size_by_traversal(tree)
    lower_bounds = numpy.full((number_of_nodes + 1), 0.0)
    upper_bounds = numpy.full((number_of_nodes + 1), max_upper_bound)
    # iterate over each node in the tree
    for node in tree.postorder_node_iter():
        node_index = node.post_order_number

        # if it has a given sampling time then fix the node with
        # that fixed sampling time equal to its upper and lower bound
        if node.sampling_time is not None:
            upper_bounds[node_index] = node.sampling_time
            lower_bounds[node_index] = node.sampling_time


        # for leaf, set the lower and upper bound to zero
        elif node.is_leaf():
            lower_bounds[node_index] = 0
            upper_bounds[node_index] = 0



    lower_bounds[0] = 1.0
    upper_bounds[0] = 1.0

    # combine all the bounds
    bounds = Bounds(lower_bounds, upper_bounds)

    return bounds


global_tree = None


# get the name of node from post order number
def get_name_from_index(node_index, tree):
    for node in tree.postorder_node_iter():
        if node.post_order_number == node_index:
            return node.label



# check how many constraints of parent child relation are broken by our output
def check_constraints_of_results(X, tree, HGT_constraint_file_name):
    # traverse each node and get calculated time from the output
    for node in tree.postorder_node_iter():
        node.calculated_time = X[node.post_order_number]
    parent_child_constraint_count = 0
    # count in how many nodes the parent node's calculated time is greater than the child node's calculated time
    for node in tree.postorder_node_iter():

        if node is not tree.seed_node and node.calculated_time - node.parent_node.calculated_time > 0.0001:
            print("Parent child Constraint Broken")
            print(node)
            parent_child_constraint_count += 1
    # print the number of constraints broken
    # print("Total parent child constraints broken: "+str(parent_child_constraint_count))

    # conut how many manual costraints got broken
    HGT_constraint_count = 0
    manual_constraints = get_contstraints_from_file(tree, HGT_constraint_file_name)
    # get HGT constraint pairs
    number_of_manual_constraints = len(manual_constraints)
    for i in range(number_of_manual_constraints):
        first_node_index = manual_constraints[i][0]
        second_node_index = manual_constraints[i][1]
        first_node_time = X[first_node_index]
        second_node_time = X[second_node_index]
        # if an HGT constraint broken then count it
        if first_node_time + 0.0001 <= second_node_time:
            HGT_constraint_count += 1
            first_node_name = get_name_from_index(first_node_index, tree)
            second_node_name = get_name_from_index(second_node_index, tree)
            # print("HGT Constraint broken between: "+str(first_node_name)+ " time: "+str(first_node_time)+" and "+str(second_node_name)+" time: "+str(second_node_time))

    # print("Total HGT constraints broken: " + str(HGT_constraint_count))
    calibration_point_count = 0



    total_constraint_broken = parent_child_constraint_count + HGT_constraint_count + calibration_point_count
    # print("Total number of calibration point constraint broken: "+ str(calibration_point_count))
    return HGT_constraint_count, calibration_point_count, total_constraint_broken


# set  bounds for each taxon where a node with a given sample time and root node of the tree
def get_lower_upper_bounds(tree):
    number_of_nodes = get_tree_size_by_traversal(tree)
    lower_bounds = numpy.full((number_of_nodes) + 1, 0.0)

    upper_bounds = numpy.full((number_of_nodes) + 1, 1.0)

    for node in tree.leaf_node_iter():
        node_index = node.post_order_number
        lower_bounds[node_index] = 1.0
        if node.sampling_time is not None:
            upper_bounds[node_index] = node.sampling_time
            lower_bounds[node_index] = node.sampling_time
        # upper_bounds[node_index] = 1.0

    root_index = tree.seed_node.post_order_number
    lower_bounds[root_index] = 0.0
    upper_bounds[root_index] = 0.0
    return lower_bounds, upper_bounds






def initialize_from_original_tree(tree):
    # get number of nodes
    number_of_nodes = get_tree_size_by_traversal(tree)
    # create an empty array with n + 1 nodes
    X = numpy.full((number_of_nodes + 1), 0.0)

    X[0] = 0.9
    for node in tree.postorder_node_iter():
        time = node.distance_from_tip()
        node_index = node.post_order_number
        X[node_index] = time
    return X


def debug_branch_length(x):
    sum = 0
    for node in global_tree.levelorder_node_iter():
        if node is not global_tree.seed_node:
            node_index = node.post_order_number
            parent_index = node.parent_node.post_order_number
            zero_val_constant = .000001
            branch_length = node.edge_length
            # print(branch_length)
            if branch_length == 0:
                branch_length = branch_length + zero_val_constant

            penalty = 0
            if x[node_index] > x[parent_index]:
                penalty = numpy.log(numpy.abs(x[0]) * (x[node_index] - x[parent_index]) / branch_length)
            elif x[parent_index] > x[node_index]:
                penalty = numpy.log(numpy.abs(x[0]) * (x[parent_index] - x[node_index]) / branch_length)
            else:
                penalty = numpy.log(zero_val_constant)

            # penalty = (numpy.abs(x[node_index]-x[parent_index])*100.0/branch_length)*numpy.abs(x[node_index]-x[parent_index])
            sum += (numpy.sqrt(branch_length + 30)) * (penalty) ** 2
            print(
                str(numpy.sqrt(branch_length + 30)) + " " + "  " + str(numpy.sqrt(branch_length)) + " penalty: " + str(
                    penalty ** 2))
            # sum += 1/(branch_length)* (branch_length-(x[0]*(x[node_index] - x[parent_index]))) ** 2
            # sum += numpy.sqrt(branch_length+.0001 ) * (numpy.log(x[0]*(x[node_index] - x[parent_index])/branch_length))**2
    return sum


# get number of HGT constraints
def get_number_of_constraints(constraints_file_name):
    file = open(constraints_file_name, "r")
    count = 0
    # initiate empty set of constraints
    constraint_pairs = []
    # for each line, get the post order number of that two nodes
    for line in file:
        count += 1
    # print(constraint_pairs)
    return count





# get all nodes except the root
def get_non_root_nodes(tree):
    node_list = []
    for node in tree.postorder_node_iter():
        if node is not tree.seed_node:
            node_list.append(node)
    return node_list


# get all nodes except the root
def get_non_tip_nodes(tree):
    node_list = []
    for node in tree.postorder_node_iter():
        if node.is_leaf() == False:
            node_list.append(node)
    return node_list





# Add all linear constraints in docplex
def add_linear_constraints_in_docplex(model, docplex_variables, tree, constraints_file_name):
    non_root_nodes = get_non_root_nodes(tree)
    for node in non_root_nodes:
        node_index = node.post_order_number
        node_parent_index = node.parent_node.post_order_number
        # the first constraint
        c1 = docplex_variables[node_parent_index - 1] - docplex_variables[node_index - 1] >= 0
        c1.set_mandatory()
        model.add_constraint(c1)

        if node.is_leaf():
            c2 = (docplex_variables[node_index - 1] == 0)
            c2.set_mandatory()
            model.add_constraint(c2)

    manual_constraints = get_contstraints_from_file(tree, constraints_file_name)
    i = 0
    for manual_constraint in manual_constraints:
        first_node_index = manual_constraint[0]
        second_node_index = manual_constraint[1]
        c3 = docplex_variables[first_node_index - 1] - docplex_variables[second_node_index - 1] >= 0
        c3.set_mandatory()
        model.add_constraint(c3)

        i += 1


# solving the problem using docplex
def add_up_soft_constraints(model, docplex_variables, tree, constraints_file_name):
    HGT_broken_count = model.integer_var(name='Soft_penalty')
    manual_constraints = get_contstraints_from_file(tree, constraints_file_name)
    i = 0
    for manual_constraint in manual_constraints:
        first_node_index = manual_constraint[0]
        second_node_index = manual_constraint[1]
        # model.add_constraint(docplex_variables[first_node_index - 1] - docplex_variables[second_node_index - 1] >= 0)

        i += 1
    return HGT_broken_count


def solve_problem_by_docplex( tree, constraints_file_name, optimization_model):
    print("Calling CPLEX library..")
    # create model
    model = Model(name='dating')

    # create an empty set of variables
    docplex_variables = []
    # iterate over all the variables in old model
    for node in tree.postorder_node_iter():
        node_name = node.label
        if node_name is None:
            node_name = node.taxon.label

        docplex_variables.append(model.continuous_var(name=node_name))
    ''' 
    #add some dummy variables for test here
    for i in range(0,1000):
        docplex_variables.append(model.continuous_var(name=str(i)))
    '''
    # get all nodes except the root
    non_root_node_list = get_non_root_nodes(tree)
    non_tip_node_list = get_non_tip_nodes(tree)

    add_linear_constraints_in_docplex(model, docplex_variables, tree, constraints_file_name)


    objective = None

    if optimization_model == 'SDD':
        objective = model.sum(1.0 / (node.distance_from_tip()) *
                              (node.distance_from_tip() - docplex_variables[node.post_order_number - 1]) *
                              (node.distance_from_tip() - docplex_variables[node.post_order_number - 1]) for node in
                              non_tip_node_list)
    elif optimization_model == 'SBD':
        # SBD
        objective = model.sum((1.0 / node.edge_length) * (
                    node.edge_length - docplex_variables[node.parent_node.post_order_number - 1] + docplex_variables[
                node.post_order_number - 1]) *
                              (node.edge_length - docplex_variables[node.parent_node.post_order_number - 1] +
                               docplex_variables[node.post_order_number - 1]) for node in non_root_node_list)

    model.context.cplex_parameters.threads = 2
    # HGT_penalty_sum = add_up_soft_constraints(model, docplex_variables,  tree, constraints_file_name)
    model.minimize(objective)
    solution = model.solve()
    #print(solution.is_feasible_solution())

    values = solution.get_values(docplex_variables)

    values.insert(0, 0)
    #print(values)
    return values




def get_edge_list_in_post_order(tree):
    edge_list = [None] * get_tree_size_by_traversal(tree)
    for node in tree.postorder_node_iter():
        node_index = node.post_order_number - 1
        edge_length = node.edge_length
        edge_list[node_index] = edge_length
    return edge_list




def solve(tree, constraints_file_name, model):
    global global_tree
    global_tree = tree
    global HGT_constraints_file_name
    HGT_constraints_file_name = constraints_file_name

    # Get input from original tree
    X2 = initialize_from_original_tree(tree)

    # get all the linear constraints
    linear_constraint = get_linear_constraints(tree, constraints_file_name)
    # Get the max and min bounds for constrained optimization
    bounds = get_default_bounds(tree)

    solution = None

    # try to solve by docplex
    if model == 'SDD' or model=='SBD':
        solution = solve_problem_by_docplex( tree, constraints_file_name, model)
    else:
        print("Calling trust-costr method of Scipy library..")
        res = minimize(rosen, X2, method='trust-constr',
                       constraints=[linear_constraint],
                       bounds=bounds,
                       options={'verbose': 1,
                                'maxiter': 10000,
                                'disp': True,

                                }
                       )
        solution = res.x




    return solution


