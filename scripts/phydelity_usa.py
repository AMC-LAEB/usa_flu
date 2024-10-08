#!/usr/bin/env python

# Phydelity
# Author: Alvin X. Han & Edyth Parker

from __future__ import division
from scipy.stats import mannwhitneyu
import re
import argparse
import numpy as np
import pandas as pd
import itertools
import time

import pyximport; pyximport.install()
from phyilpd import node_leaves_reassociation, clean_up_modules, phydelity_output
from phyilpd.tree_utils import parse_newick_tree
from phyilpd.stats_utils import qn, get_cluster_size_distribution

if __name__ == '__main__':
    # parse parameters
    version = 2.0
    parser = argparse.ArgumentParser(description='Phydelity v{}'.format(version))

    required_args = parser.add_argument_group('Required')
    required_args.add_argument('-t', '--tree', type=str, help='Input phylogenetic tree in NEWICK format.')

    analyses_options = parser.add_argument_group('Analysis options')
    analyses_options.add_argument('--wcl', type=np.float64, help='Hard WCL')
    analyses_options.add_argument('--k', type=int, help='Custom k neighbours (optional).')
    analyses_options.add_argument('--outgroup', type=str, default=False, help='Taxon (name as appeared in tree) to be set as outgroup OR type \'midpoint\' for midpoint-rooting.')
    analyses_options.add_argument('--collapse_zero_branch_length', action='store_true', help='Collapse internal nodes with zero branch length of tree before running Phydelity.')
    analyses_options.add_argument('--equivalent_zero_length', default=1e-6, type=np.float64, help='Maximum branch length to be rounded to zero if the --collapse_zero_branch_length flag is passed (default = %(default)s).')
    analyses_options.add_argument('--date_last_sample', default=False, type=np.float64, help='Last date of sampling of included taxa in decimal format')
    analyses_options.add_argument('--min_date_cluster', default=False, type=np.float64, help='Minimum date for a node to be considered a potential cluster')
    analyses_options.add_argument('--coal_time', default=False, type=np.float64, help='Time period from cluster root that must see at least the percentage of within-cluster coalescent events specified by --coal_percentage')
    analyses_options.add_argument('--coal_percentage', default=False, type=np.float64, help='Minimum percentage of coalescent events to happen within specified time period')

    solver_options = parser.add_argument_group('Solver options')
    """
    solver_options.add_argument('--solver', default='gurobi', choices=['glpk', 'gurobi'], type=str, help='Preferred ILP solver IF more than one solvers are available (default: %(default)s).')
    solver_options.add_argument('--solver_verbose', default=0, choices=[0, 1], type=int, help='ILP solver verbose (default: %(default)s)')
    solver_options.add_argument('--solver_check', action='store_true', help='Check available ILP solver(s) installed.')
    """
    solver_options.add_argument('--solver_verbose', default=0, choices=[0, 1], type=int, help='Gurobi solver verbose (default: %(default)s)')
    solver_options.add_argument('--solver_check', action='store_true', help='Check if Gurobi is installed.')

    output_options = parser.add_argument_group('Output options')
    output_options.add_argument('--pdf_tree', action='store_true', help='PDF tree output annotated with cluster results (X server required).')

    # solver_options.add_argument('--threads', type=int, help='Number of threads to use (default = all).')

    params = parser.parse_args()
    # need to change this if we ever update with glpk support
    params.solver = 'gurobi'

    print ('{}\n\n{:^72}\n{:^72}\n\n{}'.format(''.join(['-']*72), 'Phydelity', 'v{}'.format(version), ''.join(['-']*72)))

    # check solver availability
    available_solvers = {}
    """
    try:
        # check for glpsol
        cmd = ['glpsol', '--version']
        solver_version = 'glpk_{}'.format(re.search('v\d+\.\d+', subprocess.check_output(cmd)).group())
        available_solvers['glpk'] = solver_version
    except:
        pass
    """

    try:
        from gurobipy import gurobi
        solver_version = 'gurobi_v{}'.format('.'.join(map(str, gurobi.version())))
        available_solvers['gurobi'] = solver_version
    except:
        pass

    if len(available_solvers) > 0:
        # exit if only to check available ILP solver(s) installed
        if params.solver_check:
            print ('\nAvailable solvers...{}\n'.format(', '.join(available_solvers.values())))
            exit(0)
    else:
        raise Exception('\nNo supported solvers installed. See manual for details on how to download and install supported ILP solvers.\n')

    # check if tree input is given minimally
    if not params.tree:
        raise Exception('\nTree input missing. Check --help for options.\n')

    # preferred solver
    if params.solver not in available_solvers:
        print ('\nWARNING: {} is not installed.'.format(params.solver))
        params.solver = available_solvers.keys()[0]

    print('\nILP solver...{}'.format(available_solvers[params.solver]))

    """
    # limit number of threads for parallelisation
    try:
        ncpu = int(params.threads)
    except:
        from pathos.helpers import mp
        ncpu = mp.cpu_count()
    print ('Threads...{}'.format(ncpu))
    """

    print ('\nParsing tree...')
    # filenames
    treefname = re.sub('([^/]+/|\.[^.]+$)', '', params.tree)

    # parse newick tree file
    try:
        newick_tree_string = parse_newick_tree(params.tree)
    except:
        raise Exception('\nInvalid tree file.\n')

    # parse for tree distance/toplogical information
    from phyilpx import phyilpx_treeinfo
    phyilpx_obj = phyilpx_treeinfo(newick_tree_string, treefname, params.outgroup, params.collapse_zero_branch_length, params.equivalent_zero_length, params.date_last_sample)

    global_leaf_node_id_to_leafname, original_tree_string, global_node_to_child_leaves, global_node_to_child_nodes = phyilpx_obj.properties()


    # get pairwise distance array (nodepair = all nodes including leaves, leafpair = structured, just leaves)
    global_nodepair_to_dist, global_leafpair_to_distance = phyilpx_obj.get_nodepair_distance()

    # get structured array of leaf dist to node / array of node to leaves (reverse sorted by distance to node)
    global_leaf_dist_to_node, global_node_to_leaves = phyilpx_obj.get_leaf_dist_to_node()

    # structured array of node to parent node
    global_node_to_parent_node = phyilpx_obj.get_node_to_parent_node()

    global_node_to_date = phyilpx_obj.get_node_to_date()

    # ancestry relations/global_node_to_mean_pwdist are python dictionaries
    # global_node_to_mean_child_dist2anc = np.array(N,N)
    # global_node_to_pwdist = dictionary of array
    global_node_to_ancestral_nodes, global_node_to_descendant_nodes, global_leaf_to_ancestors, global_node_to_mean_child_dist2anc, global_node_to_pwdist, global_node_to_mean_pwdist = phyilpx_obj.get_ancestral_relations()

    # !-- phydelity specific --#
    print ('\nCalculating distance distribution of closely-related tips...')

    # hard WCL imposed
    if params.wcl:
        wcl = params.wcl
    else:
        # determine distance distribution of closely-related tips
        if params.k:
            k_range = [params.k]
            print ('WARNING: k fixed at {}.'.format(params.k))
        else:
            k_range = range(2, 6) # 2<= autoscaled k <= 5
            from phyilpx import p_hypotest

        # First, we determine the precision to which we should round pairwise distance to. This is determined based on the
        # median closest pair difference.
        closest_distance_diff_distribution = np.zeros(len(global_leaf_node_id_to_leafname), dtype=np.float64)

        for _, leaf in enumerate(global_leaf_node_id_to_leafname.keys()):
            # get closest neighbouring leaf of leaf
            j_array = global_leafpair_to_distance[global_leafpair_to_distance['leaf_i'] == leaf][['leaf_j', 'dist']]
            closest_leaf = np.sort(j_array, order='dist')['leaf_j'][0]

            # get difference of distance of mrca to leaf and neighbour leaf
            mrca_node = np.max(list(set(global_leaf_to_ancestors[leaf])&set(global_leaf_to_ancestors[closest_leaf])))
            mrca_leaf_dist = global_leaf_dist_to_node[(leaf, mrca_node)]
            mrca_neighbour_dist = global_leaf_dist_to_node[(closest_leaf, mrca_node)]

            closest_distance_diff_distribution[_] = np.abs(mrca_leaf_dist - mrca_neighbour_dist)

        # get median difference of closest pairwise distances
        median_closest_distance_diff = np.median(closest_distance_diff_distribution)

        # if median closest pair distance < 1, then number of demical place to the 2 significant digit will be the deimcal
        # place to round
        if median_closest_distance_diff == 0.:
            decimals_to_round = np.int64(np.abs(np.log10(np.min([_ for _ in closest_distance_diff_distribution if _ > 0.])))) + 2
        elif median_closest_distance_diff < 1:
            decimals_to_round = np.int64(np.abs(np.log10(median_closest_distance_diff))) + 2
        else:
            decimals_to_round = 2

        # construct dictionary of leaf_j(s) for every leaf_i, grouped by their distance to leaf_i
        leaf_to_dist_to_csleaf = {}
        identical_leaves = {}
        for leaf_i, leaf_j in itertools.combinations(global_leaf_node_id_to_leafname.keys(), 2):
            dist = np.round(global_nodepair_to_dist[(leaf_i, leaf_j)], decimals=decimals_to_round)

            if np.absolute(dist - 2*params.equivalent_zero_length) <= 2*params.equivalent_zero_length:
                dist = 0.

            try:
                leaf_to_dist_to_csleaf[leaf_i][dist].append(leaf_j)
            except:
                try:
                    leaf_to_dist_to_csleaf[leaf_i][dist] = [leaf_j]
                except:
                    leaf_to_dist_to_csleaf[leaf_i] = {dist:[leaf_j]}

            try:
                leaf_to_dist_to_csleaf[leaf_j][dist].append(leaf_i)
            except:
                try:
                    leaf_to_dist_to_csleaf[leaf_j][dist] = [leaf_i]
                except:
                    leaf_to_dist_to_csleaf[leaf_j] = {dist:[leaf_i]}

        # iterating over the possible k-values
        for k_strains in k_range:
            leaf_to_kth_sorted_cs_leaves = {}
            for leaf, dist_to_csleaf in leaf_to_dist_to_csleaf.items():
                # sort each leaf to every other leaf by their rounded pairwise distance up to the k-th closest neighbour
                sorted_pw_distances = sorted(dist_to_csleaf.keys())
                if sorted_pw_distances[0] == 0.:
                    sorted_pw_distances = sorted_pw_distances[1:k_strains+1]
                else:
                    sorted_pw_distances = sorted_pw_distances[:k_strains]

                leaf_to_kth_sorted_cs_leaves[leaf] = [(distance, tuple(dist_to_csleaf[distance])) for distance in sorted_pw_distances]

            core_member_pairwise_leafdist = []
            # for every leaf_i, the distance to its k-th closest neighbour leaf_j will be considered a core distance if
            # leaf_i is also <= k-th closest neighbour of leaf_i
            for leaf, sorted_kth_dist_cleaves in leaf_to_kth_sorted_cs_leaves.items():

                dist_to_add = []

                add_to_core_dist_binary = 0
                for distance, cleaves_tuple in sorted_kth_dist_cleaves:
                    dist_to_add.append(distance)

                    found_leaf_binary = 0
                    for cleaf in cleaves_tuple:
                        for (_distance, _cleaves_tuple) in leaf_to_kth_sorted_cs_leaves[cleaf]:
                            if leaf in _cleaves_tuple:
                                found_leaf_binary = 1
                                break

                        if found_leaf_binary == 1:
                            add_to_core_dist_binary = 1
                            break

                    if add_to_core_dist_binary == 1:
                        core_member_pairwise_leafdist += dist_to_add

            # Sort core member distance distribution
            sorted_core_dist = np.sort(list(core_member_pairwise_leafdist))

            # check that the i-th and (i-1)-th pairwise distance do not differ by more than a log
            diff = []
            for _d, d in enumerate(sorted_core_dist):
                if (_d == 0) or (sorted_core_dist[_d-1] <= params.equivalent_zero_length):
                    continue

                if abs(d-sorted_core_dist[_d-1])/sorted_core_dist[_d-1] <= params.equivalent_zero_length:
                    continue

                diff.append(np.log10(abs(d-sorted_core_dist[_d-1])/sorted_core_dist[_d-1]))

                if np.log10(abs(d-sorted_core_dist[_d-1])/sorted_core_dist[_d-1]) > 0.:
                    core_member_pairwise_leafdist = sorted_core_dist[:_d][:]
                    break

            """
            percentile_core = len(core_member_pairwise_leafdist)/len(sorted_core_dist)
            if percentile_core < 0.01:
                core_member_pairwise_leafdist = sorted_core_dist[:]
            """

            med_x = np.median(core_member_pairwise_leafdist)
            mad_x = qn(core_member_pairwise_leafdist)
            # Any distance in the core pairwise distance distribution that is > med_x + mad_x is subsequently removed
            # Hence, core distances of any taxa that are more distantly-related than the ensemble would be excluded.
            core_member_pairwise_leafdist = np.sort([_ for _ in core_member_pairwise_leafdist if _ <= med_x + mad_x])

            if len(k_range) > 1: # auto-scaling of k
                if k_strains > 2:

                    p_val = p_hypotest(core_member_pairwise_leafdist, prev_core_member_pairwise_distance, 1)

                    if p_val < 0.01:
                        wcl = np.amax(prev_core_member_pairwise_distance)
                        k_strains -= 1
                        print ('auto-scaling k...{}'.format(k_strains))
                        break
                    elif k_strains == 5:
                        wcl = np.amax(core_member_pairwise_leafdist)
                        print ('auto-scaling k...{}'.format(k_strains))
                        break

                prev_core_member_pairwise_distance = core_member_pairwise_leafdist[:]
            else:
                wcl = np.amax(core_member_pairwise_leafdist)
    print('MPL...{}'.format(wcl))

    # distal dissociation
    # level-order sorted list of nodes with leaves >= always 2 for transmission clusters (only reassociate such nodes)
    print ('\nDistal dissociation...')
    cs = 2 # min cluster size = pair
    curr_list_of_ancestral_node = np.sort(np.array([node for node, leaves in global_node_to_leaves.items() if len(leaves) >= cs], dtype=np.int64))

    nla_object = node_leaves_reassociation(cs, wcl, curr_list_of_ancestral_node, global_node_to_leaves, global_node_to_ancestral_nodes, global_node_to_descendant_nodes, global_node_to_mean_pwdist, global_node_to_mean_child_dist2anc, global_node_to_parent_node, global_nodepair_to_dist, global_leaf_dist_to_node, global_leaf_to_ancestors, params.equivalent_zero_length, global_leaf_node_id_to_leafname, global_node_to_child_leaves, global_leaf_node_id_to_leafname)

    curr_node_to_leaves, curr_node_to_descendant_nodes, curr_node_to_mean_pwdist = nla_object.nla_main()
        # remove any nodes with len(leaves) < cs after distal dissociation
    for node, leaves in curr_node_to_leaves.items():
        if len(leaves) < cs:
            del curr_node_to_leaves[node]
            continue

    # no nodes remaining
    if len(curr_node_to_leaves) == 0:
        raise Exception("No solution available.")

    # if two nodes have the same set of leaves, retain only the most descendant one
    deleted_nodes = []
    for node_i, node_j in itertools.combinations(curr_node_to_leaves.keys(), 2):
        if node_i in deleted_nodes or node_j in deleted_nodes:
            continue
        if set(curr_node_to_leaves[node_i]) == set(curr_node_to_leaves[node_j]):
            node_to_del = min([node_i, node_j])
            deleted_nodes.append(node_to_del)
            del curr_node_to_leaves[node_to_del]


    for node, leaves in curr_node_to_leaves.items():

        if global_node_to_date[node] < params.min_date_cluster:
            del curr_node_to_leaves[node]    
            continue

        if node in curr_node_to_descendant_nodes:
            allnodes = [[i] * (len(global_node_to_child_nodes[i])-1) for i in curr_node_to_descendant_nodes[node] if i in global_node_to_child_nodes] 
        else:
            allnodes = []
        allnodes.append((len(global_node_to_child_nodes[node]) - 1) * [node]) 
        allnodes = [i for sublist in allnodes for i in sublist]
        node_to_desc_lens = [global_nodepair_to_dist[(node,i)] for i in allnodes]
        thresh = params.coal_percentage
        if len([i for i in node_to_desc_lens if i < params.coal_time]) <  len(node_to_desc_lens) * thresh:
            del curr_node_to_leaves[node]



    # update pairwise distance distributions and ancestry relations
    print ('\nUpdating tree info to ILP model...')
    curr_leaves = [x for y in curr_node_to_leaves.values() for x in y]
    print('here4')
    curr_list_of_ancestral_node = curr_node_to_leaves.keys()[:]
    print('here5')
    curr_node_to_descendant_nodes = {k:list(set(v)&set(curr_list_of_ancestral_node)) for k,v in curr_node_to_descendant_nodes.items() if k in curr_list_of_ancestral_node and len(v) > 0}
    print('here3')
    '''
    for node, leaves in curr_node_to_leaves.items():
        for leaf in leaves:
            leaf_name = global_leaf_node_id_to_leafname[leaf]
            if re.search("NUH0018", leaf_name):
                print leaf, node
    '''

    # Build ILP model and solve
    from phyilpd.gurobi_solver import gurobi_solver
    all_solutions = gurobi_solver(curr_node_to_leaves, curr_leaves, curr_list_of_ancestral_node, curr_node_to_mean_pwdist, wcl, cs, params.solver_verbose)

    if all_solutions == 'na':
        # continue to next parameter set if no solution
        print ('\nProgram EXIT: No clustering solution found.\n')
        exit(1)

    elif len(all_solutions) > 1:
        print ('\nMultiple ({}) solutions found..'.format(len(all_solutions)))

    # analyse solution and print outputs
    for sol_index, curr_taxon_to_clusterid in all_solutions.items():
        # failed integrality
        if curr_taxon_to_clusterid == False:
            if len(all_solutions) == 1:
                print ('\nProgram EXIT: No optimal solution found (Solver failed integrality).\n')
                exit(1)
            else:
                print ('\nWARNING: No optimal solution found for solution nr. {} (Solver failed integrality).'.format(sol_index))
                continue

        curr_outfname = 'phydelity_wcl{}_pct{}_time{}_{}'.format(str(wcl), str(params.coal_percentage), str(params.coal_time),treefname)

        curr_clusterid_to_taxa = {}
        for taxon, clusterid in curr_taxon_to_clusterid.items():

            """
            taxon_name = global_leaf_node_id_to_leafname[taxon]
            if re.search("IMH0028", taxon_name):
                print taxon_name, clusterid
            """

            try:
                curr_clusterid_to_taxa[clusterid].append(taxon)
            except:
                curr_clusterid_to_taxa[clusterid] = [taxon]

        print ('\nCleaning up clusters...')

        cleanup_object = clean_up_modules(curr_node_to_descendant_nodes, global_node_to_leaves, curr_node_to_leaves, wcl, cs, global_leaf_dist_to_node, global_leaf_to_ancestors, global_node_to_parent_node, global_nodepair_to_dist, global_node_to_child_leaves, global_leaf_node_id_to_leafname)

        # ensure that there are no weird odd leaves that are further away from everyone else
        #curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.remove_odd_leaf(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

        # ensure that the most descendant-possible node-id is subtending each cluster
        curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.transmission_cleanup(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

        # leave-one-out clean up for clusters violating wcl
        curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.loo_wcl_violation(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

        # determine distribution of clusters
        curr_clusterlen_distribution = get_cluster_size_distribution(curr_clusterid_to_taxa)

        # remove cluster-size sensitivity-induced subclusters
        curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.subsume_subclusters_under_x_percentile(curr_clusterid_to_taxa, curr_taxon_to_clusterid, curr_clusterlen_distribution, 25)

        # de-cluster outliers of descendent cluster that got clustered in ancestor cluster
        curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.decluster_outlying_taxa_clustered_to_anc_clusters(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

        # delete any clusters < cs
        curr_clusterid_to_taxa, curr_taxon_to_clusterid = cleanup_object.remove_clusters_below_cs(curr_clusterid_to_taxa, curr_taxon_to_clusterid)

        # print outputs
        print ('\nWriting outputs...')

        output_obj = phydelity_output(original_tree_string, global_leaf_node_id_to_leafname, curr_taxon_to_clusterid, curr_clusterid_to_taxa, curr_outfname)

        # cluster file
        curr_modified_tree_string = output_obj.cluster_output()

        # figtree annotated tree file
        output_obj.figtree_output(curr_modified_tree_string)

        # output pdf tree
        if params.pdf_tree:
            output_obj.ete3_pdf_tree_output()


    print ('\n...done.\n')
