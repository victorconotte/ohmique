#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Victor Conotte"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Victor Conotte"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Victor Conotte"
__email__ = "victorconotte@live.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file, "r") as fic:
        for nb_line, n in enumerate(fic):
            if nb_line % 4 == 1:
                yield n.strip()


def cut_kmer(read, kmer_size):
    decal = 0
    seq_length = len(read)
    while decal + kmer_size <= seq_length:
        yield read[decal:(decal + kmer_size)]
        decal += 1


def build_kmer_dict(fastq_file, kmer_size):
    reads = read_fastq(fastq_file)
    kmer_count = {}
    for read in reads:
        kmers = cut_kmer(read, kmer_size)
        for kmer in kmers:
            if not kmer in kmer_count:
                kmer_count[kmer] = 1
            else:
                kmer_count[kmer] += 1
    return kmer_count


def build_graph(kmer_dict):
    graph = nx.DiGraph()
    for km in kmer_dict.keys():
        p = km[:-1]
        s = km[1:]
        graph.add_edge(p, s, weight=kmer_dict[km])
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        graph.remove_nodes_from(path[(not delete_entry_node):(None if delete_sink_node else -1)])
    return graph

def std(data):
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    best_weight_indexes = [i for i, weight in enumerate(weight_avg_list)
                           if weight == max(weight_avg_list)]
    best_length_and_weights = [length for i, length in enumerate(path_length)
                               if i in best_weight_indexes]
    best_path_indexes = [i for i in best_weight_indexes
                         if path_length[i] == max(best_length_and_weights)]
    best_path_index = random.choice(best_path_indexes)
    graph = remove_paths(graph, path_list[:best_path_index]+path_list[(best_path_index+1):],
                         delete_entry_node, delete_sink_node)
    return graph

def path_average_weight(graph, path):
    total = 0
    for node_1, node_2 in zip(path[:-1], path[1:]):
            total += graph[node_1][node_2]["weight"]
    return total/(len(path)-1) if total else 0

def solve_bubble(graph, ancestor_node, descendant_node):
    path_list = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    path_lengths = [len(path_list) for path in path_list]
    avg_path_weights = [path_average_weight(graph, path) for path in path_list]

    return select_best_path(graph, path_list, path_lengths, avg_path_weights)

def simplify_bubbles(graph):
    starting_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)
    for ancestor_node in starting_nodes:
        for descendant_node in sink_nodes:
            current_entry = ancestor_node
            current_exit = descendant_node
            successors = list(graph.successors(current_entry))
            predecessors = list(graph.predecessors(current_exit))
            while len(successors) < 2 and successors:
                current_entry = successors[0]
                successors = list(graph.successors(current_entry))
            while len(predecessors) < 2 and predecessors:
                current_exit = predecessors[0]
                predecessors = list(graph.predecessors(current_exit))
            if list(nx.all_simple_paths(graph, current_entry, current_exit)):
                graph = solve_bubble(graph, current_entry, current_exit)
    return graph

def solve_entry_tips(graph, starting_nodes):
    graph = simplify_bubbles(graph)
    tips = []
    for start_node in starting_nodes:
        current_node = start_node
        path = [current_node]
        successors = list(graph.successors(current_node))
        predecessors = list(graph.predecessors(current_node))
        while len(successors) < 2 and len(predecessors) < 2 and successors:
            current_node = successors[0]
            path.append(current_node)
            successors = list(graph.successors(current_node))
            predecessors = list(graph.predecessors(current_node))
        tips.append(path)
    path_lengths = [len(path) for path in tips]
    avg_path_weights = [path_average_weight(graph, path) for path in tips]
    graph = select_best_path(graph, tips, path_lengths, avg_path_weights,
                             delete_entry_node=True)

    return graph

def solve_out_tips(graph, ending_nodes):
    graph = simplify_bubbles(graph)
    tips = []
    for sink_node in ending_nodes:
        current_node = sink_node
        path = [current_node]
        successors = list(graph.successors(current_node))
        predecessors = list(graph.predecessors(current_node))
        while len(successors) < 2 and len(predecessors) < 2 and predecessors:
            current_node = predecessors[0]
            path.append(current_node)
            successors = list(graph.successors(current_node))
            predecessors = list(graph.predecessors(current_node))
        tips.append(path[::-1])

    path_lengths = [len(path) for path in tips]
    avg_path_weights = [path_average_weight(graph, path) for path in tips]
    graph = select_best_path(graph, tips, path_lengths, avg_path_weights,
                             delete_sink_node=True)
    return graph

def get_starting_nodes(graph):
    entry_node_list = []
    for node in graph.nodes():
        if not list(graph.predecessors(node)):
            entry_node_list.append(node)
    return entry_node_list

def get_sink_nodes(graph):
    exit_node_list = []
    for node in graph.nodes():
        if not list(graph.successors(node)):
            exit_node_list.append(node)
    return exit_node_list

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    for start_node in starting_nodes:
        for sink_node in ending_nodes:
            for path in nx.all_simple_paths(graph, start_node, sink_node):
                contig = path[0]
                for i in range(1, len(path)):
                    contig += path[i][-1]
                contigs.append((contig, len(contig)))
    return contigs

def save_contigs(contigs_list, output_file):
   entt = ">contig_{} len={}\n"
   with open(output_file, "w") as fichiersortie:
        for i, contig in enumerate(contigs_list):
            fichiersortie.write(entt.format(i, contig[1]))
            fichiersortie.write(fill(contig[0])+ "\n")


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    print("test")
    graph = build_graph(build_kmer_dict(args.fastq_file, args.kmer_size))
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)
    
    graph=simplify_bubbles(graph)
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    
    graph=solve_entry_tips(graph, starting_nodes)
    graph=solve_out_tips(graph, ending_nodes)
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    contigs = get_contigs(graph,starting_nodes,ending_nodes)
    save_contigs(contigs, args.output_file)
    print(starting_nodes)
    print(ending_nodes)
    if args.graphimg_file:
        draw_graph(graph, args.graphimg_file)

if __name__ == '__main__':
    main()
