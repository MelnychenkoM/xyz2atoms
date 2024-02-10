import pandas as pd
import numpy as np

def find_cycles_in_graph(graph):
    cycles = []

    def find_new_cycles(path):
        start_node = path[0]
        sub = []

        for neighbor in graph[start_node]:
            if neighbor == start_node:
                continue
            if not visited(neighbor, path):
                sub = [neighbor]
                sub.extend(path)
                find_new_cycles(sub)
            elif len(path) > 2 and neighbor == path[-1]:
                p = rotate_to_smallest(path)
                inv = invert(p)
                if is_new(p) and is_new(inv):
                    cycles.append(p)

    def invert(path):
        return rotate_to_smallest(path[::-1])

    def rotate_to_smallest(path):
        n = path.index(min(path))
        return path[n:] + path[:n]

    def is_new(path):
        return path not in cycles

    def visited(node, path):
        return node in path

    for node in graph:
        find_new_cycles([node])
    return cycles

def find_start_sugar(adj_list_sugar):
    start_node = None
    next_node = None
    for key, value in adj_list_sugar.items():
        if key.startswith('O'):
            start_node = key
        
        elif key.startswith('C') and value.str.contains('N').any():
            next_node = key
    
    return start_node, next_node

def find_second_N_purines(adj_list_2base):
    for key, value in adj_list_2base.items():
        if key[0] == 'N':
            if value.str.contains('C').sum() == 2:
                return key

def sum_of_nitrogens(naming):
    nitrogen_sum = 0
    for atom, name in naming.items():
        if atom[0] == 'N':
            nitrogen_sum += int(name[1:])
    return nitrogen_sum

def choose_best_naming(naming_list):
    min_sum = float('inf')
    best_naming = None
    for naming in naming_list:
        current_sum = sum_of_nitrogens(naming)
        if current_sum < min_sum:
            min_sum = current_sum
            best_naming = naming
    return best_naming

cpk_colors = dict(
    C="black",
    H="white",
    N="blue",
    O="red",
    P="orange"
)
cpk_color_rest = "pink"

