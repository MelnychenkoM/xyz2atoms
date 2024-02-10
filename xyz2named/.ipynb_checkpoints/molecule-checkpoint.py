import numpy as np
import pandas as pd
from xyz2named.tools import *


covalent_radius = {
    "N": 71.0,
    "O": 63.0,
    "C": 75.0,
    "H": 32.0,
    "P": 111.0
}

class AtomNames:
    
    """Contains methods necessary to cycle through graphs and name atoms"""
    
    def _name_sugar(self, graph, start_node, next_node):
        visited = {node: False for node in graph}
        carbon_index = [1]
    
        def name_sugar_dfs(current_node):
            stack = [current_node]
    
            while stack:
                current_node = stack.pop()
                if not visited[current_node]:
                    visited[current_node] = True
                    if current_node.startswith('O'):
                        self.named_list[self.elements.index(current_node)] = "O4'"
                    elif current_node.startswith('C'):
                        self.named_list[self.elements.index(current_node)] = f"C{carbon_index[0]}'"
                        carbon_index[0] += 1
                    elif current_node == next_node:
                        self.named_list[self.elements.index(current_node)] = "C1'"
    
                neighbors = [neighbor for neighbor in graph[current_node] if neighbor in graph]
    
                for neighbor in neighbors:
                    if neighbor == next_node and carbon_index[0] == 1:
                        stack.append(neighbor)
                        break
                    elif not visited[neighbor]:
                        stack.append(neighbor)

        
        name_sugar_dfs(start_node)

        # Name the only carbon left, which is C5'

        for idx, value in self.adj_list_indexes.items():
            if not self.named_list[idx] and self.elements[idx].startswith('C'):
                self.named_list[idx] = "C5'"
                break

        return 

    def _name_purines(self, base_graph, start_index=1):
        visited = {}
        atom_names = {}
        index = [start_index]
    
        def name_base_dfs(start_node, next_node):
            stack = [start_node]
    
            while stack:
                current_node = stack.pop()
                if not visited.get(current_node, False):
                    visited[current_node] = True
                    if current_node.startswith('N'):
                        atom_names[current_node] = f"N{index[0]}"
                        index[0] += 1
                    elif current_node.startswith('C'):
                        atom_names[current_node] = f"C{index[0]}"
                        index[0] += 1
                neighbors = [neighbor for neighbor in base_graph[current_node] if neighbor in base_graph]
                for neighbor in neighbors:
                    if neighbor == next_node:
                        stack.append(neighbor)
                    elif not visited.get(neighbor, False):
                        stack.append(neighbor)
    
        temp_list = []
    
        for start_node in base_graph.keys():
            for next_node in base_graph[start_node]:
                name_base_dfs(start_node, next_node)
                temp_list.append(atom_names.copy())
                visited.clear()
                atom_names.clear()
                index[0] = start_index

        best_naming = choose_best_naming(temp_list)
        self._update_names(best_naming)

    def _traverse_cycle(self, graph, start_node, direction='forward', start_index=0):
    
        visited = set()
        stack = [start_node]
        index = start_index - 1
        atom_names = {}
    
        while stack:
            current_node = stack.pop()
            if current_node not in visited:
                visited.add(current_node)
                index += 1
    
                neighbors = [neighbor for neighbor in graph[current_node] if neighbor in graph]
                
                if current_node.startswith('N'):
                    atom_names[current_node] = f"N{index}"
                elif current_node.startswith('C'):
                    atom_names[current_node] = f"C{index}"
    
                if direction == 'forward':
                    neighbors = neighbors[::-1]
                for neighbor in neighbors:
                    if neighbor not in visited:
                        stack.append(neighbor)
        return atom_names

    def _name_purines_2nd(self, graph):
        
        directions = ["forward", "reverse"]
        start_purine_node = find_second_N_purines(graph)
        
        temp_list = []
        
        for direction in directions:
            temp_list.append(self._traverse_cycle(graph, 
                                                  start_purine_node, 
                                                  direction, 
                                                  start_index=7))
            

        best_naming = choose_best_naming(temp_list)
        self._update_names(best_naming)

    def _name_pyrimidines(self, graph):
        directions = ["forward", "reverse"]
        start_nitrogen = find_second_N_purines(graph)
        
        temp_list = []
        
        for direction in directions:
            temp_list.append(self._traverse_cycle(graph, 
                                                  start_nitrogen, 
                                                  direction, 
                                                  start_index=1))
            

        best_naming = choose_best_naming(temp_list)
        self._update_names(best_naming)

    def _update_names(self, name_dict):
        for key, value in name_dict.items():
            if not self.named_list[self.elements.index(key)]:
                self.named_list[self.elements.index(key)] = value

class Molecule(AtomNames):
    def __init__(self):
        self.x = []
        self.y = []
        self.z = []
        self.elements = []
        self.adj_matrix = None
        self.named_list = None
        self.adj_list_series = {}
        self.adj_list_indexes = {}
        self.bond_lengths = {}
    
    def read_xyz(self, file_path):
        """Reads xyz file"""
        with open(file_path, 'r') as fl:
            for number, line in enumerate(fl.readlines()):
                result = line.split()
                self.x.append(float(result[1]))
                self.y.append(float(result[2]))
                self.z.append(float(result[3]))
                self.elements.append(str(result[0]) + f"{number}")
                
        self._build_adj_matrix()
        self._determine_molecule()
    
    def _build_adj_matrix(self):
        """Builds connection matrix"""
        xyz = np.stack((self.x, self.y, self.z), axis=-1)
        distances = np.linalg.norm(xyz[:, np.newaxis, :] - xyz[np.newaxis, :, :], axis=-1)
        
        covalent_radii = np.array([covalent_radius[atom[0]] for atom in self.elements])
        bond_distances = ((covalent_radii[:, np.newaxis] + covalent_radii) / 100) + 0.1
        
        self.adj_matrix = np.logical_and(0.1 < distances, distances < bond_distances).astype(int)
        graph_elements = {}
        
        for i, j in zip(*np.nonzero(self.adj_matrix)):
            graph_elements.setdefault(self.elements[i], []).append(self.elements[j])
            self.adj_list_indexes.setdefault(i, []).append(j)
            self.bond_lengths[frozenset([i, j])] = round(bond_distances[i, j], 5)

        self.adj_list_series = {key: pd.Series(value) for key, value in graph_elements.items()}

    def _determine_molecule(self):

        self.named_list = ['H' if element[0] == 'H' else '' for element in self.elements]
        
        cycles = find_cycles_in_graph(self.adj_list_indexes)

        purines = False
        if len(cycles) > 2:
            purines = True

        temp = [pd.Series([self.elements[id] for id in cycle_id]) for cycle_id in cycles]

        for cycle in temp:
            if len(cycle) == 5 and not cycle.str.contains('N').any():
                cycle = {k: v for k, v in self.adj_list_series.items() if k in cycle.values}
                
                start_sugar_node, next_sugar_node = find_start_sugar(cycle)
                self._name_sugar(cycle, start_sugar_node, next_sugar_node)

            elif purines:
                if len(cycle) == 6:
                    cycle = {k: v for k, v in self.adj_list_series.items() if k in cycle.values}
                    self._name_purines(cycle)

                elif len(cycle) == 5:
                    cycle = {k: v for k, v in self.adj_list_series.items() if k in cycle.values}
                    self._name_purines_2nd(cycle)

            else:
                if len(cycle) == 6:
                    cycle = {k: v for k, v in self.adj_list_series.items() if k in cycle.values}
                    self._name_pyrimidines(cycle)


        for key, values in self.adj_list_series.items():
            if key.startswith('C'):
                for atom in values:
                    if not self.named_list[int(atom[1:])]:
                        old_atom_index = int(atom[1:])             
                        self.named_list[old_atom_index] = atom[0] + self.named_list[int(key[1:])][1:]        

        #print(self.named_list)
    
    def __len__(self):
        return len(self.elements)

    def __getitem__(self, pos):
        if self.named_list[pos] == "":
            return self.elements[pos][0], (self.x[pos], self.y[pos], self.z[pos])
        else:
            return self.named_list[pos], (self.x[pos], self.y[pos], self.z[pos])

