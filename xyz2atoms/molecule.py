import numpy as np
import pandas as pd
import os
from xyz2atoms.tools import *


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

    def _name_purines(self, base_graph):
        
        directions = ["forward", "reverse"]
        start_purine_node = self._find_C4_purines(base_graph)

        temp_list = []
        
        for direction in directions:
            temp_list.append(self._traverse_cycle(base_graph, 
                                                  start_purine_node, 
                                                  direction, 
                                                  start_index=3, idx=1, cycle=7))
            

        best_naming = choose_best_naming(temp_list)
        self._update_names(best_naming, replace=True)

    def _traverse_cycle(self, graph, start_node, direction='forward', start_index=0, idx=1, cycle=100):
    
        visited = set()
        stack = [start_node]
        index = start_index
        atom_names = {}
    
        while stack:
            current_node = stack.pop()
            if current_node not in visited:
                visited.add(current_node)

                if index == cycle - 1:
                    index += 1

                index = (index + idx) % cycle

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
                                                  start_index=10, idx=-1))
        

        best_naming = choose_best_naming(temp_list, return_smallest=False)
        self._update_names(best_naming)

    def _name_pyrimidines(self, graph):
        directions = ["forward", "reverse"]
        start_nitrogen = find_second_N_purines(graph)
        
        temp_list = []
        
        for direction in directions:
            temp_list.append(self._traverse_cycle(graph, 
                                                  start_nitrogen, 
                                                  direction, 
                                                  start_index=0))
            

        best_naming = choose_best_naming(temp_list)
        self._update_names(best_naming)

    def _find_C4_purines(self, base_graph):
        for key, value in base_graph.items():
            if value.str.contains('N').sum() == 2:
                if self.elements[self.named_list.index('N9')] in value.to_list():
                    return key
                
    def _get_sugar_angles(self):
        atom_sets = {
            "v0": ["C4'", "O4'", "C1'", "C2'"],
            "v1": ["O4'", "C1'", "C2'", "C3'"],
            "v2": ["C1'", "C2'", "C3'", "C4'"],
            "v3": ["C2'", "O3'", "C4'", "O4'"],
            "v4": ["C3'", "C4'", "O4'", "C1'"]
        }

        self.sugar_torsion_angles = {}

        for key, atom_set in atom_sets.items():
            torsion_angle_value = torsion_angle(*[self[atom_name] for atom_name in atom_set])
            self.sugar_torsion_angles[key] = torsion_angle_value

    def _update_names(self, name_dict, replace=False):
        for key, value in name_dict.items():
            index = self.elements.index(key)
            if not self.named_list[index] or replace:
                self.named_list[index] = value


class Molecule(AtomNames):
    def __init__(self):
        self.x = []
        self.y = []
        self.z = []
        self.elements = []
        self.adj_matrix = None
        self.named_list = []
        self.adj_list_series = {}
        self.adj_list_indexes = {}
        self.bond_lengths = {}
        self.chain_identifier = []

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
        self._get_sugar_angles()

    def read_pdb(self, file_path, model=False):
        found_model = False
        with open(file_path, 'r', encoding='utf-8') as fl:
            for line in fl.readlines():
                if not found_model:
                    if line.startswith('MODEL') and int(line[10:].strip()) == model:
                        found_model = True
                else:
                    if line.startswith('ATOM'):
                        self.chain_identifier.append(float(line[25].strip()))
                        self.x.append(float(line[30:38].strip()))
                        self.y.append(float(line[38:46].strip()))
                        self.z.append(float(line[46:54].strip()))
                        self.named_list.append(line[11:16].strip())
                        self.elements.append(line[11:16].strip()[0])
                    
                    elif line.startswith('ENDMDL'):
                        self._build_adj_matrix()
                        break
        
    def write_pdb(self, file_name):
        with open(file_name, 'w') as fl:
            for i in range(len(self.elements)):
                fl.write(
                    f"ATOM  {i+1:>4}  {self.named_list[i]:<3}   A{self.named_list[i][0]:<3}"
                    f"{i+1:>4}    {self.x[i]:>8.3f}{self.y[i]:>8.3f}{self.z[i]:>8.3f}  "
                    f"1.00  0.00      {self.elements[i][0]:>2}\n"
                )
    
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
        cycles = sorted(cycles, key=len)

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

                if len(cycle) == 5:
                    cycle = {k: v for k, v in self.adj_list_series.items() if k in cycle.values}
                    self._name_purines_2nd(cycle)

                elif len(cycle) == 6:
                    cycle = {k: v for k, v in self.adj_list_series.items() if k in cycle.values}
                    self._name_purines(cycle)

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

    def __len__(self):
        return len(self.elements)

    def __getitem__(self, pos):
        if isinstance(pos, str):
            try:
                index = self.named_list.index(pos)
                return np.array([self.x[index], self.y[index], self.z[index]])
            except ValueError:
                raise KeyError(f"Atom with name {pos} is not found.")

        elif isinstance(pos, int):
            if self.named_list[pos] == "":
                return self.elements[pos][0], np.array([self.x[pos], self.y[pos], self.z[pos]])
            else:
                return self.named_list[pos], np.array([self.x[pos], self.y[pos], self.z[pos]])
        else:
            raise TypeError("Invalid index type.")
    
    def get_atom(self, name, chain=0):
        pass

    def get_distance(self, atom1, atom2):
        index1 = self.named_list.index(atom1)
        index2 = self.named_list.index(atom2)
        return self.bond_lengths.get(frozenset([index1, index2]), 0)


