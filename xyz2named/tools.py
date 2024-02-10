import pandas as pd
import numpy as np
import plotly.graph_objs as go

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
            if value.str.contains('C').sum() == 3:
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


def to_plotly_figure(graph) -> go.Figure:

    def atom_trace():
        """Creates an atom trace for the plot."""
        colors = [cpk_colors.get(element[0], cpk_color_rest) for element in graph.elements]
        markers = dict(
            color=colors,
            line=dict(color="lightgray", width=2),
            size=7,
            symbol="circle",
            opacity=0.8,
        )
        trace = go.Scatter3d(
            x=graph.x,
            y=graph.y,
            z=graph.z,
            mode="markers",
            marker=markers,
            text=graph.elements,
        )
        return trace

    def bond_trace():
        """ "Creates a bond trace for the plot."""
        trace = go.Scatter3d(
            x=[],
            y=[],
            z=[],
            hoverinfo="none",
            mode="lines",
            marker=dict(color="grey", size=7, opacity=1),
        )
        adjascent_atoms = (
            (atom, neighbour)
            for atom, neighbours in graph.adj_list_indexes.items()
            for neighbour in neighbours
        )
        for i, j in adjascent_atoms:
            trace["x"] += (graph.x[i], graph.x[j], None)
            trace["y"] += (graph.y[i], graph.y[j], None)
            trace["z"] += (graph.z[i], graph.z[j], None)
        return trace

    annotations_elements = [
        dict(text=element[0], x=x, y=y, z=z, showarrow=False, yshift=15)
        for element, (x, y, z) in graph
    ]
    annotations_named_elements = [
        dict(text=element, x=x, y=y, z=z, showarrow=False, yshift=15)
        for element, (x, y, z) in graph
    ]
    annotations_indices = [
        dict(text=number, x=x, y=y, z=z, showarrow=False, yshift=15)
        for number, (_, (x, y, z)) in enumerate(graph)
    ]

    annotations_bonds = []
    for (i, j), length in graph.bond_lengths.items():
        x = (graph.x[i] + graph.x[j]) / 2
        y = (graph.y[i] + graph.y[j]) / 2
        z = (graph.z[i] + graph.z[j]) / 2
        annotations_bonds.append(
            dict(
                text=round(length, 2),
                x=x,
                y=y,
                z=z,
                showarrow=False,
                yshift=15,
                font=dict(color="steelblue"),
            )
        )

    updatemenus = list(
        [
            dict(
                buttons=list(
                    [
                        dict(
                            label="Elements",
                            method="relayout",
                            args=[{"scene.annotations": annotations_elements}],
                        ),
                        dict(
                            label="Named Elements",
                            method="relayout",
                            args=[{"scene.annotations": annotations_named_elements}],
                        ),
                        dict(
                            label="Elements & Bond Lengths",
                            method="relayout",
                            args=[
                                {
                                    "scene.annotations": annotations_elements
                                    + annotations_bonds
                                }
                            ],
                        ),
                        dict(
                            label="Indices",
                            method="relayout",
                            args=[{"scene.annotations": annotations_indices}],
                        ),
                        dict(
                            label="Indices & Bond Lengths",
                            method="relayout",
                            args=[
                                {
                                    "scene.annotations": annotations_indices
                                    + annotations_bonds
                                }
                            ],
                        ),
                        dict(
                            label="Bond Lengths",
                            method="relayout",
                            args=[{"scene.annotations": annotations_bonds}],
                        ),
                        dict(
                            label="Hide All",
                            method="relayout",
                            args=[{"scene.annotations": []}],
                        ),
                    ]
                ),
                direction="down",
                xanchor="left",
                yanchor="top",
            ),
        ]
    )

    data = [atom_trace(), bond_trace()]
    axis_params = dict(
        showgrid=False,
        showbackground=False,
        showticklabels=False,
        zeroline=False,
        titlefont=dict(color="white"),
    )
    layout = dict(
        scene=dict(
            xaxis=axis_params,
            yaxis=axis_params,
            zaxis=axis_params,
            annotations=annotations_elements,
        ),
        margin=dict(r=0, l=0, b=0, t=0),
        showlegend=False,
        updatemenus=updatemenus,
    )
    figure = go.Figure(data=data, layout=layout)

    return figure