#!/usr/bin/env python3

import argparse
import networkx as nx
import sys


GRAY_ATTRS = 'fillcolor=lightgray style=filled'

STATES_CONN_NODE  = [('AZ','CO'), ('NM','UT')]
STATES_WEST_COAST = ['WA', 'OR', 'CA']

STATES_DOTFILE_PLANNED = 'states_planned.dot'
STATES_DOTFILE_SPRINGS = 'states_springs.dot'

STATES_FILE_COORDS  = 'states.coords'
STATES_FILE_EDGES   = 'states.edges'
STATES_FILE_LATLONG = 'states.latlong'


def read_graph(nodes_file, edges_file, order='xy'):
    def read_nodes(nodes_file):
        with open(nodes_file, 'r') as f:
            for line in f.readlines():
                items = line.strip().split('#')[0].split()
                if len(items) == 0:
                    continue
                assert(len(items) == 3)
                name, x, y = items
                if order == 'yx':
                    x, y = y, x
                yield (name, float(x), float(y))

    def read_edges(edges_file):
        with open(edges_file, 'r') as f:
            for line in f.readlines():
                vals = line.strip().split('#')[0].split(':')
                if len(vals) == 0:
                    continue
                assert(len(vals) == 2)
                name = vals[0]
                nbrs = vals[1].split() if len(vals) > 1 else []
                yield (name, nbrs)

    result = nx.Graph()
    for name, x, y in read_nodes(nodes_file):
        result.add_node(name, x=x, y=y)
    for node, nbrs in read_edges(edges_file):
        for nbr in nbrs:
            if nbr > node:
                result.add_edge(node, nbr)
    return result


def write_dotfile_planned(g, dotfile_name):
    def write_nodes(f):
        # writeln(f, 1, 'subgraph cluster_states {')
        for name, attrs in sorted(g.nodes(data=True)):
            scale_factor = 50
            x = attrs['x']
            y = attrs['y']
            pos_spec = f'pos="{scale_factor * x:.1f},{scale_factor * y:.1f}"'
            writeln(f, 1, f'{name} [{pos_spec}]')
        # writeln(f, 1, '}')

    def write_edges(f):
        for a, b in sorted(g.edges(data=False)):
            if b > a:
                is_conn_node = (a, b) in STATES_CONN_NODE
                color = 'color=red'
                style = 'style=dashed'
                attrs = f'[{color} {style}]' if is_conn_node else ''
                writeln(f, 1, f'{a} -- {b} {attrs}')

    with open(dotfile_name, 'w') as f:
        writeln(f, 0, 'strict graph States {')
        writeln(f, 1, 'node [fixedsize=true fontsize=10 height=0.45 width=0.45]')
        writeln(f)
        write_nodes(f)
        writeln(f)
        write_edges(f)
        writeln(f, 0, '}')


def write_dotfile_springs(g, dotfile_name):
    name2attrs = {name:attrs for name,attrs in g.nodes(data=True)}

    def write_nodes(f):
        for name, attrs in sorted(g.nodes(data=True)):
            if name == 'AK':
                x = name2attrs['CA']['x']
                y = name2attrs['TX']['y']
            elif name == 'HI':
                x = name2attrs['NM']['x']
                y = name2attrs['TX']['y']
            else:
                x = attrs['x']
                y = attrs['y']
            pos_attrs = f'[pos="{x:.2f},{y:.2f}"]'
            writeln(f, 1, f'{name} {pos_attrs}')

    def write_edges(f):
        writeln(f, 1, 'AK -- HI [style=invis]')  # Place nearby, but unlinked
        writeln(f, 1, 'AK -- CA [style=invis]')  # Tweak to place AK near CA
        writeln(f, 1, 'HI -- TX [style=invis]')  # Tweak to place HI near TX
        writeln(f, 1, 'NJ -- CT [style=invis]')  # Tweak to place NJ
        writeln(f, 1, 'NJ -- MD [style=invis]')  # Tweak to place NJ
        writeln(f)
        for a, b in sorted(g.edges(data=False)):
            if b > a:  # Avoid duplication
                is4 = (a, b) in STATES_CONN_NODE
                edge_attrs = '[color=red style=dashed]' if is4 else ''
                writeln(f, 1, f'{a} -- {b} {edge_attrs}')

    with open(dotfile_name, 'w') as f:
        writeln(f, 0, 'strict graph States {')
        writeln(f, 1, 'rankdir=LR')
        writeln(f, 1, 'splines=line')
        writeln(f, 1, 'node [fixedsize=true fontsize=10 height=0.45 width=0.45]')
        writeln(f)
        write_nodes(f)
        writeln(f)
        write_edges(f)
        writeln(f, 0, '}')


def writeln(f, indent_level=0, s=''):
    f.write(' ' * 4 * indent_level)
    f.write(s + '\n')

# ----------------------------------------

def test_graph(state_graph):
    assert(len(state_graph.nodes(data=False)) == 50)
    nbr_count = 2.0 * len(state_graph.edges(data=False))
    node_count = len(state_graph.nodes(data=False))
    avg_nbr_count = nbr_count / node_count
    # print(f'avg_nbr_count={avg_nbr_count}')
    assert(4.25 < avg_nbr_count < 4.3)

# ----------------------------------------

def make_dotfile_planned():
    src = read_graph( nodes_file=STATES_FILE_COORDS
                    , edges_file=STATES_FILE_EDGES
                    )
    test_graph(src)
    write_dotfile_planned( g=src
                         , dotfile_name=STATES_DOTFILE_PLANNED
                         )


def make_dotfile_springs():
    src = read_graph( nodes_file=STATES_FILE_LATLONG
                    , edges_file=STATES_FILE_EDGES
                    , order='yx'
                    )
    test_graph(src)
    write_dotfile_springs(g=src, dotfile_name=STATES_DOTFILE_SPRINGS)


if __name__ == '__main__':
    prog='states_viz'
    parser = argparse.ArgumentParser( prog=prog
                                    , description='Visualize layout of US states'
                                    )
    parser.add_argument('-p', '--planned', action='store_true')
    parser.add_argument('-s', '--springs', action='store_true')
    args = parser.parse_args(sys.argv[1:])

    if (not args.planned) and (not args.springs):
        print(f'{prog}: At least one output flag required.')
        parser.print_help(sys.stderr)
        sys.exit(1)
    if args.planned:
        make_dotfile_planned()
    if args.springs:
        make_dotfile_springs()
