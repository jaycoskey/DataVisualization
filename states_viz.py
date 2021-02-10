#!/usr/bin/env python3

import argparse
from collections import defaultdict
from copy import deepcopy
import networkx as nx
import sys


GRAY_FILL = 'style=filled fillcolor=lightgray'

STATES_CONN_NODE  = [('AZ','CO'), ('NM','UT')]
STATES_WEST_COAST = ['WA', 'OR', 'CA']

STATES_DOTFILE_GRIDIFY = 'states_gridify.dot'
STATES_DOTFILE_PLANNED = 'states_planned.dot'
STATES_DOTFILE_SPRINGS = 'states_springs.dot'

STATES_FILE_COORDS  = 'states.coords'
STATES_FILE_EDGES   = 'states.edges'

# Note: Use caution when crossing International Date Line
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


# Write to named dotfile after aligning lat/long coords to grid
# Start with the original grid, and find the x and y ordinals for each of the nodes.
# Repeatedly modify the (originally ordinal) coordinates by merging x-coords or y-coords.
# Finish and return modified graph when reaching a stopping condition.
#
# Details:
#   * At each step, consider the percentage gap between x coordinates and y coordinates,
#     and merge the two (x or y) coordinates with the smallest percentage gap,
#     among those merges that would not cause a constraint violation.
#   * Constraints:
#       - No two nodes share the same x & y ordinal coords.
#       - Graph is planar (motivated by example of adjacency in a geographical map)
#       - Oblique edges (neither horiz nor vert) match N-S and W-E ordering
#       - x ordinals start at 0, and form an integral set without gaps. Same for y.
# # Note:
#   * This is a proof of concept with 50 nodes. Efficiency and conciseness are not priorites.
#
def write_dotfile_gridify(g, dotfile, verbose=False):
    name2attrs = {name:attrs for name,attrs in g.nodes(data=True)}
    s2x = lambda s: name2attrs[s]['x']
    s2y = lambda s: name2attrs[s]['y']

    def avg(items):
        return sum(items) / len(items)

    def get_orddict(d, val2sortkey):
        result = { name:ordinal
                   for (ordinal, name) in enumerate(sorted(d.keys(), key=lambda k: val2sortkey(d[k])))
                 }
        return result

    def gridify(g, maxiter=100):
        # ----------------------------------------
        xs = [attrs['x'] for attrs in name2attrs.values()]
        xspan = max(xs) - min(xs)
        # ----------
        ys = [attrs['y'] for attrs in name2attrs.values()]
        yspan = max(ys) - min(ys)
        # ----------------------------------------

        s2xord = get_orddict(name2attrs, lambda attrs: attrs['x'])
        s2yord = get_orddict(name2attrs, lambda attrs: attrs['y'])

        gord = deepcopy(g)
        for s in name2attrs.keys():
            nx.set_node_attributes(gord, values={s:s2xord[s]}, name='x')
            nx.set_node_attributes(gord, values={s:s2yord[s]}, name='y')

        itercount = 0
        while True:
            if itercount >= maxiter:
                break
            itercount += 1
            verbose and print(f'# ========================================')
            verbose and print(f'INFO: Starting iteration {itercount}')

            # ----------------------------------------
            xord2ss = invert_dict(s2xord)
            xords = sorted(xord2ss.keys())
            xord2xavg = {xord:avg(list(map(s2x, xord2ss[xord]))) for xord in xords}
            # ----------
            yord2ss = invert_dict(s2yord)
            yords = sorted(yord2ss.keys())
            yord2yavg = {yord:avg(list(map(s2y, yord2ss[yord]))) for yord in yords}
            # ----------------------------------------

            # Percentage Longitude/latitude gap between ordinal keys i & i+1
            xord_pgaps = { ('x', xord) : 100.0 * (xord2xavg[xord+1] - xord2xavg[xord]) / xspan
                           for xord in range(max(xords))
                         }
            yord_pgaps = { ('y', yord) : 100.0 * (yord2yavg[yord+1] - yord2yavg[yord]) / yspan
                           for yord in range(max(yords))
                         }

            pgaps = deepcopy(xord_pgaps)
            pgaps.update(yord_pgaps)
            pgaps = sorted(pgaps.items(), key=lambda item: item[1])

            cand_count = 0
            for pgap_candidate in pgaps: 
                cand_count += 1
                # Unify the closest pair of ordinal coords.
                # TODO: Implement constraint-checking.
                # TODO: Use a temporary graph as a candidate. Commit changes on constraint failure.
                # TODO: When a merging of ordinal coords would cause a constraint failure,
                # TODO:   search down the prioritized list for a suitable candidate.
                pgap_axis = pgap_candidate[0][0]
                pgap_ord  = pgap_candidate[0][1]
                pgap_val  = pgap_candidate[1]
    
                max_ord = max(xords) if pgap_axis == 'x' else max(yords)
    
                if verbose:
                    print(f'INFO: Gap is smallest @ {pgap_axis}ord={pgap_ord}: pgap={pgap_val:.4f}')
                    if pgap_axis == 'x':
                        print(f'INFO:   --> Between {xord2ss[pgap_ord]} and {xord2ss[pgap_ord + 1]}')
                    else:
                        print(f'INFO:   --> Between {yord2ss[pgap_ord]} and {yord2ss[pgap_ord + 1]}')
    
                if pgap_axis == 'x':
                    # Unify two ordinal coords, and shift in more distal ones
                    for ordk in range(pgap_ord + 1, max_ord + 1):
                        for s in xord2ss[ordk]:
                            # verbose and print(f'INFO: Shifting west @ xord={ordk}: {s}')
                            s2xord[s] -= 1
                            # gord.nodes[s]['x'] = s2xord[s]  # Not supported. Use next line instead.
                            nx.set_node_attributes(gord, values={s:s2xord[s]}, name='x')
                else:
                    # Unify two ordinal coords, and shift in more distal ones
                    for ordk in range(pgap_ord + 1, max_ord + 1):
                        for s in yord2ss[ordk]:
                            # verbose and print(f'INFO: Shifting south @ yord={ordk}: {s}')
                            s2yord[s] -= 1
                            # gord.nodes[s]['y'] = s2yord[s]  # Not supported. Use next line instead.
                            nx.set_node_attributes(gord, values={s:s2yord[s]}, name='y')
                break

        verbose and print(f'# ========================================')
        verbose and print(f'INFO: Completed {itercount} iterations')

        nx.set_node_attributes(gord, values={'AK':int(gord.nodes(data=True)['CA']['x'])}, name='x')
        nx.set_node_attributes(gord, values={'AK':int(gord.nodes(data=True)['TX']['y'])}, name='y')

        nx.set_node_attributes(gord, values={'HI':int(gord.nodes(data=True)['AZ']['x'])}, name='x')
        nx.set_node_attributes(gord, values={'HI':int(gord.nodes(data=True)['TX']['y'])}, name='y')

        return gord

    def invert_dict(d):
        result = defaultdict(list)
        for (k,v) in d.items():
            result[v].append(k)
        return result

    def xpos_avg(names):
        return avg(lambda name: name2attrs[name]['x'], names)

    def ypos_avg(names):
        return avg(lambda name: name2attrs[name]['y'], names)

    gridified = gridify(g, maxiter=50)
    write_dotfile_planned(gridified, dotfile)


# Write to named dotfile using coordinates in graph argument
#
def write_dotfile_planned(g, dotfile):
    def write_nodes(f):
        for name, attrs in sorted(g.nodes(data=True)):
            scale_factor = 50
            x = attrs['x']
            y = attrs['y']
            pos_spec = f'pos="{scale_factor * x:.1f},{scale_factor * y:.1f}"'
            color_spec = f' {GRAY_FILL}' if name in ['AK', 'HI'] else ''
            writeln(f, 1, f'{name} [{pos_spec}{color_spec}]')

    def write_edges(f):
        for a, b in sorted(g.edges(data=False)):
            if b > a:
                is_conn_node = (a, b) in STATES_CONN_NODE
                color = 'color=red'
                style = 'style=dashed'
                attrs = f'[{color} {style}]' if is_conn_node else ''
                writeln(f, 1, f'{a} -- {b} {attrs}')

    with open(dotfile, 'w') as f:
        writeln(f, 0, 'strict graph States {')
        writeln(f, 1, 'node [fixedsize=true fontsize=10 height=0.45 width=0.45]')
        writeln(f)
        write_nodes(f)
        writeln(f)
        write_edges(f)
        writeln(f, 0, '}')


# Write to named dotfile using Graphviz layout algorithm
#
def write_dotfile_springs(g, dotfile):
    name2attrs = {name:attrs for name,attrs in g.nodes(data=True)}

    def write_nodes(f):
        for name, attrs in sorted(g.nodes(data=True)):
            x = attrs['x']
            y = attrs['y']
            color_attrs = ''
            if name in ['AK', 'HI']:
                color_attrs = f' {GRAY_FILL}'
                if name == 'AK':
                    x = name2attrs['CA']['x']
                    y = name2attrs['TX']['y']
                elif name == 'HI':
                    x = name2attrs['NM']['x']
                    y = name2attrs['TX']['y']
            pos_attrs = f'pos="{x:.2f},{y:.2f}"'
            writeln(f, 1, f'{name} [{pos_attrs}{color_attrs}]')

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

    with open(dotfile, 'w') as f:
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

def make_dotfile_gridify():
    src = read_graph( nodes_file=STATES_FILE_LATLONG
                    , edges_file=STATES_FILE_EDGES
                    , order='yx'
                    )
    test_graph(src)
    write_dotfile_gridify( g=src
                         , dotfile=STATES_DOTFILE_GRIDIFY
                         , verbose=False
                         )

def make_dotfile_planned():
    src = read_graph( nodes_file=STATES_FILE_COORDS
                    , edges_file=STATES_FILE_EDGES
                    )
    test_graph(src)
    write_dotfile_planned( g=src
                         , dotfile=STATES_DOTFILE_PLANNED
                         )

def make_dotfile_springs():
    src = read_graph( nodes_file=STATES_FILE_LATLONG
                    , edges_file=STATES_FILE_EDGES
                    , order='yx'
                    )
    test_graph(src)
    write_dotfile_springs( g=src
                         , dotfile=STATES_DOTFILE_SPRINGS
                         )


if __name__ == '__main__':
    prog='states_viz'
    parser = argparse.ArgumentParser( prog=prog
                                    , description='Visualize layout of US states'
                                    )
    parser.add_argument('-g', '--gridify', action='store_true')
    parser.add_argument('-p', '--planned', action='store_true')
    parser.add_argument('-s', '--springs', action='store_true')
    args = parser.parse_args(sys.argv[1:])

    if (not args.gridify) and (not args.planned) and (not args.springs):
        print(f'{prog}: At least one output flag required.')
        parser.print_help(sys.stderr)
        sys.exit(1)
    if args.gridify:
        make_dotfile_gridify()
    if args.planned:
        make_dotfile_planned()
    if args.springs:
        make_dotfile_springs()
