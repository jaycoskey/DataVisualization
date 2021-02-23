#!/usr/bin/env python3
# Copyright Jay M. Coskey, 2021
# Note: This is a proof of concept. Efficiency and conciseness are not priorites.

import argparse
from collections import defaultdict
from copy import deepcopy
import networkx as nx
import sys
import unittest


GRAY_FILL = 'style=filled fillcolor=lightgray'
MAX_ITER_COUNT = 100
MAX_PGAP = 30.0

STATES_REPOSITIONS = [('AK', ('CA', 'TX')), ('HI', ('AZ', 'TX'))]
STATES_REPOSITIONS_GRIDIFY = [('NH', ('MA', 'MN'))]

STATES_DISCONNECTED = ['AK','HI']
STATES_FOUR_CORNERS = [('AZ','CO'), ('NM','UT')]

STATES_DOTFILE_GRIDIFY = 'states_gridify.dot'
STATES_DOTFILE_PLANNED = 'states_planned.dot'
STATES_DOTFILE_SPRINGS = 'states_springs.dot'

STATES_GRIDIFIED_COORDS = 'states_gridified.coords'
STATES_INFILE_COORDS    = 'states.coords'
STATES_INFILE_EDGES     = 'states.edges'

# Note: Use caution when crossing International Date Line
STATES_INFILE_LATLONG = 'states.latlong'


# Used to place the disconnected states: AK & HI
def customize_graph(g, repositions):
    for reposition in repositions:
        s = reposition[0]
        x_src_name, y_src_name = reposition[1]

        x_src = int(g.nodes(data=True)[x_src_name]['x'])
        y_src = int(g.nodes(data=True)[y_src_name]['y'])

        nx.set_node_attributes(g, values={s:x_src}, name='x')
        nx.set_node_attributes(g, values={s:y_src}, name='y')


def geom_has_general_intersection(p1, q1, p2, q2):
    or_p2 = geom_orientation(p1, q1, p2)
    or_q2 = geom_orientation(p1, q1, q2)
    or_p1 = geom_orientation(p2, q2, p1)
    or_q1 = geom_orientation(p2, q2, q1)

    if or_p2 * or_q2 == -1 and or_p1 * or_q1 == -1:
        return True

    return False


def geom_has_open_segment_violation(p1, q1, p2, q2):
    if (  geom_open_segment_ab_contain_c(p1, q1, p2)
       or geom_open_segment_ab_contain_c(p1, q1, q2)
       or geom_open_segment_ab_contain_c(p2, q2, p1)
       or geom_open_segment_ab_contain_c(p2, q2, q1)
       ):
        return True
    return False


def geom_box_ab_contain_c(a, b, c):
    if (   c[0] <= max(a[0], b[0])
       and c[0] >= min(a[0], b[0])
       and c[1] <= max(a[1], b[1])
       and c[1] >= min(a[1], b[1])
       ):
        return True
    return False


def geom_closed_segment_ab_contain_c(a, b, c):
    result = ( geom_orientation(a, b, c) == 0
             and geom_box_ab_contain_c(a, b, c)
             )
    return result


def geom_open_segment_ab_contain_c(a, b, c):
    result = (c != a) and (c != b) and geom_closed_segment_ab_contain_c(a, b, c)
    return result


def geom_orientation(a, b, c):
    """Returns 0 if points are colinear, 1 if clockwise, -1 if counter-clockwise
       Does so by computing determinant of matrix row vectors a, b, c, extended w/ 1.
    """
    det = (b[1] - a[1]) * (c[0] - b[0]) - (b[0] - a[0]) * (c[1] - b[1])
    return sgn(det)


def get_gord(g, s2xy):
    gord    = deepcopy(g)
    s2xord  = get_orddict(s2xy, lambda xy: xy['x'])
    s2yord  = get_orddict(s2xy, lambda xy: xy['y'])
    for s in s2xy.keys():
        nx.set_node_attributes(gord, values={s:s2xord[s]}, name='x')
        nx.set_node_attributes(gord, values={s:s2yord[s]}, name='y')
    return gord


def get_orddict(d, val2sortkey):
    result = { s:ordinal
               for (ordinal, s) in enumerate(sorted(d.keys(), key=lambda k: val2sortkey(d[k])))
             }
    return result


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


def sgn(x):
    return 1 if x > 0 else (-1 if x < 0 else 0)


def signature(g, do_check_and_update=True):
    xys = [(xy['x'], xy['y']) for s,xy in g.nodes(data=True)]
    xsum = sum([xy[0] for xy in xys])
    ysum = sum([xy[1] for xy in xys])
    if do_check_and_update and (xsum, ysum) == signature.prev:
        print(f'WARN: Graph signature equals last computed signature: ({xsum}, {ysum})')
    if do_check_and_update:
        signature.prev = (xsum, ysum)
    return xsum, ysum

signature.prev = None  # Could be implemented via decorator instead.


# Unify two ordinal coords, and shift in more distal ones
def squish(g, s2xord, s2yord, xord2ss, yord2ss, pgap_axis, pgap_ord, max_ord, verbose=False):
    xord_min = min(s2xord.values())
    xord_max = max(s2xord.values())
    yord_min = min(s2yord.values())
    yord_max = max(s2yord.values())

    if verbose:
        prefix = f'\tINFO: squish: Range=({xord_min}-{xord_max}, {yord_min}-{yord_max})'
        print(f'{prefix}: Shifting {pgap_axis} down from {pgap_ord + 1}')
    for ordk in range(pgap_ord + 1, max_ord + 1):
        if pgap_axis == 'x':
            new_xord = ordk - 1
            state_count = len(xord2ss[ordk])
            for s in xord2ss[ordk]:
                assert(g.nodes[s]["x"] == ordk)
                nx.set_node_attributes(g, values={s:new_xord}, name='x')
        else:
            new_yord = ordk - 1
            state_count = len(yord2ss[ordk])
            for s in yord2ss[ordk]:
                assert(g.nodes[s]["y"] == ordk)
                nx.set_node_attributes(g, values={s:new_yord}, name='y')


# Write to named dotfile after aligning lat/long coords to grid.
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
#       - TODO: Ensure that oblique edges (neither horiz nor vert) match N-S and W-E ordering
#       - TODO: Improve heuristics to make result more compact while preserving or improving aesthetics.
#       - TODO: Add flexibility to evoluationary algorithm.
#
def write_dotfile_gridify(g, dotfile, verbose=False):
    do_write_gridified = True

    def avg(items):
        return sum(items) / len(items)

    def get_ord_info(g, s2x, s2y):
        s2xyords = {s:xyords for s,xyords in g.nodes(data=True)}

        s2xord   = {s:xyords['x'] for (s, xyords) in s2xyords.items()}
        xord2ss  = invert_dict(s2xord)

        s2yord   = {s:xyords['y'] for (s, xyords) in s2xyords.items()}
        yord2ss  = invert_dict(s2yord)

        xords = sorted(set(xord2ss.keys()))
        yords = sorted(set(yord2ss.keys()))

        xord2xavg = {xord:avg(list(map(s2x, xord2ss[xord]))) for xord in xords}
        yord2yavg = {yord:avg(list(map(s2y, yord2ss[yord]))) for yord in yords}

        xs = [s2x(s) for s in s2xord.keys()]
        ys = [s2y(s) for s in s2yord.keys()]

        xspan = max(xs) - min(xs)
        yspan = max(ys) - min(ys)

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

        assert_ords(s2xord, s2yord, xord2ss, yord2ss)
        return s2xord, s2yord, xord2ss, yord2ss, pgaps

    def gridify(g, max_iter_count=100, verbose=False):
        s2xy = {s:xy for s,xy in g.nodes(data=True)}

        s2x = lambda s: s2xy[s]['x']
        s2y = lambda s: s2xy[s]['y']

        gord = get_gord(g, s2xy)

        iter_count = 0
        while True:
            if iter_count >= max_iter_count:
                break
            iter_count += 1
            verbose and print('='*20, f'INFO: Starting iteration {iter_count}')

            s2xord, s2yord, xord2ss, yord2ss, pgaps = get_ord_info(gord, s2x, s2y)
            pgap_candidates = list(filter(lambda pgap_cand: pgap_cand[1] <= MAX_PGAP, pgaps))

            # From all the candidates, unify the closest viable pair of (ordinal) coords.
            for cand_count in range(len(pgap_candidates)):
                print(' '*4, '-'*10, f'INFO: (iter,cand)=({iter_count},{cand_count})')
                pgap_candidate = pgap_candidates[cand_count]
                pgap_axis = pgap_candidate[0][0]
                pgap_ord  = pgap_candidate[0][1]
                pgap_val  = pgap_candidate[1]

                gcand = deepcopy(gord)
                s2xord, s2yord, xord2ss, yord2ss, pgaps = get_ord_info(gcand, s2x, s2y)
                max_ord = ( max([gcand.nodes[s]['x'] for s in s2xord.keys()]) if pgap_axis == 'x'
                            else max([gcand.nodes[s]['y'] for s in s2yord.keys()])
                          )

                if verbose:
                    print(f'\tINFO: gridify: Gap is smallest @ {pgap_axis}ord={pgap_ord}: pgap={pgap_val:.4f}')
                    if pgap_axis == 'x':
                        print(f'\tINFO:   --> Between {xord2ss[pgap_ord]} and {xord2ss[pgap_ord + 1]}')
                    else:
                        print(f'\tINFO:   --> Between {yord2ss[pgap_ord]} and {yord2ss[pgap_ord + 1]}')

                squish( gcand
                      , s2xord,  s2yord
                      , xord2ss, yord2ss
                      , pgap_axis, pgap_ord, max_ord, verbose
                      )
                prefix = f'\tINFO: gridify: {pgap_axis}ord={pgap_ord}'
                verbose and print(f'{prefix}: Post-squish: sig(gcand)={signature(gcand, True)}')

                is_cand_valid_nodes = is_valid_distinct_nodes(gcand)
                is_cand_valid_edges = is_valid_intersections(gcand)
                is_cand_valid = is_cand_valid_nodes and is_cand_valid_edges
                if is_cand_valid:
                    verbose and print(f'\tINFO: gridify: Updating graph')
                    gord = deepcopy(gcand)
                    break  # Using this candidate. Go to next iteration
                else:
                    if verbose:
                        if is_cand_valid_nodes:
                            print(f'\tINFO: gridify: Graph not updated: Candidate nodes invalid')
                        if is_cand_valid_edges:
                            print(f'\tINFO: gridify: Graph not updated: Candidate edges invalid')
            else:
                verbose and print(f'INFO: (iter,cand)=({iter_count},{cand_count}): Exhausted all candidates. Exiting iteration loop.')
                break  # Did not act on any candidates. Do not proceed to next iteration.

        verbose and print(f'# ========================================')
        verbose and print(f'INFO: Completed {iter_count} iterations')

        return gord

    def invert_dict(d):
        result = defaultdict(list)
        for (k,v) in d.items():
            result[v].append(k)
        return result

    def is_valid_distinct_nodes(g):
        coords = [(attrs['x'], attrs['y']) for _,attrs in g.nodes(data=True)]
        return len(set(coords)) == len(coords)

    def is_valid_intersections(g):
        for edgei in range(len(g.edges())):
            (ep1, eq1) = list(g.edges())[edgei]
            for edgej in range(edgei + 1, len(g.edges())):
                (ep2, eq2) = list(g.edges())[edgej]
                coords = { s:(xyords['x'], xyords['y'])
                           for s,xyords in g.nodes(data=True)
                           if s in [ep1, eq1, ep2, eq2]
                         }
                p1 = coords[ep1]
                q1 = coords[eq1]
                p2 = coords[ep2]
                q2 = coords[eq2]
                if (  geom_has_general_intersection(p1, q1, p2, q2)
                   or geom_has_open_segment_violation(p1, q1, p2, q2)
                   ):
                    print(f'\tINFO: Found violation @ {ep1}={p1},{eq1}={q1},{ep2}={p2},{eq2}={q2}')
                    return False
        return True

    def xpos_avg(s2xy, ss):
        return avg(lambda s: s2xy[s]['x'], ss)

    def ypos_avg(s2xy, ss):
        return avg(lambda s: s2xy[s]['y'], ss)

    gridified = gridify(g, MAX_ITER_COUNT, verbose)

    if do_write_gridified:
        with open(STATES_GRIDIFIED_COORDS, 'w') as f:
            writeln(f, 0, f'# Written by states_viz when creating "gridify" visualization.')
            s2xyords = {s:xyords for s,xyords in gridified.nodes(data=True)}
            for s, xyords in s2xyords.items():
                x = xyords["x"]
                y = xyords["y"]
                writeln(f, 0, f'{s} {x:>2} {y:>2}')
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
            color_spec = f' {GRAY_FILL}' if name in STATES_DISCONNECTED else ''
            writeln(f, 1, f'{name} [{pos_spec}{color_spec}]')

    def write_visible_edges(f):
        for a, b in sorted(g.edges(data=False)):
            if b > a:
                writeln(f, 1, f'{a} -- {b}')
        for (a,b) in STATES_FOUR_CORNERS:
            color = 'color=red'
            style = 'style=dashed'
            attrs = f'[{color} {style}]'
            writeln(f, 1, f'{a} -- {b} {attrs}')

    with open(dotfile, 'w') as f:
        writeln(f, 0, 'strict graph States {')
        writeln(f, 1, 'node [fixedsize=true fontsize=10 height=0.30 width=0.30]')
        writeln(f)
        write_nodes(f)
        writeln(f)
        write_visible_edges(f)
        writeln(f)
        writeln(f, 0, '}')


# Write to named dotfile using Graphviz layout algorithm
#
def write_dotfile_springs(g, dotfile):
    s2attrs = {s:attrs for s,attrs in g.nodes(data=True)}

    def write_nodes(f):
        for s,attrs in sorted(g.nodes(data=True)):
            x = attrs['x']
            y = attrs['y']
            pos_attrs   = f'pos="{x:.2f},{y:.2f}"'
            color_attrs = f' {GRAY_FILL}' if s in STATES_DISCONNECTED else ''
            writeln(f, 1, f'{s} [{pos_attrs}{color_attrs}]')

    def write_visible_edges(f):
        for a, b in sorted(g.edges(data=False)):
            if b > a:  # Avoid duplication
                writeln(f, 1, f'{a} -- {b}')
        for (a,b) in STATES_FOUR_CORNERS:
            edge_attrs = '[color=red style=dashed]'
            writeln(f, 1, f'{a} -- {b} {edge_attrs}')

    with open(dotfile, 'w') as f:
        writeln(f, 0, 'strict graph States {')
        writeln(f, 1, 'rankdir=LR')
        writeln(f, 1, 'splines=line')
        writeln(f, 1, 'node [fixedsize=true fontsize=10 height=0.45 width=0.45]')
        writeln(f)
        write_nodes(f)
        writeln(f)
        write_visible_edges(f)
        writeln(f)

        new_edges = [ ('AK', 'HI')  # Place near each other, but without visible edge'
                    , ('AK', 'CA')  # Place AK near CA
                    , ('HI', 'TX')  # Place HI near TX
                    , ('NJ', 'CT'), ('NJ', 'MD')  # Orient the northestern states properly
                    ]
        write_invisible_edges(f, new_edges)
        writeln(f, 0, '}')


def write_invisible_edges(f, edges):
    for edge in edges:
        edge_str = f'{edge[0]} -- {edge[1]} [style=invis]'
        writeln(f, 1, edge_str)


def writeln(f, indent_level=0, s=''):
    f.write(' ' * 4 * indent_level)
    f.write(s + '\n')

# ----------------------------------------

def assert_graph_sanity(state_graph):
    assert(len(state_graph.nodes(data=False)) == 50)
    nbr_count = 2.0 * len(state_graph.edges(data=False))
    node_count = len(state_graph.nodes(data=False))
    avg_nbr_count = nbr_count / node_count
    # print(f'avg_nbr_count={avg_nbr_count}')
    assert(4.15 < avg_nbr_count < 4.25)

def assert_ords(s2xord, s2yord, xord2ss, yord2ss):
    for s in s2xord.keys():
        assert(s in xord2ss[s2xord[s]])
        assert(s in yord2ss[s2yord[s]])
    for xord in set(s2xord.values()):
        for s in xord2ss[xord]:
            assert(xord == s2xord[s])
    for yord in set(s2yord.values()):
        for s in yord2ss[yord]:
            assert(yord == s2yord[s])

class TestGeomInvalidIntersection(unittest.TestCase):
    def test_general_intersection_bend(self):
        p1 = (0,0)
        q1 = (1,0)
        p2 = (0,0)
        q2 = (0,1)
        self.assertFalse(geom_has_general_intersection(p1, q1, p2, q2))

    def test_general_intersection_cross(self):
        p1 = (0,1)
        q1 = (2,1)
        p2 = (1,0)
        q2 = (1,3)
        self.assertTrue(geom_has_general_intersection(p1, q1, p2, q2))

    def test_general_intersection_straight(self):
        p1 = (0,0)
        q1 = (1,0)
        p2 = (1,0)
        q2 = (2,0)
        self.assertFalse(geom_has_general_intersection(p1, q1, p2, q2))

    def test_open_segment(self):
        a = (0, 0)
        b = (2, 2)
        c = (1, 1)
        self.assertTrue(geom_orientation(a, b, c) == 0)
        self.assertTrue(geom_box_ab_contain_c(a, b, c))
        self.assertTrue(geom_open_segment_ab_contain_c(a, b, c))

# ----------------------------------------

def make_dotfile_gridify():
    g = read_graph( nodes_file=STATES_INFILE_LATLONG
                    , edges_file=STATES_INFILE_EDGES
                    , order='yx'
                    )
    # Prevent the ordinalified ME-NH connection from crossing the MA-VT one.
    customize_graph(g, STATES_REPOSITIONS + STATES_REPOSITIONS_GRIDIFY)
    assert_graph_sanity(g)
    write_dotfile_gridify( g
                         , dotfile=STATES_DOTFILE_GRIDIFY
                         , verbose=True
                         )

def make_dotfile_planned(coords_file=None):
    if coords_file:
        print(f'Using custom coords file: <<{coords_file}>>')
    else:
        coords_file = STATES_INFILE_COORDS
        print(f'Using default coords file: {coords_file}')
    g = read_graph( nodes_file=coords_file
                  , edges_file=STATES_INFILE_EDGES
                  )
    customize_graph(g, STATES_REPOSITIONS)
    assert_graph_sanity(g)
    write_dotfile_planned( g
                         , dotfile=STATES_DOTFILE_PLANNED
                         )

def make_dotfile_springs():
    g = read_graph( nodes_file=STATES_INFILE_LATLONG
                  , edges_file=STATES_INFILE_EDGES
                  , order='yx'
                  )
    customize_graph(g, STATES_REPOSITIONS)
    assert_graph_sanity(g)
    write_dotfile_springs( g
                         , dotfile=STATES_DOTFILE_SPRINGS
                         )


if __name__ == '__main__':
    # unittest.main(argv=[''])
    prog='states_viz'
    parser = argparse.ArgumentParser( prog=prog
                                    , description='Visualize layout of US states'
                                    )
    parser.add_argument('-g', '--gridify',  action='store_true')
    parser.add_argument('-p', '--planned',  action='store_true')
    parser.add_argument('-s', '--springs',  action='store_true')
    parser.add_argument('-f', '--file',     action='store')
    args = parser.parse_args(sys.argv[1:])

    if (not args.gridify) and (not args.planned) and (not args.springs):
        print(f'{prog}: At least one output flag required.')
        parser.print_help(sys.stderr)
        sys.exit(1)
    if args.gridify:
        make_dotfile_gridify()
    if args.planned:
        make_dotfile_planned(coords_file=args.file)
    if args.springs:
        make_dotfile_springs()
