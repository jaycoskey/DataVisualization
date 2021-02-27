#!/usr/bin/env python3
# Copyright Jay M. Coskey, 2021
# Note: This is a proof of concept. Efficiency and conciseness are not priorites.

import argparse
from collections import defaultdict, namedtuple
from copy import deepcopy
import networkx as nx
import sys
import unittest

# ----------------------------------------

GRAY_FILL = 'style=filled fillcolor=lightgray'
MAX_ITER_COUNT = 1000
MAX_PGAP = 10  # In gridify(), merge columns/rows only if within X% of total span.

DEFAULT_COORDS_STATES_V1 = 'states_v1.coords'
DEFAULT_COORDS_STATES_V2 = 'states_v2.coords'
DEFAULT_COORDS_STATES_V3 = 'states_v3.coords'

DEFAULT_DOTFILE_NATIONS_GRIDIFY = 'nations_gridify.dot'
DEFAULT_DOTFILE_NATIONS_PLANNED = 'nations_planned.dot'

DEFAULT_DOTFILE_STATES_GRIDIFY  = 'states_gridify.dot'
DEFAULT_DOTFILE_STATES_PLANNED  = 'states_planned.dot'
DEFAULT_DOTFILE_STATES_SPRINGS  = 'states_springs.dot'

DEFAULT_PLANNED_SCALE_FACTOR = 50
DEFAULT_PLANNED_FONTSIZE     = 10
DEFAULT_PLANNED_HEIGHT       = 0.45
DEFAULT_PLANNED_WIDTH        = 0.45

GRIDIFIED_COORDS_LOGFILE = 'gridified.coords'  # TODO: Control via args/config

INFILE_EDGES_NATIONS = 'nations.edges'
INFILE_EDGES_STATES  = 'states.edges'

INFILE_LATLONG_NATIONS = 'nations.latlong'
INFILE_LATLONG_STATES  = 'states.latlong'  # Note: International Date Line

STATES_DISCONNECTED = ['AK','HI']
STATES_FOUR_CORNERS_EDGES = [('AZ','CO'), ('NM','UT')]
STATES_GRIDIFY_REPOSITIONS = [('NH', ('MA', 'MN')), ('SC', ('NC', 'AZ'))]
STATES_SPRINGS_EDGES = [ ('AK', 'HI')  # Place near each other, but without visible edge'
                       , ('AK', 'CA')  # Place AK near CA
                       , ('HI', 'TX')  # Place HI near TX
                       , ('NJ', 'CT'), ('NJ', 'MD')  # Orient the northestern states properly
                       ]
STATES_SPRINGS_REPOSITIONS = [('AK', ('CA', 'TX')), ('HI', ('AZ', 'TX'))]

# ----------------------------------------

InvalidIntersection = namedtuple('InvalidIntersection', ['type_', 'names', 'coords'])

# ----------------------------------------

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

def has_items(iterable):
    return iterable and len(iterable) > 0

def maybe_dict_get(d, key, default):
    if d:
        return d.get(key, default)
    else:
        return default

def read_graph(nodes_file, edges_file, order='xy', has_groups=False, verbose=False):
    def read_nodes(nodes_file):
        with open(nodes_file, 'r') as f:
            for line in f.readlines():
                items = line.strip().split('#')[0].split()
                if len(items) == 0:
                    continue
                if has_groups:
                    assert(len(items) >= 4)
                    _, name, x, y = items[0:4]  # TODO: Support block display of groups (e.g., continents).
                else:
                    assert(len(items) == 3)
                    name, x, y = items
                if order == 'yx':  # TODO: Swap based on whether reading lat/long or coords file.
                    x, y = y, x
                yield (name, float(x), float(y))

    def read_edges(edges_file):
        with open(edges_file, 'r') as f:
            for line in f.readlines():
                vals = line.strip().split('#')[0].split(':')
                if len(vals) == 0:
                    continue
                assert(len(vals) == 2)
                name = vals[0].split(' ')[0]
                nbrs = vals[1].split() if len(vals) > 1 else []
                yield (name, nbrs)

    result = nx.Graph()

    verbose and print(f'Reading from nodes file: "{nodes_file}"')
    for name, x, y in read_nodes(nodes_file):
        result.add_node(name, x=x, y=y)

    verbose and print(f'Reading from edges file: "{edges_file}"')
    for node, nbrs in read_edges(edges_file):
        for nbr in nbrs:
            if nbr > node:
                result.add_edge(node, nbr)

    return result


# Used to place the disconnected states: AK & HI
def reposition_nodes(g, repositions):
    for reposition in repositions:
        s = reposition[0]
        x_src_name, y_src_name = reposition[1]
        x_src = int(g.nodes(data=True)[x_src_name]['x'])
        y_src = int(g.nodes(data=True)[y_src_name]['y'])

        nx.set_node_attributes(g, values={s:x_src}, name='x')
        nx.set_node_attributes(g, values={s:y_src}, name='y')


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
            for s in xord2ss[ordk]:
                assert(g.nodes[s]["x"] == ordk)
                nx.set_node_attributes(g, values={s:new_xord}, name='x')
        else:
            new_yord = ordk - 1
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
def write_dotfile_gridify( g
                         , dotfile
                         , custom_visible_edges=[]
                         , custom_visible_edge_attrs=''
                         , custom_invisible_edges=[]
                         , do_allow_invalid_intersections=False
                         , settings=None
                         , verbose=False):
    do_write_gridified = False

    def avg(items):
        return sum(items) / len(items)

    # TODO: Allow intersections that are present in initial ordinalized graph
    #   (The initial state of the ordinalized nations graph has 22 invalid intersections.)
    # Or these could somehow be resolved (e.g., using reposition_nodes) to make the graph planar.
    #
    def get_invalid_intersections(g):
        result = []

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
                if geom_has_general_intersection(p1, q1, p2, q2):
                    result.append(InvalidIntersection( 'general'
                                                     , f'({ep1},{eq1}),({ep2},{eq2})'
                                                     , f'({p1},{q1}),({p2},{q2})'
                                                     )
                                 )
                elif geom_has_open_segment_violation(p1, q1, p2, q2):
                    result.append(InvalidIntersection( 'openseg'
                                                     , f'({ep1},{eq1}),({ep2},{eq2})'
                                                     , f'({p1},{q1}),({p2},{q2})'
                                                     )
                                 )
        return result

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

        # TODO: Consider prioritizing squishing coords around large gaps in pre-ordinalified locations (e.g., Atlantic/Pacific oceans)
        pgaps = deepcopy(xord_pgaps)
        pgaps.update(yord_pgaps)
        pgaps = sorted(pgaps.items(), key=lambda item: item[1])

        assert_states_ords(s2xord, s2yord, xord2ss, yord2ss)
        return s2xord, s2yord, xord2ss, yord2ss, pgaps

    def gridify(g, max_iter_count=MAX_ITER_COUNT, do_allow_invalid_intersections=False, verbose=False):
        do_print_initial_invalid_intersections = False
        s2xy = {s:xy for s,xy in g.nodes(data=True)}

        s2x = lambda s: s2xy[s]['x']
        s2y = lambda s: s2xy[s]['y']

        gord = get_gord(g, s2xy)
        if do_print_initial_invalid_intersections:  # TODO: Control via args/config
            initial_invalid_intersections = get_invalid_intersections(gord)
            print('Initial invalid intersections:')
            for invalid_intersection in initial_invalid_intersections:
                print(f'\t{invalid_intersection}')

        iter_count = 0
        while True:
            if iter_count >= max_iter_count:
                break
            iter_count += 1

            s2xord, s2yord, xord2ss, yord2ss, pgaps = get_ord_info(gord, s2x, s2y)
            pgap_candidates = list(filter(lambda pgap_cand: pgap_cand[1] <= MAX_PGAP, pgaps))
            num_cands = len(pgap_candidates)
            verbose and print('='*20, f'INFO: Starting iteration {iter_count} with {num_cands} candidates.')

            # From all the candidates, unify the closest viable pair of (ordinal) coords.
            for cand_count in range(len(pgap_candidates)):
                print(' '*4, '-'*10, f'INFO: (iter,cand)=({iter_count},{cand_count}/{num_cands})')
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
                is_cand_valid_edges = (  do_allow_invalid_intersections
                                      or is_valid_intersections(gcand)
                                      )
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

    gridified = gridify( g
                       , max_iter_count=MAX_ITER_COUNT
                       , do_allow_invalid_intersections=do_allow_invalid_intersections
                       , verbose=verbose)

    if do_write_gridified:
        with open(GRIDIFIED_COORDS_LOGFILE, 'w') as f:
            writeln(f, 0, f'# Written by states_viz when creating "gridify" visualization.')
            s2xyords = {s:xyords for s,xyords in gridified.nodes(data=True)}
            for s, xyords in s2xyords.items():
                x = xyords["x"]
                y = xyords["y"]
                writeln(f, 0, f'{s} {x:>2} {y:>2}')
    write_dotfile_planned( gridified
                         , dotfile=dotfile
                         , custom_visible_edges=custom_visible_edges
                         , custom_visible_edge_attrs=custom_visible_edge_attrs
                         , custom_invisible_edges=custom_invisible_edges
                         , settings=settings
                         , verbose=verbose
                         )


# Write to named dotfile using coordinates in graph argument
#
def write_dotfile_planned( g
                         , dotfile
                         , custom_visible_edges=[]
                         , custom_visible_edge_attrs=''
                         , custom_invisible_edges=[]
                         , verbose=False
                         , settings=None  # dict w/ scale_factor, fontsize, height, width
                         ):

    def write_custom_invisible_edges(f):
        if custom_invisible_edges:
            for (a,b) in custom_invisible_edges:
                attrs = f' [style=invis]'
                writeln(f, 1, f'{a} -- {b}{attrs}')

    def write_custom_visible_edges(f):
        if custom_visible_edges:
            for (a,b) in custom_visible_edges:
                attrs = f' [{custom_visible_edge_attrs}]' if custom_visible_edge_attrs else ''
                writeln(f, 1, f'{a} -- {b}{attrs}')

    def write_edges(f):
        for a, b in sorted(g.edges(data=False)):
            if b > a:
                writeln(f, 1, f'{a} -- {b}')
        writeln(f)
        if has_items(custom_visible_edges):
            writeln(f)
        write_custom_visible_edges(f)
        if has_items(custom_invisible_edges):
            writeln(f)
        write_custom_invisible_edges(f)

    def write_nodes(f, scale_factor=50):
        for name, attrs in sorted(g.nodes(data=True)):
            x = attrs['x']
            y = attrs['y']
            pos_spec = f'pos="{scale_factor * x:.1f},{scale_factor * y:.1f}"'
            color_spec = f' {GRAY_FILL}' if name in STATES_DISCONNECTED else ''
            writeln(f, 1, f'{name} [{pos_spec}{color_spec}]')

    with open(dotfile, 'w') as f:
        do_write_sig = False  # Enable "signing" resulting viz.

        writeln(f, 0, 'strict graph States {')
        fontsize = maybe_dict_get(settings, 'fontsize', DEFAULT_PLANNED_FONTSIZE)
        height   = maybe_dict_get(settings, 'height',   DEFAULT_PLANNED_HEIGHT)
        width    = maybe_dict_get(settings, 'width',    DEFAULT_PLANNED_WIDTH)
        writeln(f, 1, f'node [fixedsize=true fontsize={fontsize} height={height} width={width}]')
        writeln(f, 1, 'overlap=false')
        writeln(f, 1, 'penwidth=7')
        writeln(f)
        write_nodes(f, maybe_dict_get(settings, 'scale_factor', 50))

        # TODO: Control "signature" initials, placement, etc. with args/config.
        if do_write_sig:
            writeln(f)
            write_sig(f, 'ABC', pos='525.0,25.0', color='lightblue', fontsize=14, height=0.6, width=0.6)
        writeln(f)
        write_edges(f)
        writeln(f, 0, '}')


# Write to named dotfile using Graphviz layout algorithm
#
def write_dotfile_springs( g
                         , dotfile
                         , custom_visible_edges=[]
                         , custom_visible_edge_attrs=''
                         , custom_invisible_edges=[]
                         , verbose=False
                         ):
    s2attrs = {s:attrs for s,attrs in g.nodes(data=True)}

    def write_custom_invisible_edges(f):
        for (a,b) in custom_invisible_edges:
            attrs = f' [style=invis]'
            writeln(f, 1, f'{a} -- {b}{attrs}')

    def write_custom_visible_edges(f):
        for (a,b) in custom_visible_edges:
            attrs = f' [{custom_visible_edge_attrs}]' if custom_visible_edge_attrs else ''
            writeln(f, 1, f'{a} -- {b}{attrs}')

    def write_edges(f):
        for a, b in sorted(g.edges(data=False)):
            if b > a:  # Avoid duplication
                writeln(f, 1, f'{a} -- {b}')
        if has_items(custom_visible_edges):
            writeln(f)
        write_custom_visible_edges(f)
        if has_items(custom_invisible_edges):
            writeln(f)
        write_custom_invisible_edges(f)

    def write_nodes(f):
        for s,attrs in sorted(g.nodes(data=True)):
            x = attrs['x']
            y = attrs['y']
            pos_attrs   = f'pos="{x:.2f},{y:.2f}"'
            color_attrs = f' {GRAY_FILL}' if s in STATES_DISCONNECTED else ''
            writeln(f, 1, f'{s} [{pos_attrs}{color_attrs}]')

    with open(dotfile, 'w') as f:
        writeln(f, 0, 'strict graph States {')
        writeln(f, 1, 'rankdir=LR')
        writeln(f, 1, 'splines=line')
        writeln(f, 1, 'node [fixedsize=true fontsize=10 height=0.45 width=0.45]')
        writeln(f, 1, 'penwidth=5')
        writeln(f)
        write_nodes(f)
        writeln(f)
        write_edges(f)

        writeln(f, 0, '}')


def write_invisible_edges(f, edges):
    for edge in edges:
        edge_str = f'{edge[0]} -- {edge[1]} [style=invis]'
        writeln(f, 1, edge_str)

def write_sig(f, initials, pos, color, fontsize, height, width):
    attrs = f'style=filled fillcolor={color} fontsize={fontsize} height={height} width={width}'
    writeln(f, 1, f'{initials} [pos="{pos}" {attrs}]')

def writeln(f, indent_level=0, s=''):
    f.write(' ' * 4 * indent_level)
    f.write(s + '\n')

# ----------------------------------------

def assert_states_graph_sanity(g):
    assert(len(g.nodes(data=False)) == 50)
    nbr_count = 2.0 * len(g.edges(data=False))
    node_count = len(g.nodes(data=False))
    avg_nbr_count = nbr_count / node_count
    # print(f'avg_nbr_count={avg_nbr_count}')
    assert(4.15 < avg_nbr_count < 4.25)

def assert_states_ords(s2xord, s2yord, xord2ss, yord2ss):
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

def make_dotfile_nations_gridify(dotfile=DEFAULT_DOTFILE_NATIONS_GRIDIFY):
    g = read_graph( nodes_file=INFILE_LATLONG_NATIONS
                  , edges_file=INFILE_EDGES_NATIONS
                  , order='yx'
                  , has_groups=True
                  )
    write_dotfile_gridify( g
                         , dotfile=dotfile
                         , custom_visible_edges=None
                         , custom_visible_edge_attrs=None
                         , custom_invisible_edges=None
                         , do_allow_invalid_intersections=True
                         , settings={"scale_factor":40, "fontsize":10, "height":0.35, "width":0.35}
                         , verbose=True
                         )

def make_dotfile_nations_planned(latlong_file=INFILE_LATLONG_NATIONS, dotfile=DEFAULT_DOTFILE_NATIONS_PLANNED):
    g = read_graph( nodes_file=INFILE_LATLONG_NATIONS
                  , edges_file=INFILE_EDGES_NATIONS
                  , order='yx'
                  , has_groups=True
                  )
    write_dotfile_planned( g
                         , dotfile=dotfile
                         , settings={"scale_factor":40, "fontsize":60, "height":2, "width":2}
                         )

def make_dotfile_states_gridify(dotfile=DEFAULT_DOTFILE_STATES_GRIDIFY):
    g = read_graph( nodes_file=INFILE_LATLONG_STATES
                  , edges_file=INFILE_EDGES_STATES
                  , order='yx'
                  )
    # Prevent the ordinalified ME-NH connection from crossing the MA-VT one.
    reposition_nodes(g, STATES_SPRINGS_REPOSITIONS + STATES_GRIDIFY_REPOSITIONS)
    assert_states_graph_sanity(g)
    write_dotfile_gridify( g
                         , dotfile=dotfile
                         , custom_visible_edges=STATES_FOUR_CORNERS_EDGES
                         , custom_visible_edge_attrs='color=red style=dashed'
                         , custom_invisible_edges=None
                         , do_allow_invalid_intersections=False
                         , settings=None # {"scale_factor":40, "fontsize":10, "height":1.5, "width":1.5}
                         , verbose=True
                         )

def make_dotfile_states_planned( coords_file=DEFAULT_COORDS_STATES_V3
                               , dotfile=DEFAULT_DOTFILE_STATES_PLANNED
                               ):
    print(f'Using coords file:     "{coords_file}"')
    print(f'Will write to dotfile: "{dotfile}"')
    g = read_graph( nodes_file=coords_file if coords_file else DEFAULT_COORDS_STATES_V3
                  , edges_file=INFILE_EDGES_STATES
                  )
    reposition_nodes(g, STATES_SPRINGS_REPOSITIONS)
    assert_states_graph_sanity(g)
    write_dotfile_planned( g
                         , dotfile=dotfile if dotfile else DEFAULT_DOTFILE_STATES_PLANNED
                         , custom_visible_edges=STATES_FOUR_CORNERS_EDGES
                         , custom_visible_edge_attrs='color=red style=dashed'
                         , settings={"height":0.35, "width":0.35}
                         )

def make_dotfile_states_springs(dotfile=DEFAULT_DOTFILE_STATES_SPRINGS):
    g = read_graph( nodes_file=INFILE_LATLONG_STATES
                  , edges_file=INFILE_EDGES_STATES
                  , order='yx'
                  )
    reposition_nodes(g, STATES_SPRINGS_REPOSITIONS)
    assert_states_graph_sanity(g)
    write_dotfile_springs( g
                         , dotfile=dotfile
                         , custom_visible_edges=STATES_FOUR_CORNERS_EDGES
                         , custom_visible_edge_attrs='color=red style=dashed'
                         , custom_invisible_edges=STATES_SPRINGS_EDGES
                         )

# The options form an inelegant set: gridify is a verb, planned is an adj; springs is a noun. Oh well.
if __name__ == '__main__':
    prog='states_viz'
    description='Visualize diagrammatic layout of nations or US states using different layouts.'
    parser = argparse.ArgumentParser(prog=prog , description=description)
    parser.add_argument('-c', '--coords_file', action='store',      help='Specify input coordinates file (x,y or lat,long)')
    parser.add_argument('-d', '--dotfile',     action='store',      help='Specify output dotfile')

    parser.add_argument('-n', '--nations',     action='store_true', help='Visualizing nations rather than US states')

    parser.add_argument('-g', '--gridify',     action='store_true', help='Align coordinated to a dense grid')
    parser.add_argument('-p', '--planned',     action='store_true', help='Plot graph using predetermined coordinates')
    parser.add_argument('-s', '--springs',     action='store_true', help='Use spring-based layout (Graphviz standard)')

    parser.add_argument('-t', '--test',        action='store_true', help='Run unit tests')

    args = parser.parse_args(sys.argv[1:])

    if args.test:
        unittest.main(argv=[''])
    if (not args.gridify) and (not args.planned) and (not args.springs):
        print(f'{prog}: At least one output flag required.')
        parser.print_help(sys.stderr)
        sys.exit(1)
    if args.gridify:
        if args.nations:
            make_dotfile_nations_gridify(dotfile=args.dotfile)
        else:
            make_dotfile_states_gridify(dotfile=args.dotfile)
    if args.planned:
        if args.nations:
            make_dotfile_nations_planned(latlong_file=args.coords_file, dotfile=args.dotfile)
        else:
            make_dotfile_states_planned(coords_file=args.coords_file, dotfile=args.dotfile)
    if args.springs:
        make_dotfile_states_springs(dotfile=args.dotfile)
