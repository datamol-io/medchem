from typing import cast
from typing import Union
from typing import List
from typing import Any

import datamol as dm
import numpy as np
import networkx as nx

from rdkit.Chem import SaltRemover

DEFAULT_NODE_ATTR = ["atomic_num", "degree", "ring_atom", "is_aromatic"]
DEFAULT_EDGE_ATTR = ["bond_type"]


def automorphism(
    mol: Union[str, dm.Mol],
    standardize: bool = True,
    node_attrs: List[str] = DEFAULT_NODE_ATTR,
    edge_attrs: List[str] = DEFAULT_EDGE_ATTR,
):
    """Compute automorphism in a molecular graph

    Args:
        mol: input molecular graph
        standardize: whether to standardize the compound or not
        node_attrs: list of categorical atom attributes/properties to consider for node matching
        edge_attrs: list of categorical bond attributes/properties to consider for edge matching
    """

    # there is no need to kekulize the molecule, in order to preserve
    # aromaticity instead of switching to double bonds
    with dm.without_rdkit_log():
        mol = dm.to_mol(mol, kekulize=False)
        # important to remove hydrogen to reduce automorphism
        # similarly there is no point in considering stereochemistry
        mol = dm.remove_hs(mol)

        if standardize:
            mol = dm.standardize_mol(mol, disconnect_metals=True, uncharge=True)
            # we also remove salts
            remover = SaltRemover.SaltRemover()
            mol = remover.StripMol(mol, dontRemoveEverything=True)

    # the graph is undirected by default
    graph = dm.graph.to_graph(mol)
    graph_copy = graph.copy()

    node_match, edge_match = None, None
    if node_attrs:
        node_match = nx.algorithms.isomorphism.categorical_node_match(node_attrs, [0] * len(node_attrs))
    if edge_attrs:
        edge_match = nx.algorithms.isomorphism.categorical_edge_match(edge_attrs, [0] * len(node_attrs))

    graph_matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph, graph_copy, node_match=node_match, edge_match=edge_match
    )
    matches = list(graph_matcher.match())
    return dict(graph=graph, matches=matches, mol=mol)


def score_symmetry(
    mol: Union[dm.Mol, str], exclude_self_mapped_edged: bool = False, **automorphism_kwargs: Any
):
    """Provide a symmetry score for a given input molecule

    !!! note
        This is an heuristic and our definition of symmetry is pretty loose.
        We define symmetry according to any (set of) plans dividing the molecule into two very similar subgraph.
        We include both edge and vertex transitivity. For example the star-molecular graph
        (e.g neopentane) is symmetrical here, although it's not vertex-transitive.

    Args:
        mol: inputs molecules
        exclude_self_mapped_edged: Whether to exclude edges that matches to themselves in automorphism.
        automorphism_kwargs: keyword for determining automorphism
    """

    def _is_self_mapped(mapping):
        return all(node1 == node2 for node1, node2 in mapping.items())

    mol_automorphism = automorphism(mol, **automorphism_kwargs)

    # We check and filter for identity mapping x->x for all nodes
    graph = mol_automorphism["graph"]
    graphmol = mol_automorphism["mol"]
    matches = mol_automorphism["matches"]
    matches = [x for x in matches if not _is_self_mapped(x)]

    # Always make the type explicit
    graphmol = cast(dm.Mol, graphmol)
    graph = cast(nx.Graph, graph)
    matches = cast(List[dict], matches)

    # Step 1: We set the symmetry score to 0 when there is not automorphism
    if len(matches) == 0:
        return 0

    # Step 2: let's find the coverage of edge satisfying the symmetric graph notion
    # https://en.wikipedia.org/wiki/Symmetric_graph

    arc_transitivity = []
    for bond in graphmol.GetBonds():
        edge = tuple(sorted([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]))
        is_covered = False
        for mapping in matches:
            iso_edge = (mapping.get(edge[0]), mapping.get(edge[1]))
            if graph.has_edge(*iso_edge) and (iso_edge != edge):
                is_covered = True
                break
            elif set(iso_edge) == set(edge) and not exclude_self_mapped_edged:
                # we only accept those when at least some neighbors of the nodse in the edge map to each other
                neighbors1 = set(graph.neighbors(edge[0]))
                neighbors2 = set(graph.neighbors(edge[1]))
                if any((mapping.get(x) in neighbors1 - {x}) for x in neighbors1) and any(
                    (mapping.get(x) in neighbors2 - {x}) for x in neighbors2
                ):
                    is_covered = True
                    break
        arc_transitivity.append(is_covered)
    return np.mean(arc_transitivity)
