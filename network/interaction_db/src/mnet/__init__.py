"""mnet - metabolite extension of the STRING human protein network.

Builds ONE undirected edge-list (edges.csv) plus a nodes.csv that adds
metabolites to the existing STRING protein-protein network:

  * protein-protein  : the STRING build output, unchanged (edge_type = ppi)
  * protein-metabolite : an enzyme catalyses a Rhea reaction using the
                         metabolite (edge_type = catalysis)
  * lipid -> lipid-class : lipid_is_a hierarchy edges

This package sits next to ``string_network`` and reuses its helpers
(config, download, load). It never changes ``string_network`` behaviour.

This is the Phase 0 scaffold: modules exist and the CLI is wired, but no data
processing is implemented yet.
"""

__version__ = "0.0.1"
