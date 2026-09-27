"""Reproducible pipeline to rebuild the human STRING v12.0 PPI network.

Downloads public STRING and UniProt data, maps ENSP -> UniProt, deduplicates
into an undirected network, filters by combined score and writes a parquet file
matching the internal reference format.
"""

__version__ = "0.1.0"
