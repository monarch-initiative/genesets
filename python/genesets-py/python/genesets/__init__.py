"""Python bindings for genesets-rs.

``mgsa`` runs Bayesian model-based gene set analysis (MGSA; Bauer, Gagneur &
Robinson 2010), following the collapsed model of Holmes & Mungall's wtfgenes.
See the genesets-rs Statistics docs for the model.
"""

from genesets._genesets import __version__, mgsa

__all__ = ["__version__", "mgsa"]
