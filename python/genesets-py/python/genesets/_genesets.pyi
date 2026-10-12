from collections.abc import Iterable
from typing import Any

__version__: str

def mgsa(
    annotations: Iterable[tuple[str, str]],
    study: Iterable[str],
    background: Iterable[str] | None = None,
    closure: Iterable[tuple[str, str]] | None = None,
    *,
    samples_per_term: int = 100,
    burn_per_term: int = 10,
    term_prob: float = 0.5,
    term_count: float = 0.0,
    false_pos_prob: float = 0.5,
    false_pos_count: float = 0.0,
    false_neg_prob: float = 0.5,
    false_neg_count: float = 0.0,
    flip_rate: float = 1.0,
    step_rate: float = 1.0,
    jump_rate: float = 1.0,
    seed: int = 123456789,
) -> dict[str, Any]:
    """Run MGSA for one study set.

    ``annotations`` are (gene, term) pairs already propagated to ancestors.
    ``background`` defaults to every annotated gene. ``closure`` is optional
    (child, ancestor) pairs used for neighbor moves and to pick the most
    specific of several equivalent terms.

    Returns ``{"terms": [...], "stats": {...}, "study_size": int,
    "background_size": int}``. Each term has ``id``, ``posterior``,
    ``overlap``, ``term_size``, ``p_value`` (term-for-term hypergeometric),
    ``equivalent_ids``, and ``genes`` (study genes in the term); terms are
    sorted by posterior, highest first.
    """
