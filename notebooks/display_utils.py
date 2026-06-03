"""Display helpers for genesets-rs tutorial notebooks."""

from __future__ import annotations

import subprocess
import html
from pathlib import Path
from typing import Iterable

import pandas as pd
from IPython.display import HTML, display


def repo_root() -> Path:
    return Path(
        subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            text=True,
        ).strip()
    )


def metric_cards(items: Iterable[tuple[str, str]]) -> None:
    cards = "".join(
        f'<div class="notebook-metric"><span class="notebook-metric-value">{value}</span>'
        f'<span class="notebook-metric-label">{label}</span></div>'
        for label, value in items
    )
    display(HTML(f'<div class="notebook-metrics">{cards}</div>'))


def display_table(
    df: pd.DataFrame,
    *,
    columns: list[str] | None = None,
    caption: str | None = None,
    max_rows: int | None = 25,
    formatters: dict | None = None,
) -> None:
    view = df if columns is None else df.loc[:, columns]
    if max_rows is not None:
        view = view.head(max_rows)
    view = view.copy()
    formatted_columns = set(formatters or {})
    for column in view.columns:
        if column not in formatted_columns:
            view[column] = view[column].astype(object).where(pd.notna(view[column]), "")
    safe_formatters = {
        column: (lambda value, formatter=formatter: "" if pd.isna(value) else formatter(value))
        for column, formatter in (formatters or {}).items()
    }
    table = view.to_html(
        index=False,
        classes="notebook-table",
        border=0,
        escape=True,
        justify="left",
        formatters=safe_formatters,
        na_rep="",
    )
    if caption:
        table = f"<figure><figcaption>{html.escape(caption)}</figcaption>{table}</figure>"
    display(HTML(table))


def callout(message: str) -> None:
    display(HTML(f'<div class="notebook-callout">{message}</div>'))
