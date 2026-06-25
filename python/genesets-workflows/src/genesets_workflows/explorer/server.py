from __future__ import annotations

import argparse
import sys
import threading
import webbrowser
from pathlib import Path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="genesets-workflows explore",
        description="Launch a local browser for genesets-rs workflow report bundles.",
    )
    parser.add_argument(
        "bundles",
        nargs="*",
        type=Path,
        help="Report directories or summary.yaml files. Defaults to notebooks/generated/*/summary.yaml.",
    )
    parser.add_argument("--host", default="127.0.0.1", help="Host interface for the local server.")
    parser.add_argument("--port", type=int, default=8765, help="Port for the local server.")
    parser.add_argument("--open", action="store_true", help="Open the explorer in the default browser.")
    parser.add_argument("--reload", action="store_true", help="Enable uvicorn reload for UI/backend development.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        import uvicorn  # type: ignore

        from genesets_workflows.explorer.app import create_app
    except ImportError as exc:
        raise SystemExit(
            "The explorer requires optional web dependencies. Run with:\n"
            "  uv run --project python/genesets-workflows --extra explorer "
            "genesets-workflows explore <bundle>\n"
            "or install:\n"
            "  python3 -m pip install -e 'python/genesets-workflows[explorer]'"
        ) from exc

    app = create_app(args.bundles)
    url = f"http://{args.host}:{args.port}"
    print(f"genesets-rs explorer: {url}", file=sys.stderr)
    if args.open:
        threading.Timer(0.8, lambda: webbrowser.open(url)).start()
    uvicorn.run(app, host=args.host, port=args.port, reload=args.reload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

