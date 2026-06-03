from __future__ import annotations

import json
from typing import TextIO


def yaml_scalar(value: object) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    text = str(value)
    if text == "" or any(ch in text for ch in ":#[]{}&,*!|>'\"%@`") or text.strip() != text:
        return json.dumps(text)
    return text


def write_yaml_value(handle: TextIO, value: object, indent: int = 0) -> None:
    pad = " " * indent
    if isinstance(value, dict):
        for key, item in value.items():
            if isinstance(item, (dict, list)):
                handle.write(f"{pad}{key}:\n")
                write_yaml_value(handle, item, indent + 2)
            else:
                handle.write(f"{pad}{key}: {yaml_scalar(item)}\n")
    elif isinstance(value, list):
        for item in value:
            if isinstance(item, (dict, list)):
                handle.write(f"{pad}-\n")
                write_yaml_value(handle, item, indent + 2)
            else:
                handle.write(f"{pad}- {yaml_scalar(item)}\n")
    else:
        handle.write(f"{pad}{yaml_scalar(value)}\n")
