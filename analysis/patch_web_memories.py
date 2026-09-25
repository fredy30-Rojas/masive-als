#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Parche: limita las memorias de la web a 20k chars por peticion."""
import pathlib

p = pathlib.Path("/home/fredy/claude-bot/web_server.py")
lines = p.read_text(encoding="utf-8").splitlines()

idx = next(k for k, l in enumerate(lines) if l.startswith("def load_memories():"))
new_fn = [
    "MAX_MEMORY_CHARS = 20000  # presupuesto de memorias por peticion (evita 400 por prompt gigante)",
    "",
    "def load_memories():",
    "    if not MEMORY_DIR.exists():",
    '        return ""',
    "    parts = []",
    "    total = 0",
    '    for f in sorted(MEMORY_DIR.glob("*.md")):',
    '        txt = f.read_text(encoding="utf-8", errors="replace")',
    "        if total + len(txt) > MAX_MEMORY_CHARS:",
    "            restante = MAX_MEMORY_CHARS - total",
    "            if restante > 2000:",
    '                parts.append(f"{f.stem}: {txt[:restante]}")',
    "            break",
    '        parts.append(f"{f.stem}: {txt}")',
    "        total += len(txt)",
    '    return "\\n\\n".join(parts)',
]
j = idx + 1
while j < len(lines) and not lines[j].startswith("def "):
    j += 1
lines[idx:j] = new_fn
p.write_text("\n".join(lines) + "\n", encoding="utf-8")
print("PARCHE_OK lineas:", idx + 1, "-", idx + len(new_fn))
