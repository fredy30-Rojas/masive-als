#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Ajusta MAX_MEMORY_CHARS a 8000 y max_tokens de la cadena de chat a 2048."""
import pathlib, re

base = pathlib.Path("/home/fredy/claude-bot")

# 1) web_server.py: memoria 20k -> 8k
p = base / "web_server.py"
t = p.read_text(encoding="utf-8")
t2 = t.replace("MAX_MEMORY_CHARS = 20000", "MAX_MEMORY_CHARS = 8000")
t2 = t2.replace('("Cerebras", cerebras_client, "gemma-4-31b", 1024)', '("Cerebras", cerebras_client, "gemma-4-31b", 2048)')
t2 = t2.replace('("Groq", groq_client, "openai/gpt-oss-120b", 1024)', '("Groq", groq_client, "openai/gpt-oss-120b", 2048)')
if t2 != t:
    p.write_text(t2, encoding="utf-8")
    print("web_server.py OK")
else:
    print("web_server.py sin cambios")

# 2) claude_bot.py y userbot.py: max_tok 1024 -> 2048 en la cadena
for f in ("claude_bot.py", "userbot.py"):
    p = base / f
    t = p.read_text(encoding="utf-8")
    t2 = t.replace('"openai/gpt-oss-120b", 1024)', '"openai/gpt-oss-120b", 2048)')
    t2 = t2.replace('"gpt-oss-120b", 1024)', '"gpt-oss-120b", 2048)')
    if t2 != t:
        p.write_text(t2, encoding="utf-8")
        print(f, "OK")
    else:
        print(f, "sin cambios")

# 3) wa_bot_limpio.js: max_tokens 1024 -> 2048 en la cadena de chat
p = base / "whatsapp-bot" / "wa_bot_limpio.js"
t = p.read_text(encoding="utf-8")
t2 = t.replace("max_tokens: 1024, temperature: 0.7", "max_tokens: 2048, temperature: 0.7")
if t2 != t:
    p.write_text(t2, encoding="utf-8")
    print("wa_bot_limpio.js OK")
else:
    print("wa_bot_limpio.js sin cambios")

# verificacion
for f, pat in [("web_server.py", r"MAX_MEMORY_CHARS = 8000"),
               ("web_server.py", r"gpt-oss-120b\", 2048"),
               ("claude_bot.py", r"gpt-oss-120b\", 2048"),
               ("userbot.py", r"gpt-oss-120b\", 2048"),
               ("whatsapp-bot/wa_bot_limpio.js", r"max_tokens: 2048")]:
    txt = (base / f).read_text(encoding="utf-8")
    print("VERIF", f, bool(re.search(pat, txt)))
