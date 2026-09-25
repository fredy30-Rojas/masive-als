#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Ajuste final: memorias 4k chars y max_tokens 1500 (límite 8k TPM de Groq)."""
import pathlib, re

base = pathlib.Path("/home/fredy/claude-bot")

p = base / "web_server.py"
t = p.read_text(encoding="utf-8")
t2 = t.replace("MAX_MEMORY_CHARS = 8000", "MAX_MEMORY_CHARS = 4000")
t2 = t2.replace('("Cerebras", cerebras_client, "gemma-4-31b", 2048)', '("Cerebras", cerebras_client, "gemma-4-31b", 1500)')
t2 = t2.replace('("Groq", groq_client, "openai/gpt-oss-120b", 2048)', '("Groq", groq_client, "openai/gpt-oss-120b", 1500)')
p.write_text(t2, encoding="utf-8")
print("web_server.py OK")

for f in ("claude_bot.py", "userbot.py"):
    p = base / f
    t = p.read_text(encoding="utf-8")
    t2 = t.replace('"openai/gpt-oss-120b", 2048)', '"openai/gpt-oss-120b", 1500)')
    t2 = t2.replace('"gpt-oss-120b", 2048)', '"gpt-oss-120b", 1500)')
    p.write_text(t2, encoding="utf-8")
    print(f, "OK")

p = base / "whatsapp-bot" / "wa_bot_limpio.js"
t = p.read_text(encoding="utf-8")
t2 = t.replace("max_tokens: 2048, temperature: 0.7", "max_tokens: 1500, temperature: 0.7")
p.write_text(t2, encoding="utf-8")
print("wa_bot_limpio.js OK")

for f, pat in [("web_server.py", r"MAX_MEMORY_CHARS = 4000"),
               ("web_server.py", r"gpt-oss-120b\", 1500"),
               ("claude_bot.py", r"gpt-oss-120b\", 1500"),
               ("userbot.py", r"gpt-oss-120b\", 1500"),
               ("whatsapp-bot/wa_bot_limpio.js", r"max_tokens: 1500")]:
    txt = (base / f).read_text(encoding="utf-8")
    print("VERIF", f, bool(re.search(pat, txt)))
