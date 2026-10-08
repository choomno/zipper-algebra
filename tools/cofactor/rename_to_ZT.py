#!/usr/bin/env python3
"""Compatibility shim.

Historically this script translated R_22/R_25/R_26/R_28 to Z_2/T_2/T_2_j2/T_3
in the cofactor .txt output. After the canonical-labeling refactor (see
relations_naming.py), the generator scripts emit Z/T labels directly, so this
script reduces to a copy of cofactor-X3-n4-Qqt-three-rels.txt into
cofactor-X3-n4-Qqt-ZT.txt for backwards compatibility with downstream
consumers that hard-code the ZT filename.
"""
from __future__ import annotations
import os
import shutil

ROOT = "/home/choom/braid-research/zipper-algebra"
SRC = os.path.join(ROOT, "logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-three-rels.txt")
DST = os.path.join(ROOT, "logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-ZT.txt")

shutil.copyfile(SRC, DST)
print(f"Copied {SRC} -> {DST}")
print("(Three-rels file already uses Z_2/T_2/T_3 labels; copy is a no-op rename.)")
