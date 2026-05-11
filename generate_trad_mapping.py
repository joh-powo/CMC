# -*- coding: utf-8 -*-
"""
Generate trad_mapping.csv: traditional character -> wubi code extension table.
Two-pass approach:
  Pass 1: For each simplified char, use s2t to find traditional form
  Pass 2: Scan all CJK chars, use t2s to find simplified form in dict (catches reverse misses)
"""
import csv, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import opencc

converter_s2t = opencc.OpenCC('s2t')
converter_t2s = opencc.OpenCC('t2s')

# Load existing dictionary
dict_chars = set()
char_to_code = {}

with open('character.csv', 'r', encoding='utf-8') as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    for row in reader:
        code = row[0].strip()
        for col_idx in range(1, min(5, len(row))):
            ch = row[col_idx].strip()
            if ch and len(ch) == 1:
                dict_chars.add(ch)
                char_to_code[ch] = code

print(f"Loaded {len(dict_chars)} single characters from dictionary")

# Pass 1: s2t direction
trad_map = {}  # trad_char -> (code, simp_char)

for simp_ch, code in char_to_code.items():
    trad_ch = converter_s2t.convert(simp_ch)
    if (len(trad_ch) == 1 and 
        trad_ch != simp_ch and 
        trad_ch not in dict_chars and 
        trad_ch not in trad_map):
        trad_map[trad_ch] = (code, simp_ch)

print(f"Pass 1 (s2t): {len(trad_map)} mappings")

# Pass 2: scan all CJK codepoints, use t2s to find simplified in dict
added_pass2 = 0
for cp in range(0x4E00, 0xA000):
    trad_ch = chr(cp)
    if trad_ch in dict_chars or trad_ch in trad_map:
        continue  # already handled
    simp_ch = converter_t2s.convert(trad_ch)
    if (len(simp_ch) == 1 and 
        simp_ch != trad_ch and 
        simp_ch in char_to_code):
        code = char_to_code[simp_ch]
        trad_map[trad_ch] = (code, simp_ch)
        added_pass2 += 1

print(f"Pass 2 (t2s scan): +{added_pass2} mappings")
print(f"Total: {len(trad_map)} traditional character mappings")

# Sort by wubi code
trad_entries = [(code, trad_ch, simp_ch) for trad_ch, (code, simp_ch) in trad_map.items()]
trad_entries.sort(key=lambda x: x[0])

# Write CSV
with open('trad_mapping.csv', 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(['五笔码', '字1', '字2', '字3', '字4'])
    for code, trad_ch, simp_ch in trad_entries:
        writer.writerow([code, trad_ch, '', '', ''])

print(f"Written to trad_mapping.csv")

# Verification: check 出師表 coverage
from docx import Document
doc = Document('data_raw/繁体版出師表.docx')
text = ''.join(p.text for p in doc.paragraphs)

all_chars = dict_chars | set(trad_map.keys())
missing_cjk = []
missing_punct = []
for ch in set(text):
    if ord(ch) > 127 and ch not in all_chars:
        # Check if it's punctuation (handled by punctMap)
        simp = converter_t2s.convert(ch)
        if simp == ch and ch not in char_to_code:
            # Not a CJK char with mapping - likely punctuation
            if ord(ch) < 0x4E00 or ord(ch) > 0x9FFF:
                missing_punct.append(ch)
            else:
                missing_cjk.append(ch)

print(f"\n出師表 coverage check:")
print(f"  CJK chars missing: {len(missing_cjk)}")
for ch in missing_cjk:
    print(f"    {ch} (U+{ord(ch):04X})")
print(f"  Punctuation (handled by punctMap): {len(missing_punct)}")
for ch in missing_punct:
    print(f"    {ch} (U+{ord(ch):04X})")
