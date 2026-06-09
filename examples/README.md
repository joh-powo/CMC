# CMC Example Inputs

This directory provides a curated set of small Chinese-text inputs for quick CMC
tests before users run the full benchmark suite. The examples are intentionally
lightweight: they let users **choose an input type that matches their own needs**,
check text-handling behaviour, and confirm that the encode → DNA → decode workflow
is configured correctly.

## Why These Examples

In response to a reviewer suggestion to *provide more test cases so that potential
users can choose examples suited to their own needs*, this set has been expanded to
**20 inputs spanning the common Chinese-text DNA-storage scenarios**: smoke tests,
several Simplified Chinese genres, classical poetry, numeric/data text, punctuation
restoration (half- and full-width), mixed Chinese–English–digit content, multi-line
input, a compression-stress case, and Traditional/classical text for the direct
Traditional encoding path.

> **Note on robustness.** Any character that is not found in the Wubi-86 dictionary
> (punctuation, Latin letters, digits, line breaks, rare glyphs) is recorded
> verbatim in the punctuation/collision map and restored on decode. Every example
> therefore round-trips losslessly; characters outside the Wubi table are simply
> stored literally rather than densely encoded.

## Files by Category

### Smoke tests
| File | Length | Purpose |
| --- | --- | --- |
| `smoke_single_char.txt` | 1 char | Smallest possible round-trip check. |
| `quick_test.txt` | short | Minimal sentence to confirm the pipeline starts. |

### Simplified Chinese (by genre)
| File | Length | Purpose |
| --- | --- | --- |
| `simplified_modern.txt` | short | Modern prose, ordinary text storage. |
| `simplified_news.txt` | short | Journalistic / report style. |
| `simplified_technical.txt` | medium | Scientific-abstract style with domain terms. |
| `simplified_dialogue.txt` | short | Conversational text with quotation marks. |
| `simplified_poetry_classical.txt` | short | Classical Tang poem (public domain, Li Bai). |
| `simplified_numbers_data.txt` | short | Number- and percentage-heavy text. |

### Length scaling & compression
| File | Length | Purpose |
| --- | --- | --- |
| `medium_simplified_paragraph.txt` | medium | Moderate-length compression test. |
| `long_simplified_paragraph.txt` | long | Larger input for compression & constraint metrics. |
| `repeated_chars.txt` | short | Highly repetitive text (compression behaviour). |
| `multiparagraph_mixed.txt` | short | Multi-line input; tests newline/paragraph restoration. |

### Punctuation restoration
| File | Length | Purpose |
| --- | --- | --- |
| `punctuation_mixed.txt` | short | Brackets, quotes, and digits. |
| `punctuation_fullwidth.txt` | medium | Full-width Chinese punctuation coverage. |

### Mixed content
| File | Length | Purpose |
| --- | --- | --- |
| `mixed_alnum_symbols.txt` | short | Chinese + English words + digits + symbols. |
| `mixed_english_terms.txt` | short | Chinese interleaved with English technical terms. |
| `mixed_dates_units.txt` | short | Dates, times, room numbers, currency units. |

### Traditional Chinese (direct encoding path)
| File | Length | Purpose |
| --- | --- | --- |
| `traditional_classical.txt` | short | Classical-style Traditional text. |
| `traditional_variants.txt` | short | Common Traditional characters / script preservation. |
| `traditional_culture.txt` | short | Modern Traditional sentence on cultural preservation. |

The same information is available in machine-readable form in
[`manifest.csv`](manifest.csv) (columns: `filename, category, script_form,
length_class, purpose, copyright_note`).

## Suggested Use

1. Start with `smoke_single_char.txt` or `quick_test.txt` to confirm your local
   environment (MATLAB + Python Brotli bridge) is configured.
2. Then pick an example that matches your intended application:
   - Traditional Chinese → `traditional_classical.txt`, `traditional_variants.txt`,
     `traditional_culture.txt`
   - Punctuation recovery → `punctuation_mixed.txt`, `punctuation_fullwidth.txt`
   - Mixed Chinese/English/digits → `mixed_*` files
   - Larger compression behaviour → `medium_*`, `long_*`, `repeated_chars.txt`

## Running the Examples

A MATLAB runner is provided. It auto-discovers every file listed in
`manifest.csv`, runs the full pipeline, and prints per-example metrics plus a
round-trip PASS/FAIL summary:

```matlab
% from the project root
run('examples/run_examples.m')
```

Expected per-example output includes input length, Wubi length, DNA length,
GC content, maximum homopolymer run, information density (bits/nt), and the
round-trip result. To add your own case, drop a UTF-8 `.txt` file in this folder
and add one row to `manifest.csv`.

> The full manuscript benchmark remains based on the curated `data_raw/` and
> `data_normalized/` records. These examples are intended for usability testing
> and demonstration, not for replacing the reported benchmark results.
