# CMC: Constrained Mixed-Radix Coding for Chinese Text DNA Storage

CMC is a reproducible code and data package for encoding Chinese text into DNA sequences under synthesis- and sequencing-oriented constraints. The pipeline combines Wubi-based text preprocessing, Brotli compression, constrained mixed-radix DNA mapping, dynamic GC/AT base reordering, and online homopolymer suppression.

This repository is intended for academic research on source-aware DNA storage coding. It provides the encoder/decoder implementation, mapping tables, metric summaries, and publication figure generation scripts.

# Link to the paper
Paper describing this work has been received in Frontiers of Computer Science（FCS） special column “Code & Data”.
Cited as: Liuming SHI, Ziqi LI, YanfenZHENG, XueLI, Ziniu Mu, Shihua ZHOU, Ben CAO, Bin WANG, EdmundMing KitLAI. 
CMC: Constrained Mixed-radix Coding Tools for DNA Storage of Chinese Text. 
Front. Comput. Sci., 2026, DOI: 10.1007/s11704-026-60971-6

## Benchmark Texts

The original benchmark documents are not shipped with this repository due to copyright and redistribution considerations. To reproduce the end-to-end document tests, prepare your own `.docx` files and place them under a local `data_raw/` directory, or update the paths in:

- `run_docx_test.m`
- `run_docx_test_trad_direct.m`
- `run_e2e_verification.py`

The reported experiments used two Chinese prose benchmarks: Zhu Ziqing's *Bei Ying* and a Traditional Chinese version of *Chu Shi Biao*. The repository keeps the code, mapping tables, metrics, and figures, but not the raw `.docx` input files.

## Example Inputs

To address different user needs, the repository provides a lightweight `examples/` directory with **20 small test cases grouped by category**. These files are not used as the formal benchmark in the manuscript; they are small public test cases for choosing an input style, checking the local environment, and verifying text-handling behavior before running the full workflow. Any character outside the Wubi-86 dictionary (punctuation, Latin letters, digits, and line breaks) is recorded verbatim and restored on decode, so every example round-trips losslessly.

| File | Use case | What it checks |
|:--|:--|:--|
| `examples/smoke_single_char.txt` | Smoke test (minimal) | Single-character input for the smallest possible round-trip check |
| `examples/quick_test.txt` | Smoke test | Minimal Simplified Chinese input for confirming that the pipeline starts correctly |
| `examples/simplified_modern.txt` | Simplified Chinese prose | Ordinary modern Chinese text storage |
| `examples/simplified_news.txt` | Simplified Chinese (news) | Journalistic / report-style sentence |
| `examples/simplified_technical.txt` | Simplified Chinese (technical) | Scientific-abstract style text with domain terminology |
| `examples/simplified_dialogue.txt` | Simplified Chinese (dialogue) | Conversational text with full-width quotation marks |
| `examples/simplified_poetry_classical.txt` | Classical poetry | Public-domain Tang poem with short lines |
| `examples/simplified_numbers_data.txt` | Numeric / data text | Number- and percentage-heavy content and numeric restoration |
| `examples/medium_simplified_paragraph.txt` | Medium-length input | Moderate compression and sequence-quality behavior |
| `examples/long_simplified_paragraph.txt` | Longer input | Compression behavior and sequence-quality checks on a longer paragraph |
| `examples/repeated_chars.txt` | Compression stress | Highly repetitive text and its compression behavior |
| `examples/multiparagraph_mixed.txt` | Multi-line input | Newline and paragraph restoration across lines |
| `examples/punctuation_mixed.txt` | Punctuation-rich input | Punctuation maps, brackets, quotation marks, and digits |
| `examples/punctuation_fullwidth.txt` | Full-width punctuation | Restoration coverage for the full-width Chinese punctuation set |
| `examples/mixed_alnum_symbols.txt` | Mixed content | Chinese text with English words, digits, percent signs, and symbols |
| `examples/mixed_english_terms.txt` | Mixed content (terms) | Chinese interleaved with English technical terms and base letters |
| `examples/mixed_dates_units.txt` | Mixed content (data) | Dates, times, room numbers, and currency units |
| `examples/traditional_classical.txt` | Traditional/classical-style input | Direct Traditional Chinese handling |
| `examples/traditional_variants.txt` | Traditional variants | Common Traditional characters and script preservation |
| `examples/traditional_culture.txt` | Traditional (modern) | Modern Traditional sentence for the direct Traditional encoding path |

The machine-readable `examples/manifest.csv` summarizes the same test cases. A MATLAB runner is also provided:

```matlab
run('examples/run_examples.m')
```

## Key Results

| Metric | CMC result | Notes |
|:--|:--|:--|
| Storage density | 1.94-2.17 bits/nt | Range across the tested direct-Traditional, conversion-assisted, and Simplified-input settings |
| Global GC content | 50.0% | Dynamic GC/AT reordering keeps the final composition balanced |
| Maximum homopolymer | 2 nt | Enforced during nucleotide generation |
| Undesired motif content | 0.000% | The predefined forbidden 3-mer set is avoided by construction |
| Minimum Hamming distance | 76-81 | Measured on 120 nt fragments for the CMC settings |
| Mean MFE | -27.8 to -31.1 kcal/mol | Moderate fragment-level MFE, avoiding over-stable secondary structures |
| Round-trip fidelity | 100% lossless | Verified for the included benchmark protocol |

The comparison tables in `metrics/` provide the full numerical values for CMC and the baseline methods. MFE values are reported for representative 120 nt fragments; large negative values in baseline rows should be interpreted under the same fragment-level calculation protocol rather than as per-nucleotide normalized quantities.

## Repository Contents

```text
CMC/
|-- Core MATLAB encoder/decoder
|   |-- cmc_brotli_encoder.m
|   |-- cmc_brotli_decoder.m
|   |-- wubi_encoder.m
|   |-- wubi_decoder.m
|   |-- get_allowed_bases.m
|   `-- reorder_bases_gc.m
|
|-- Mapping tables and utilities
|   |-- character.csv
|   |-- trad_mapping.csv
|   |-- cultural_extension.csv
|   |-- generate_trad_mapping.py
|   |-- tool_search_wubi.m
|   |-- export_dna_for_plots.m
|   `-- phrase_Wubi86.m
|
|-- Experiment and demo scripts
|   |-- brotli_compress.py
|   |-- run_cmc_demo.m
|   |-- run_brotli_demo.m
|   |-- run_cultural_demo.m
|   |-- run_docx_test.m
|   |-- run_docx_test_trad_direct.m
|   |-- run_e2e_verification.py
|   `-- run_full_experiment.py
|
|-- examples/
|   |-- README.md
|   |-- manifest.csv
|   |-- run_examples.m
|   `-- lightweight Chinese text test cases
|
|-- metrics/
|   `-- CSV files for density, GC, motif, homopolymer, Hamming distance, MFE, and summary tables
|
`-- plots/
    |-- generate_publication_figures.py
    |-- generate_comparison_charts.py
    |-- generate_supplementary_charts.py
    |-- generate_cmc_vs_drrc_chart.py
    `-- fig1-fig6 PNG/SVG outputs
```

The following local folders are intentionally not part of the public code package: `data_raw/`, `legacy/`, `analysis_results/`, `logs/`, `data_normalized/`, `cmrc_out/`, `external/`, `baseline_out/`, manuscript folders, submission packages, and Word document drafts.

## Installation

### Python

Use Python 3.8 or newer.

```bash
pip install -r requirements.txt
```

The Python-side dependencies are used for Brotli compression, document parsing, Traditional/Simplified conversion, metric aggregation, and figure generation.

### MATLAB

MATLAB R2019a or newer is recommended. The MATLAB scripts call Python packages, so configure MATLAB's Python bridge before running the demos:

```matlab
pyenv('Version', 'C:\path\to\python.exe')
```

On older MATLAB releases, use `pyversion` instead of `pyenv`.

## Quick Start

Run a short in-memory demo:

```matlab
run_cmc_demo
```

Run the Brotli-based long-text demo:

```matlab
run_brotli_demo
```

Run Traditional Chinese conversion-assisted encoding:

```matlab
run_cultural_demo
```

Run the lightweight example inputs:

```matlab
run('examples/run_examples.m')
```

For document-level tests, first place your own `.docx` files under `data_raw/` or edit the paths in the scripts, then run:

```matlab
run_docx_test
run_docx_test_trad_direct
```

Generate verification metrics and comparison outputs:

```bash
python run_e2e_verification.py
python run_full_experiment.py
python plots/generate_publication_figures.py
```

## Encoding Pipeline

```text
Chinese text
  -> Wubi preprocessing
  -> punctuation and collision-map recording
  -> Brotli compression
  -> constrained mixed-radix DNA mapping
  -> dynamic GC/AT base reordering
  -> online homopolymer and motif control
  -> DNA sequence
```

The decoding path reverses the same steps and restores punctuation and collision-disambiguated characters through the recorded maps.

## Constraint System

At each DNA position, `get_allowed_bases.m` removes bases that would violate local constraints. The remaining base set is then reordered by `reorder_bases_gc.m` to balance GC and AT counts dynamically.

| Constraint | Implementation |
|:--|:--|
| Motif avoidance | Rejects bases that would create predefined undesired 3-mers |
| Homopolymer control | Prevents runs longer than 2 nt |
| GC balancing | Reorders candidate bases according to current GC/AT counts |
| Mixed-radix capacity | Uses 2-, 3-, and 4-base candidate sets to encode bits efficiently |

## Metrics and Figures

The `metrics/` directory contains machine-readable CSV summaries, including:

- `main_result_table.csv`
- `cmc_e2e_verification.csv`
- `storage_density_summary.csv`
- `gc_summary.csv`
- `gc_detail_table.csv`
- `homopolymer_summary.csv`
- `motif_summary.csv`
- `hamming_summary.csv`
- `mfe_summary.csv`
- `sequence_quality_table.csv`
- `input_manifest.csv`
- `cmc_vs_drrc_summary.csv`

The `plots/` directory contains the final `fig1`-`fig6` PNG/SVG outputs and the Python scripts used to regenerate them.

## Version Notes

| Version | Main strategy | Density range | Key addition |
|:--|:--|:--|:--|
| v1.0 | Dual-table mapping | about 1.6 bits/nt | Initial constrained mapping |
| v2.0 | Adaptive mixed radix | about 1.7 bits/nt | Constraint-aware encoding |
| v3.0 | Brotli + mixed radix | about 1.94 bits/nt | Compression pipeline |
| v4.0 | Brotli + ternary grouping + GC balancing | 1.94-2.17 bits/nt | Direct Traditional support and improved constraint control |

## Citation

If you use this repository, please cite the accompanying manuscript:

> CMC: Constrained Mixed-Radix Coding for High-Density Chinese Text DNA Storage. Manuscript in preparation.

## License

This project is provided for academic research use. A formal license file has not yet been added; please contact the authors before redistribution or commercial use.
