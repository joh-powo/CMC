# CMC Example Inputs

This directory provides small Chinese text inputs for quick CMC tests before users run the full benchmark suite. The examples are intentionally lightweight: they help users choose a suitable input type, check text handling behavior, and confirm that the encoding workflow is configured correctly.

## Why These Examples

The reviewer asked for more test cases so that users can choose examples matching their own needs. We therefore include examples that cover common Chinese text DNA-storage scenarios:

- short smoke test input;
- Simplified Chinese modern prose;
- Traditional/classical-style Chinese text;
- punctuation-rich Chinese text;
- mixed Chinese, English, digits, and symbols;
- longer Simplified paragraph for compression-oriented tests;
- common Traditional-character variants for script-preservation checks.

## Files

| File | Category | Purpose |
| --- | --- | --- |
| `quick_test.txt` | Smoke test | Minimal input for checking whether the pipeline starts correctly. |
| `simplified_modern.txt` | Simplified Chinese | Modern prose for ordinary Chinese text storage tests. |
| `traditional_classical.txt` | Traditional/classical | Classical-style text for direct Traditional processing. |
| `punctuation_mixed.txt` | Punctuation restoration | Tests punctuation maps, brackets, quotation marks, and digits. |
| `mixed_alnum_symbols.txt` | Mixed content | Tests Chinese text with English words, digits, percent signs, and symbols. |
| `long_simplified_paragraph.txt` | Longer input | Useful for observing compression and sequence-quality behavior. |
| `traditional_variants.txt` | Traditional variants | Checks common Traditional characters and script preservation. |
| `run_examples.m` | MATLAB runner | Runs all example inputs through the lightweight encode/decode check. |

The same information is provided in machine-readable form in `manifest.csv`.

## Suggested Use

Use `quick_test.txt` first to confirm the local environment. Then choose an example that matches the intended application: `traditional_classical.txt` or `traditional_variants.txt` for Traditional Chinese, `punctuation_mixed.txt` for punctuation recovery, `mixed_alnum_symbols.txt` for mixed-content input, and `long_simplified_paragraph.txt` for a slightly larger compression test.

To run all examples from MATLAB:

```matlab
run('examples/run_examples.m')
```

The full manuscript benchmark remains based on the curated `data_raw/` and `data_normalized/` records. These examples are intended for usability testing and demonstration, not for replacing the reported benchmark results.
