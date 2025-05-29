[![CI Status](https://github.com/aulyxair/protfetch/actions/workflows/ci.yml/badge.svg)](https://github.com/aulyxair/protfetch/actions/workflows/ci.yml)
[![PyPI version](https://img.shields.io/pypi/v/protfetch.svg)](https://pypi.org/project/protfetch/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/protfetch.svg)](https://pypi.org/project/protfetch/)
[![License](https://img.shields.io/pypi/l/protfetch.svg)](https://github.com/yourusername/protfetch/blob/main/LICENSE)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# protfetch

**protfetch** is a command-line tool designed to fetch and process protein FASTA sequences and associated metadata from NCBI's Entrez databases. Given a list of gene symbols (or protein names and gene symbols), it retrieves relevant protein sequences, filters them based on various criteria (identity, similarity, fragments), and outputs curated FASTA files and a metadata CSV.

## Features

- **Flexible Input**: Accepts gene lists in two formats:
  - `Protein Name | GENE_SYMBOL` (e.g., "Mitofusin 1 | MFN1")
  - `GENE_SYMBOL` (e.g., "MFN1")

- **NCBI Entrez Integration**: Uses Biopython to query NCBI Entrez for gene UIDs, linked protein UIDs, and FASTA sequences.

- **Keyword Filtering**: Filters fetched protein sequences based on keywords derived from the input (protein name or gene symbol) found in FASTA headers. This step can be skipped.

- **Sequence Processing & Filtering**:
  - Parses FASTA headers to extract accession and other identifiers
  - Removes duplicate entries by accession (first seen is kept)
  - Filters out identical sequences (keeps one representative, typically by lexicographical accession)
  - Filters near-identical sequences using Levenshtein distance (configurable threshold)
  - Removes fragment sequences (shorter sequences that are substrings of longer ones)

- **Concurrent Fetching**: Utilizes multiple workers to fetch data for different genes concurrently, speeding up processing for large lists.

- **Organized Output**:
  - Generates combined FASTA files (short headers and full headers) and a combined metadata CSV file
  - Optionally saves processed files for each gene individually
  - Combined files are deduplicated to ensure unique protein entries

## Installation

### Prerequisites
- Python 3.8 or higher
- pip (Python package installer)

```bash
pip install protfetch
```
## Usage

```bash
protfetch <input_gene_list_file> [OPTIONS]
```

### Required Arguments:
- `input_gene_list_file`: Path to the input file containing the list of genes.
  - Each line should be either `Protein Name | GENE_SYMBOL` or just `GENE_SYMBOL`
  - Lines starting with `#` are treated as comments and ignored

### NCBI Configuration:
- `--entrez-email YOUR_EMAIL`: **Strongly recommended.** Your email address for NCBI Entrez. NCBI requires this for reliable access and to contact you if there are issues with your queries.
- `--entrez-api-key YOUR_API_KEY`: *(**Optional, but highly recommended**)* Your NCBI API key for higher request rates.

### Key Options:
- `-o, --output-dir DIR`: Directory to save output files (default: `protfetch_results`)
- `--max-dist INT`: Max Levenshtein distance for filtering near-identical sequences (default: 4; 0 to disable)
- `--max-workers INT`: Maximum number of concurrent workers for fetching data (default: 5)
- `--save-individual-files`: Save processed FASTA and CSV for each gene individually
- `--skip-keyword-filter`: Skip filtering FASTA sequences by keyword
- `--timeout SECONDS`: Timeout for NCBI requests (default: 60)
- `--retries INT`: Number of retries for NCBI requests (default: 3)
- `--debug`: Enable detailed debug logging
- `-v, --version`: Show program's version number and exit
- `-h, --help`: Show help message and exit

### Example:

Create a sample gene list file (e.g., `genes.txt`):
```
# Mitofusin 1 | MFN1
# Calreticulin | CALR
# CANX
# PDIA3
```

Run protfetch:
```bash
protfetch genes.txt --entrez-email your.name@example.com -o my_protein_data --save-individual-files
```

This command will:
1. Read `genes.txt`
2. Fetch data for MFN1, CALR, CANX, and PDIA3 from NCBI
3. Filter and process the sequences
4. Save combined results (e.g., `genes_combo_short.fasta`, `genes_combo_full.fasta`, `genes_combo_meta.csv`) in the `my_protein_data` directory
5. Save individual files for each gene in `my_protein_data/individual_gene_files/`

## Output Files

In the specified output directory:

- `{input_file_stem}_combo_short.fasta`: Combined FASTA file with short headers (e.g., `>ACCESSION`), deduplicated by accession
- `{input_file_stem}_combo_full.fasta`: Combined FASTA file with full original headers, deduplicated by full header
- `{input_file_stem}_combo_meta.csv`: Combined metadata CSV file (accession, gene_input, identifier_from_header), deduplicated by accession

If `--save-individual-files` is used, a subdirectory (default: `individual_gene_files`) will contain:
- `{GENE_SYMBOL}_filtered_short.fasta`
- `{GENE_SYMBOL}_filtered_full.fasta`
- `{GENE_SYMBOL}_filtered_meta.csv`

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
