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

- **CD-HIT Integration**:
  - `--enable-cdhit (-ec)`: Enables adaptive redundancy reduction using CD-HIT. Uses automatic thresholds unless `--cdhit-fixed-threshold` is set.
  - `--cdhit-fixed-threshold <float>`: Specify a fixed CD-HIT identity threshold (e.g., 0.9 for 90%). Overrides auto-thresholding.
  - Automatic thresholds per gene/class based on sequence count:
    - \> 1000 seqs: ~70% identity
    - 500-1000 seqs: ~80% identity
    - 250-500 seqs: ~85% identity
    - 100-250 seqs: ~90% identity
    - < 100 seqs: ~95% identity
  - Requires cd-hit executable to be installed and in system PATH.
    
- **Solo Redundancy Filtering Mode**:
  - `--solo-redundancy-filtering (-srf)`: Operates on existing FASTA and CSV files to apply CD-HIT based redundancy filtering.
  - Requires `--input-fasta-srf`, `--input-csv-srf`, and CD-HIT to be enabled (via `--enable-cdhit` or `--cdhit-fixed-threshold`).
  - `--cdhit-group-workers INT`: Number of parallel workers for CD-HIT in SRF mode (default: 6).
  - Bypasses NCBI fetching.

- **Concurrent Fetching**: Utilizes multiple workers to fetch data for different genes concurrently, speeding up processing for large lists.

- **Organized Output**:
  - Generates combined FASTA files (short headers and full headers) and a combined metadata CSV file
  - Optionally saves processed files for each gene individually
  - Combined files are deduplicated to ensure unique protein entries

## Installation

### Prerequisites
- Python 3.8 or higher
- pip (Python package installer)
- **CD-HIT**: If using CD-HIT features (`--enable-cdhit`, `--cdhit-fixed-threshold`, or `--solo-redundancy-filtering`), the cd-hit executable must be installed and accessible in your system's PATH (See CD-HIT website for installation, or use `conda install bioconda::cd-hit`).
```bash
pip install protfetch
```
## Usage
### Standard Fetching Mode:

```bash
protfetch <input_gene_list_file> [OPTIONS]
```

**Example: Fetch, process, and apply CD-HIT auto-thresholding:**

```bash
protfetch genes.txt --entrez-email your@email.com -o results --enable-cdhit
```

**To rely more on CD-HIT for similarity, disable protfetch's Levenshtein filter:**

```bash
protfetch genes.txt --entrez-email your@email.com -o results --enable-cdhit --max-dist 0
```

**Use a fixed CD-HIT threshold of 98%:**

```bash
protfetch genes.txt --entrez-email your@email.com -o results --enable-cdhit --cdhit-fixed-threshold 0.98
```

### Solo Redundancy Filtering Mode:

```bash
protfetch --solo-redundancy-filtering --input-fasta-srf <path_to_fasta> --input-csv-srf <path_to_csv> --enable-cdhit -o <output_dir> [OPTIONS]
```

**Example: Apply CD-HIT auto-thresholding to existing files using 4 workers for CD-HIT groups:**

```bash
protfetch --solo-redundancy-filtering \
          --input-fasta-srf my_proteins.fasta \
          --input-csv-srf my_metadata.csv \
          --enable-cdhit \
          --cdhit-group-workers 4 \
          -o filtered_results
```

### Required Arguments:
- `input_gene_list_file`: Path to the input file containing the list of genes (required unless in -srf mode).
  - Each line should be either `Protein Name | GENE_SYMBOL` or just `GENE_SYMBOL`
  - Lines starting with `#` are treated as comments and ignored

### NCBI Configuration:
- `--entrez-email YOUR_EMAIL`: **Strongly recommended.** Your email address for NCBI Entrez. NCBI requires this for reliable access and to contact you if there are issues with your queries.
- `--entrez-api-key YOUR_API_KEY`: *(**Optional, but highly recommended**)* Your NCBI API key for higher request rates.

### Key Options:
- `-o, --output-dir DIR`: Output directory (default: protfetch_results).
- `--entrez-email EMAIL`: Strongly recommended for NCBI access.
- `--entrez-api-key KEY`: (Optional) NCBI API key.
- `--max-dist INT`: Max Levenshtein distance for internal filter (default: 4; 0 to disable).
- `--enable-cdhit, -ec`: Enable CD-HIT. Uses auto-thresholds unless `--cdhit-fixed-threshold` is set.
- `--cdhit-fixed-threshold FLOAT`: Use a fixed CD-HIT identity threshold (0.4-1.0). Implies CD-HIT is enabled.
- `--solo-redundancy-filtering, -srf`: Enable solo mode.
- `--input-fasta-srf PATH`: Input FASTA for -srf mode.
- `--input-csv-srf PATH`: Input CSV for -srf mode (needs 'identifier', 'gene' columns).
- `--cdhit-group-workers INT`: Parallel workers for CD-HIT in -srf mode (default: 6).
- `--max-workers INT`: Max concurrent workers for NCBI fetching (standard mode, default: 5).
- `--save-individual-files`: Save individual and intermediate files (raw, keyword-filtered, pre-CD-HIT).
- `--skip-keyword-filter`: Skip keyword-based FASTA filtering.
- `--debug`: Enable detailed debug logging.
- `-v, --version`: Show version.
- `-h, --help`: Show help.

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

Output files are saved in the specified output directory.

### Combined files:
- `*_combo_short.fasta`, `*_combo_full.fasta`, `*_combo_meta.csv`.

### Individual files (if `--save-individual-files`):
Saved in `output_dir/individual_gene_files/`.

- `GENE_0_raw_ncbi.fasta`: Raw sequences from NCBI.
- `GENE_1_keyword_filtered.fasta`: Sequences after keyword filtering.
- `GENE_1.5_pre_cdhit.fasta`: Sequences just before CD-HIT (if CD-HIT enabled).
- `GENE_2_final_short.fasta`, `GENE_2_final_full.fasta`, `GENE_2_final_meta.csv`: Final processed sequences.

### Solo mode (-srf) combined files:
Suffix `_srf_cdhit_combined` (e.g., `input_srf_cdhit_combined_short.fasta`).

### Solo mode (-srf) individual group files (if `--save-individual-files`):
Saved in `output_dir/individual_gene_files/GENE_GROUP_KEY/`.

- `GENE_GROUP_KEY_1.5_pre_cdhit.fasta`
- `GENE_GROUP_KEY_srf_final_short.fasta`, etc.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
