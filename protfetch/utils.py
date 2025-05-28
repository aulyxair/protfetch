# protfetch/utils.py
import logging
import sys
from pathlib import Path
from typing import List, Union

DEFAULT_MAX_DIST = 4
DEFAULT_ENTREZ_EMAIL = "your.email@example.com"
DEFAULT_REQUEST_TIMEOUT = 60
DEFAULT_REQUEST_RETRIES = 3
DEFAULT_MAX_WORKERS = 5

OUTPUT_SUBDIR_INDIVIDUAL = "individual_gene_files"
COMBINED_FASTA_SHORT_SUFFIX = "_combo_short.fasta"
COMBINED_FASTA_FULL_SUFFIX = "_combo_full.fasta"
COMBINED_CSV_SUFFIX = "_combo_meta.csv"


def setup_logging(level=logging.INFO):
    logger = logging.getLogger("protfetch")
    if logger.hasHandlers():
        logger.handlers.clear()
    logger.setLevel(level)
    logger.propagate = False
    ch = logging.StreamHandler(sys.stderr)
    ch.setLevel(level)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.debug(
        f"Logging for 'protfetch' initialized with level {logging.getLevelName(level)}."
    )
    return logger


log = logging.getLogger("protfetch")
if not log.handlers:
    log.addHandler(logging.NullHandler())


def ensure_output_dir(path_str: str) -> Path:
    path = Path(path_str)
    path.mkdir(parents=True, exist_ok=True)
    log.debug(f"Ensured output directory exists: {path}")
    return path


class GeneInput:
    def __init__(self, line: str):
        self.original_line = line.strip()
        self.gene_symbol: str = ""
        self.query_keyword: str = ""
        self.protein_name: Union[str, None] = None

        if not self.original_line:
            raise ValueError("Input line cannot be empty.")

        if "|" in self.original_line:
            parts = self.original_line.split("|", 1)
            self.protein_name = parts[0].strip()
            self.gene_symbol = parts[1].strip()
            self.query_keyword = self.protein_name
        else:
            self.gene_symbol = self.original_line
            self.query_keyword = self.gene_symbol

        if not self.gene_symbol:
            raise ValueError(
                f"Could not parse gene symbol from line: '{self.original_line}'"
            )
        if not self.query_keyword:
            raise ValueError(
                f"Could not determine query keyword for line: '{self.original_line}'"
            )

    def __repr__(self):
        return f"GeneInput(gene_symbol='{self.gene_symbol}', query_keyword='{self.query_keyword}', protein_name='{self.protein_name}')"


def parse_gene_list_file(file_path: str) -> List[GeneInput]:
    log.debug(f"Attempting to parse gene list file: {file_path}")
    genes_to_process: List[GeneInput] = []
    try:
        with open(file_path, "r") as f:
            for i, line_content in enumerate(f, 1):
                line_content = line_content.strip()
                log.debug(f"Reading line {i}: '{line_content}'")
                if not line_content or line_content.startswith("#"):
                    log.debug(f"Skipping line {i} (empty or comment).")
                    continue
                try:
                    gene_input = GeneInput(line_content)
                    genes_to_process.append(gene_input)
                    log.debug(f"Parsed line {i}: {gene_input}")
                except ValueError as e:
                    log.warning(
                        f"Skipping invalid line {i} in '{file_path}': {line_content}. Error: {e}"
                    )
        log.info(
            f"Successfully parsed {len(genes_to_process)} gene inputs from {file_path}."
        )
    except FileNotFoundError:
        log.error(f"Input gene list file not found: {file_path}")
        raise
    except Exception as e:
        log.error(f"Error reading gene list file '{file_path}': {e}", exc_info=True)
        raise
    return genes_to_process
