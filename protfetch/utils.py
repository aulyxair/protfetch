# protfetch/utils.py
import logging
import sys
from pathlib import Path
from typing import List, Union

# --- Constants ---
DEFAULT_MAX_DIST = 4
DEFAULT_ENTREZ_EMAIL = "your.email@example.com"  # User should override
DEFAULT_REQUEST_TIMEOUT = 60  # seconds
DEFAULT_REQUEST_RETRIES = 3
DEFAULT_MAX_WORKERS = 5  # For concurrent fetching

OUTPUT_SUBDIR_INDIVIDUAL = "individual_gene_files"
COMBINED_FASTA_SHORT_SUFFIX = "_combo_short.fasta"
COMBINED_FASTA_FULL_SUFFIX = "_combo_full.fasta"
COMBINED_CSV_SUFFIX = "_combo_meta.csv"


# --- Logging Setup ---
def setup_logging(level=logging.INFO):
    """
    Sets up basic logging.
    """
    logger = logging.getLogger("protfetch")
    logger.setLevel(level)

    ch = logging.StreamHandler(sys.stderr)
    ch.setLevel(level)

    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    ch.setFormatter(formatter)

    if not logger.handlers:
        logger.addHandler(ch)
    return logger


log = logging.getLogger("protfetch")
if not log.handlers:
    log.addHandler(logging.NullHandler())


# --- File System Utilities ---
def ensure_output_dir(path_str: str) -> Path:
    """
    Ensures the output directory exists.
    Creates it if it doesn't.
    Returns a Path object.
    """
    path = Path(path_str)
    path.mkdir(parents=True, exist_ok=True)
    return path


# --- Input Parsing ---
class GeneInput:
    """
    Represents a single gene input, handling the two formats.
    """

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
    """
    Parses the input gene list file.
    Handles two formats:
    1. Protein Name | GENE_SYMBOL
    2. GENE_SYMBOL
    """
    genes_to_process: List[GeneInput] = []
    try:
        with open(file_path, "r") as f:
            for i, line_content in enumerate(f, 1):
                line_content = line_content.strip()
                if not line_content or line_content.startswith("#"):
                    continue
                try:
                    genes_to_process.append(GeneInput(line_content))
                except ValueError as e:
                    log.warning(
                        f"Skipping invalid line {i} in '{file_path}': {line_content}. Error: {e}"
                    )
    except FileNotFoundError:
        log.error(f"Input gene list file not found: {file_path}")
        raise
    except Exception as e:
        log.error(f"Error reading gene list file '{file_path}': {e}")
        raise
    return genes_to_process
