import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def temp_output_dir(tmp_path: Path) -> Path:
    """Provides a temporary directory for test outputs."""
    output_subdir = tmp_path / "test_outputs"
    output_subdir.mkdir(exist_ok=True)
    return output_subdir


@pytest.fixture
def sample_gene_list_content() -> str:
    """Provides sample content for a gene list file."""
    return "Protein Kinase A | PRKACA\nTP53\n# This is a comment\nBRCA1 | BRCA1_GENE"


@pytest.fixture
def sample_gene_list_file(tmp_path: Path, sample_gene_list_content: str) -> Path:
    """Creates a temporary sample gene list file and returns its path."""
    file_path = tmp_path / "sample_genes.txt"
    file_path.write_text(sample_gene_list_content)
    return file_path
