tests/conftest.py
import pytest
from pathlib import Path
import tempfile

@pytest.fixture
def temp_output_dir():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)

@pytest.fixture
def sample_gene_list_file(tmp_path):
    content = "Protein Kinase A | PRKACA\nTP53\n# This is a comment\nBRCA1 | BRCA1_HUMAN"
    file_path = tmp_path / "sample_genes.txt"
    file_path.write_text(content)
    return file_path

