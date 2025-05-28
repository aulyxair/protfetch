from pathlib import Path

import pytest

from protfetch.utils import GeneInput, parse_gene_list_file


def test_parse_gene_list_file_formats(sample_gene_list_file: Path):
    parsed_genes = parse_gene_list_file(str(sample_gene_list_file))
    assert len(parsed_genes) == 3

    assert isinstance(parsed_genes[0], GeneInput)
    assert parsed_genes[0].gene_symbol == "PRKACA"
    assert parsed_genes[0].query_keyword == "Protein Kinase A"
    assert parsed_genes[0].protein_name == "Protein Kinase A"

    assert isinstance(parsed_genes[1], GeneInput)
    assert parsed_genes[1].gene_symbol == "TP53"
    assert parsed_genes[1].query_keyword == "TP53"
    assert parsed_genes[1].protein_name is None

    assert isinstance(parsed_genes[2], GeneInput)
    assert parsed_genes[2].gene_symbol == "BRCA1_GENE"
    assert parsed_genes[2].query_keyword == "BRCA1"
    assert parsed_genes[2].protein_name == "BRCA1"
