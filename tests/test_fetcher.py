from unittest.mock import MagicMock, patch

import pytest
from Bio import Entrez

from protfetch.fetcher import configure_entrez, fetch_protein_fasta_for_gene
from protfetch.utils import GeneInput


@pytest.fixture(autouse=True)
def auto_configure_entrez_for_tests():
    configure_entrez("ci_test@example.com")


@patch("protfetch.fetcher.time.sleep")
@patch("Bio.Entrez.read")
@patch("protfetch.fetcher.Entrez.efetch")
@patch("protfetch.fetcher.Entrez.elink")
@patch("protfetch.fetcher.Entrez.esearch")
def test_fetch_protein_fasta_for_gene_success(
    mock_esearch, mock_elink, mock_efetch, mock_entrez_read, mock_sleep
):
    mock_entrez_read.side_effect = [
        {"IdList": ["12345"]},
        [{"LinkSetDb": [{"Link": [{"Id": "NP_001"}, {"Id": "XP_002"}]}]}],
    ]

    mock_esearch_handle = MagicMock()
    mock_esearch.return_value = mock_esearch_handle

    mock_elink_handle = MagicMock()
    mock_elink.return_value = mock_elink_handle

    mock_efetch_handle = MagicMock()
    mock_efetch_handle.read.return_value = (
        ">NP_001 Protein1\nACGT\n>XP_002 Protein2\nTGCA"
    )
    mock_efetch.return_value = mock_efetch_handle

    gene_input = GeneInput("GENESYM")
    fasta_content = fetch_protein_fasta_for_gene(gene_input, timeout=10, retries=1)

    assert fasta_content is not None, "fasta_content should not be None"
    assert ">NP_001 Protein1" in fasta_content
    mock_esearch.assert_called_once()
    mock_elink.assert_called_once()
    mock_efetch.assert_called_once()
    assert mock_entrez_read.call_count == 2


@patch("protfetch.fetcher.time.sleep")
@patch("Bio.Entrez.read")
@patch("protfetch.fetcher.Entrez.esearch")
def test_fetch_protein_fasta_no_gene_uid(mock_esearch, mock_entrez_read, mock_sleep):
    mock_entrez_read.side_effect = [
        {"IdList": []},  # Result for the first esearch call (e.g., with [Gene Symbol])
        {"IdList": []},  # Result for the second esearch call (e.g., with [sym])
    ]

    mock_esearch_handle = MagicMock()
    mock_esearch.return_value = mock_esearch_handle

    gene_input = GeneInput("NONEXISTENT_GENE")
    fasta_content = fetch_protein_fasta_for_gene(gene_input, timeout=10, retries=1)
    assert fasta_content is None
    assert mock_esearch.call_count == 2
    assert mock_entrez_read.call_count == 2


@patch("protfetch.fetcher.time.sleep")
@patch("Bio.Entrez.read")
@patch(
    "protfetch.fetcher.Entrez.elink",
    side_effect=RuntimeError("NCBI elink failed intentionally"),
)
@patch("protfetch.fetcher.Entrez.esearch")
def test_entrez_retry_call_failure(
    mock_esearch, mock_elink_fails, mock_entrez_read, mock_sleep
):
    mock_entrez_read.return_value = {"IdList": ["12345"]}

    mock_esearch_handle = MagicMock()
    mock_esearch.return_value = mock_esearch_handle

    gene_input = GeneInput("FAILING_GENE")
    fasta_content = fetch_protein_fasta_for_gene(gene_input, timeout=10, retries=2)

    assert fasta_content is None
    assert mock_elink_fails.call_count == 2
