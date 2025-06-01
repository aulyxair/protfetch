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
        {  # Result for the first Entrez.esearch call
            "Count": "1",
            "IdList": ["12345"],
            "WebEnv": "SAMPLE_WEBENV_ESEARCH_GENE",
            "QueryKey": "1",
        },
        [  # Result for the Entrez.elink call
            {
                "LinkSetDb": [
                    {
                        "DbTo": "protein",
                        "LinkName": "gene_protein",
                        "Link": [{"Id": "NP_001"}, {"Id": "XP_002"}],
                    }
                ]
            }
        ],
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
        {  # Result for the first esearch call (e.g., with [symbol])
            "Count": "0",
            "IdList": [],
            "WebEnv": "WEBENV1",
            "QueryKey": "1",
        },
        {  # Result for the second esearch call (e.g., with [Gene Symbol])
            "Count": "0",
            "IdList": [],
            "WebEnv": "WEBENV2",
            "QueryKey": "2",
        },
    ]

    mock_esearch_handle = MagicMock()  # This handle will be used for both calls
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
    # Mock for esearch result - needs to include "Count", "WebEnv", "QueryKey"
    mock_entrez_read.return_value = {
        "Count": "1",
        "IdList": ["12345"],
        "WebEnv": "SAMPLE_WEBENV_ESEARCH_FAIL",
        "QueryKey": "1",
    }

    mock_esearch_handle = MagicMock()
    mock_esearch.return_value = mock_esearch_handle

    gene_input = GeneInput("FAILING_GENE")
    fasta_content = fetch_protein_fasta_for_gene(gene_input, timeout=10, retries=2)

    assert fasta_content is None
    assert mock_elink_fails.call_count == 2
    mock_esearch.assert_called_once()
    assert mock_entrez_read.call_count == 1
