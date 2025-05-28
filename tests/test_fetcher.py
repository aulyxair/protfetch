import pytest
from unittest.mock import patch, MagicMock
from protfetch.fetcher import fetch_protein_fasta_for_gene, configure_entrez
from protfetch.utils import GeneInput

@pytest.fixture(autouse=True)
def auto_configure_entrez():
    # Ensure Entrez is configured for all fetcher tests to avoid errors
    # In actual tests, you might want more control or to mock Entrez.email itself.
    configure_entrez("test@example.com")


@patch('protfetch.fetcher.Entrez.esearch')
@patch('protfetch.fetcher.Entrez.elink')
@patch('protfetch.fetcher.Entrez.efetch')
def test_fetch_protein_fasta_for_gene_success(mock_efetch, mock_elink, mock_esearch):
    # Mock Entrez calls to return expected data
    mock_esearch_handle = MagicMock()
    mock_esearch_handle.read.return_value = "<IdList><Id>12345</Id></IdList>" # Simplified XML or use Entrez.read result
    Entrez.read = MagicMock(return_value={'IdList': ['12345']}) # Mock Entrez.read for esearch
    mock_esearch.return_value = mock_esearch_handle

    mock_elink_handle = MagicMock()
    # Entrez.read for elink returns a list of LinkSetDb dictionaries
    Entrez.read.side_effect = [ # First call for esearch, second for elink
        {'IdList': ['12345']}, 
        [{'LinkSetDb': [{'Link': [{'Id': 'NP_001'}, {'Id': 'XP_002'}]}]}]
    ]
    mock_elink.return_value = mock_elink_handle
    
    mock_efetch_handle = MagicMock()
    mock_efetch_handle.read.return_value = ">NP_001 Protein1\nACGT\n>XP_002 Protein2\nTGCA"
    mock_efetch.return_value = mock_efetch_handle

    gene_input = GeneInput("GENESYM")
    fasta_content = fetch_protein_fasta_for_gene(gene_input, timeout=10, retries=1)

    assert fasta_content is not None
    assert ">NP_001 Protein1" in fasta_content
    assert "ACGT" in fasta_content
    mock_esearch.assert_called_once()
    mock_elink.assert_called_once()
    mock_efetch.assert_called_once()


@patch('protfetch.fetcher.Entrez.esearch')
def test_fetch_protein_fasta_no_gene_uid(mock_esearch):
    mock_esearch_handle = MagicMock()
    Entrez.read = MagicMock(return_value={'IdList': []}) # No IDs found
    mock_esearch.return_value = mock_esearch_handle

    gene_input = GeneInput("NONEXISTENT_GENE")
    fasta_content = fetch_protein_fasta_for_gene(gene_input, timeout=10, retries=1)
    assert fasta_content is None
