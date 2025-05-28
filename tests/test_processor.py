from io import StringIO
from protfetch.processor import process_fasta_stream, ProcessedProtein

MOCK_FASTA_CONTENT = """>sp|P01111|RAS_HUMAN GTPase HRas OS=Homo sapiens OX=9606 GN=HRAS PE=1 SV=1
MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMNCKCVLS
>tr|A0A023|FAKE_GENE Hypothetical protein OS=Mus musculus OX=10090 GN=Fake PE=4 SV=1
MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMNCKCVLS
>sp|P01112|RAS_MOUSE GTPase HRas OS=Mus musculus OX=10090 GN=Hras PE=1 SV=1
MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLPSRTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMNCKCVLS
>fragment|FRAG01|FRAG_HUMAN Short fragment
MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMNCKCVLSKLNPPDESGPGCMNCKCVLS
>fragment_actual|FRAG02|FRAG_HUMAN Shorter fragment
SAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMNCKCVLS
"""

def test_process_fasta_stream_basic():
    fasta_stream = StringIO(MOCK_FASTA_CONTENT)
    gene_symbol = "TESTGENE"
    max_dist = 4 # Levenshtein distance
    
    # This is a placeholder, actual test would need more specific assertions
    # and potentially mocking of rapidfuzz if it's slow or complex to test.
    proteins, stats = process_fasta_stream(fasta_stream, gene_symbol, max_dist)
    
    assert isinstance(proteins, list)
    assert all(isinstance(p, ProcessedProtein) for p in proteins)
    assert "final_sequences_kept" in stats
    
    # Example: Check if identical sequence filter worked (P01111 vs A0A023)
    # P01111 should be kept as its accession is lexicographically smaller than A0A023 for identical seq.
    accessions_kept = {p.accession for p in proteins}
    # Based on the logic:
    # P01111 and A0A023 are identical. P01111 is kept.
    # P01112 is slightly different (near identical to P01111). If dist <= 4, P01112 might be removed.
    # FRAG02 is a substring of P01111, so FRAG02 should be removed.
    # FRAG01 is longer than P01111, so P01111 is not a fragment of FRAG01.
    # This needs careful stepping through the filter logic.

    # For a real test, you'd construct input where you know the exact output.
    # E.g., assert "P01111" in accessions_kept
    # assert "A0A023" not in accessions_kept (if identical filter works as expected)
    # assert "FRAG02" not in accessions_kept (if fragment filter works)
