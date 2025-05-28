# protfetch/processor.py
import csv
from collections import defaultdict
import time
from io import StringIO
from typing import TextIO, Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rapidfuzz.distance import Levenshtein # Now a core dependency

from .utils import log # Use the centralized logger

# --- Data Structures ---
class ProcessedProtein:
    """
    Holds data for a single processed protein sequence.
    """
    def __init__(self, accession: str, gene_symbol_input: str, identifier_from_header: str,
                 sequence: str, full_header: str):
        self.accession = accession
        self.gene_symbol_input = gene_symbol_input.lower() # Store gene symbol in lowercase as in original
        self.identifier_from_header = identifier_from_header # This is the 'species' or specific ID from header
        self.sequence = sequence
        self.full_header = full_header.strip()

    def get_short_header_line(self) -> str:
        """
        Generates the short FASTA header line (e.g., >P12345).
        """
        return f">{self.accession}"

    def get_full_header_line(self) -> str:
        """
        Generates the full FASTA header line (e.g., >sp|P12345|GENE_HUMAN Protein description).
        """
        return f">{self.full_header}"

    def get_csv_row(self) -> list[str]:
        """
        Returns a list of strings for the metadata CSV row.
        Corresponds to ['identifier', 'gene', 'species'] from original script.
        Here, 'identifier' is accession, 'gene' is input gene_symbol, 'species' is identifier_from_header.
        """
        return [self.accession, self.gene_symbol_input, self.identifier_from_header]

    def __repr__(self):
        return f"ProcessedProtein(accession='{self.accession}', gene='{self.gene_symbol_input}', seq_len={len(self.sequence)})"

# --- FASTA Header Parsing Logic (adapted from original script) ---
def _parse_fasta_header(header_content: str) -> tuple[str | None, str | None, str]:
    """
    Parses a FASTA header line to extract accession and identifier.
    Returns (accession, identifier, original_full_header_content).
    'identifier' is often a species tag or a more specific ID.
    """
    header_content = header_content.strip()
    accession: str | None = None
    identifier: str | None = None

    if not header_content:
        log.warning(f"Empty header content found.")
        return None, None, header_content

    # Try UniProt style: sp|ACCESSION|ID or tr|ACCESSION|ID
    if header_content.startswith('sp|') or header_content.startswith('tr|'):
        parts = header_content.split('|')
        if len(parts) >= 2:
            accession = parts[1].strip()
            if len(parts) >= 3:
                identifier = parts[2].split(None, 1)[0].strip() # First part of the description field
            else: # Should ideally not happen for well-formed UniProt
                identifier = accession
        else:
            log.warning(f"Could not parse UniProt-style header: {header_content}")
            # Fallback: try to get accession from the first part
            accession = header_content.split(None, 1)[0].strip() # Get the first word as accession
            identifier = accession # Default identifier to accession if specific parsing fails
    else:
        # Try general NCBI or other styles: >ACCESSION [Organism Name] or >ACCESSION description
        parts = header_content.split(None, 1)
        accession = parts[0].strip()

        # Try to find identifier in square brackets (e.g., [Homo sapiens])
        start_bracket = header_content.rfind('[')
        end_bracket = header_content.rfind(']')
        if start_bracket != -1 and end_bracket != -1 and end_bracket > start_bracket:
            identifier = header_content[start_bracket + 1 : end_bracket].replace(' ', '_').strip()
        else:
            # If no brackets, use accession as identifier, or try to get it from description if available
            if len(parts) > 1:
                # A simple heuristic: if the description contains the gene symbol pattern (e.g., _HUMAN, _MOUSE)
                # This is very heuristic and might need refinement.
                # For now, let's stick to a simpler identifier if not in brackets.
                # The original script used accession if brackets not found.
                identifier = accession # Default to accession if no specific identifier found
            else:
                identifier = accession

    if not accession: # Should be caught by earlier checks, but as a safeguard
        log.warning(f"Failed to parse accession from header: {header_content}")
        # Fallback to using the whole header (without '>') as accession if desperate
        accession = header_content.split(None,1)[0] if header_content else "unknown_accession"
        identifier = identifier or accession


    return accession, identifier, header_content


# --- Core FASTA Processing and Filtering ---
def process_fasta_stream(
    fasta_stream: TextIO,
    input_gene_symbol: str,
    max_levenshtein_distance: int = 4
) -> tuple[list[ProcessedProtein], dict[str, Any]]:
    """
    Processes a FASTA stream: parses, filters, and returns processed protein data.
    This function encapsulates the logic from the original `parse_fasta_to_dict` and `filter_sequences`.

    Args:
        fasta_stream: A file-like object (e.g., open file, StringIO) containing FASTA data.
        input_gene_symbol: The gene symbol associated with this FASTA data.
        max_levenshtein_distance: Max Levenshtein distance for near-identical filtering. 0 to disable.

    Returns:
        A tuple: (list_of_ProcessedProtein_objects, statistics_dict)
    """
    log.info(f"Processing FASTA for gene '{input_gene_symbol}'...")
    stats = {
        "headers_encountered": 0,
        "parsing_skipped_or_incomplete": 0,
        "duplicates_accession_skipped": 0,
        "initial_unique_sequences_parsed": 0,
        "removed_identical_sequences": 0,
        "removed_near_identical_sequences": 0,
        "removed_fragment_sequences": 0,
        "final_sequences_kept": 0,
    }

    # Step 1: Parse FASTA and initial duplicate accession removal
    # Store as: accession -> ProcessedProtein object
    parsed_proteins_map: dict[str, ProcessedProtein] = {}
    seen_accessions_in_file = set()

    for record in SeqIO.parse(fasta_stream, "fasta"):
        stats["headers_encountered"] += 1
        raw_header_content = record.description # Biopython's record.description is the full header line after '>'
        
        accession, identifier_from_header, full_header_for_storage = _parse_fasta_header(raw_header_content)

        if not accession or not identifier_from_header:
            log.warning(f"Gene {input_gene_symbol}: Skipping header due to parsing issue: '{raw_header_content}'")
            stats["parsing_skipped_or_incomplete"] += 1
            continue
        
        sequence_str = str(record.seq).upper() # Ensure sequence is uppercase
        if not sequence_str:
            log.warning(f"Gene {input_gene_symbol}: Skipping accession '{accession}' due to empty sequence.")
            stats["parsing_skipped_or_incomplete"] += 1
            continue

        if accession in seen_accessions_in_file:
            log.debug(f"Gene {input_gene_symbol}: Duplicate accession '{accession}' in this file. Keeping first seen.")
            stats["duplicates_accession_skipped"] += 1
            continue # Keep the first one encountered with this accession
        
        seen_accessions_in_file.add(accession)
        parsed_proteins_map[accession] = ProcessedProtein(
            accession=accession,
            gene_symbol_input=input_gene_symbol,
            identifier_from_header=identifier_from_header,
            sequence=sequence_str,
            full_header=full_header_for_storage
        )

    stats["initial_unique_sequences_parsed"] = len(parsed_proteins_map)
    if not parsed_proteins_map:
        log.info(f"Gene {input_gene_symbol}: No valid sequences parsed from input.")
        return [], stats

    # Step 2: Filtering
    # Convert map to list of ProcessedProtein objects for easier filtering
    current_proteins = list(parsed_proteins_map.values())

    # Filter 1: Identical sequences (different accessions, same sequence)
    # Keep the one with the lexicographically smallest accession
    log.debug(f"Gene {input_gene_symbol}: [Filter 1/3] Identifying identical sequences among {len(current_proteins)} entries...")
    start_time_filter1 = time.time()
    sequence_to_proteins_map = defaultdict(list)
    for protein in current_proteins:
        sequence_to_proteins_map[protein.sequence].append(protein)
    
    proteins_after_identity_check = []
    for _, protein_group in sequence_to_proteins_map.items():
        if len(protein_group) > 1:
            protein_group.sort(key=lambda p: p.accession) # Sort by accession
            kept_protein = protein_group[0]
            proteins_after_identity_check.append(kept_protein)
            stats["removed_identical_sequences"] += (len(protein_group) - 1)
            log.debug(f"Gene {input_gene_symbol}: Kept '{kept_protein.accession}', removed {len(protein_group)-1} identical sequence proteins: {[p.accession for p in protein_group[1:]]}")
        else:
            proteins_after_identity_check.append(protein_group[0])
    current_proteins = proteins_after_identity_check
    log.debug(f"Gene {input_gene_symbol}:  - Removed {stats['removed_identical_sequences']} identical sequences. Remaining: {len(current_proteins)}. Time: {time.time() - start_time_filter1:.2f}s")


    # Filter 2: Near-identical sequences (Levenshtein distance)
    # Keep longer sequences, or lexicographically first accession for ties.
    if max_levenshtein_distance > 0:
        log.debug(f"Gene {input_gene_symbol}: [Filter 2/3] Identifying near-identical (Levenshtein <= {max_levenshtein_distance}) among {len(current_proteins)}...")
        start_time_filter2 = time.time()
        # Sort by sequence length (descending), then by accession (ascending) for tie-breaking
        current_proteins.sort(key=lambda p: (-len(p.sequence), p.accession))
        
        representatives: list[ProcessedProtein] = []
        proteins_to_remove_near_identity = set() # Set of accessions to remove

        for protein in current_proteins:
            if protein.accession in proteins_to_remove_near_identity:
                continue
            is_near_duplicate = False
            for rep_protein in representatives:
                # Optimization: if length difference is already > max_dist, skip Levenshtein
                if abs(len(protein.sequence) - len(rep_protein.sequence)) > max_levenshtein_distance:
                    continue
                
                distance = Levenshtein.distance(protein.sequence, rep_protein.sequence, score_cutoff=max_levenshtein_distance)
                if distance <= max_levenshtein_distance:
                    # protein is near-duplicate of rep_protein. Since current_proteins is sorted by length (desc)
                    # rep_protein is either longer or same length but earlier accession.
                    # So, we remove the current 'protein'.
                    proteins_to_remove_near_identity.add(protein.accession)
                    is_near_duplicate = True
                    log.debug(f"Gene {input_gene_symbol}: Near-Identical: Removing '{protein.accession}' (dist {distance} to '{rep_protein.accession}')")
                    break
            if not is_near_duplicate:
                representatives.append(protein)
        
        num_removed = len(proteins_to_remove_near_identity)
        stats["removed_near_identical_sequences"] = num_removed
        current_proteins = [p for p in current_proteins if p.accession not in proteins_to_remove_near_identity]
        log.debug(f"Gene {input_gene_symbol}:  - Removed {num_removed} near-identical. Remaining: {len(current_proteins)}. Time: {time.time() - start_time_filter2:.2f}s")
    else:
        log.debug(f"Gene {input_gene_symbol}: [Filter 2/3] Skipping near-identical check (max_dist=0).")


    # Filter 3: Fragment sequences (substrings)
    # If seqA is a substring of seqB, remove seqA.
    log.debug(f"Gene {input_gene_symbol}: [Filter 3/3] Identifying fragment sequences among {len(current_proteins)}...")
    start_time_filter3 = time.time()
    # Sort by length (descending) so longer sequences are checked first as potential containers
    current_proteins.sort(key=lambda p: (-len(p.sequence), p.accession))
    
    proteins_to_remove_fragments = set() # Set of accessions to remove
    n_candidates = len(current_proteins)

    for i in range(n_candidates):
        protein1 = current_proteins[i]
        if protein1.accession in proteins_to_remove_fragments:
            continue
        
        for j in range(n_candidates): # Compare with all others (including those already processed as protein1)
            if i == j:
                continue
            protein2 = current_proteins[j]
            if protein2.accession in proteins_to_remove_fragments:
                continue

            # If protein2.sequence is a substring of protein1.sequence, and they are different,
            # and protein2 is shorter, then protein2 is a fragment.
            if len(protein2.sequence) < len(protein1.sequence) and protein2.sequence in protein1.sequence:
                proteins_to_remove_fragments.add(protein2.accession)
                log.debug(f"Gene {input_gene_symbol}: Fragment: Removing '{protein2.accession}' (substring of '{protein1.accession}')")

    num_removed_frags = len(proteins_to_remove_fragments)
    stats["removed_fragment_sequences"] = num_removed_frags
    final_proteins_kept = [p for p in current_proteins if p.accession not in proteins_to_remove_fragments]
    log.debug(f"Gene {input_gene_symbol}:  - Removed {num_removed_frags} fragments. Remaining: {len(final_proteins_kept)}. Time: {time.time() - start_time_filter3:.2f}s")

    # Sort final results by accession for consistent output
    final_proteins_kept.sort(key=lambda p: p.accession)
    stats["final_sequences_kept"] = len(final_proteins_kept)

    log.info(f"Gene {input_gene_symbol}: Processing complete. Kept {stats['final_sequences_kept']} sequences.")
    return final_proteins_kept, stats


# --- Output Writing Helpers (used by main CLI module) ---
def write_processed_proteins_to_fasta(
    proteins: list[ProcessedProtein],
    output_file_path: str,
    use_full_header: bool
):
    """Writes a list of ProcessedProtein objects to a FASTA file."""
    try:
        with open(output_file_path, "w") as outfile:
            for protein in proteins:
                header_line = protein.get_full_header_line() if use_full_header else protein.get_short_header_line()
                outfile.write(f"{header_line}\n")
                # Write sequence in lines of, e.g., 60 characters
                seq = protein.sequence
                for i in range(0, len(seq), 60):
                    outfile.write(seq[i:i+60] + "\n")
        log.debug(f"FASTA ({'full' if use_full_header else 'short'} header) saved to: {output_file_path}")
    except IOError as e:
        log.error(f"Error writing FASTA file {output_file_path}: {e}")
        raise

def write_processed_proteins_to_csv(
    proteins: list[ProcessedProtein],
    output_file_path: str
):
    """Writes metadata from a list of ProcessedProtein objects to a CSV file."""
    try:
        with open(output_file_path, "w", newline='') as outfile_csv:
            csv_writer = csv.writer(outfile_csv)
            csv_writer.writerow(['accession', 'gene_input', 'identifier_from_header']) # Header row
            for protein in proteins:
                csv_writer.writerow(protein.get_csv_row())
        log.debug(f"Metadata CSV saved to: {output_file_path}")
    except IOError as e:
        log.error(f"Error writing CSV file {output_file_path}: {e}")
        raise

# --- Deduplication for combined files ---
def deduplicate_processed_proteins(
    proteins: list[ProcessedProtein],
    dedup_key: str = "accession" # 'accession' or 'full_header'
) -> list[ProcessedProtein]:
    """
    Deduplicates a list of ProcessedProtein objects based on a key.
    Keeps the first encountered.
    """
    seen_keys = set()
    deduplicated_list = []
    for protein in proteins:
        key_val = ""
        if dedup_key == "accession":
            key_val = protein.accession
        elif dedup_key == "full_header":
            key_val = protein.full_header
        else:
            raise ValueError(f"Invalid dedup_key: {dedup_key}. Must be 'accession' or 'full_header'.")

        if key_val not in seen_keys:
            deduplicated_list.append(protein)
            seen_keys.add(key_val)
        else:
            log.debug(f"Deduplication: Removing protein with {dedup_key} '{key_val}' as it was already seen.")
    return deduplicated_list


