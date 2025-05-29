# protfetch/processor.py
import csv
import time
from collections import defaultdict
from io import StringIO
from typing import (  # Added Union, Tuple, List, Dict, Set for older Python compatibility
    Any,
    Dict,
    List,
    Set,
    TextIO,
    Tuple,
    Union,
)

from Bio import SeqIO

# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
from rapidfuzz.distance import Levenshtein

from .utils import log


class ProcessedProtein:
    """
    Holds data for a single processed protein sequence.
    """

    def __init__(
        self,
        accession: str,
        gene_symbol_input: str,
        identifier_from_header: str,
        sequence: str,
        full_header: str,
    ):
        self.accession = accession
        self.gene_symbol_input = gene_symbol_input.lower()
        self.identifier_from_header = identifier_from_header
        self.sequence = sequence
        self.full_header = full_header.strip()

    def get_short_header_line(self) -> str:
        return f">{self.accession}"

    def get_full_header_line(self) -> str:
        return f">{self.full_header}"

    def get_csv_row(self) -> List[str]:  # Changed to List[str]
        return [self.accession, self.gene_symbol_input, self.identifier_from_header]

    def __repr__(self):
        return f"ProcessedProtein(accession='{self.accession}', gene='{self.gene_symbol_input}', seq_len={len(self.sequence)})"


def _parse_fasta_header(
    header_content: str,
) -> Tuple[Union[str, None], Union[str, None], str]:
    """
    Parses a FASTA header line to extract accession and identifier.
    Returns (accession, identifier, original_full_header_content).
    'identifier' is often a species tag or a more specific ID.
    """
    header_content = header_content.strip()
    accession: Union[str, None] = None
    identifier: Union[str, None] = None

    if not header_content:
        log.warning(f"Empty header content found.")
        return None, None, header_content

    if header_content.startswith("sp|") or header_content.startswith("tr|"):
        parts = header_content.split("|")
        if len(parts) >= 2:
            accession = parts[1].strip()
            if len(parts) >= 3:
                identifier = parts[2].split(None, 1)[0].strip()
            else:
                identifier = accession
        else:
            log.warning(f"Could not parse UniProt-style header: {header_content}")
            accession = header_content.split(None, 1)[0].strip()
            identifier = accession
    else:
        parts = header_content.split(None, 1)
        accession = parts[0].strip()

        start_bracket = header_content.rfind("[")
        end_bracket = header_content.rfind("]")
        if start_bracket != -1 and end_bracket != -1 and end_bracket > start_bracket:
            identifier = (
                header_content[start_bracket + 1 : end_bracket]
                .replace(" ", "_")
                .strip()
            )
        else:
            identifier = accession

    if not accession:
        log.warning(f"Failed to parse accession from header: {header_content}")
        accession = (
            header_content.split(None, 1)[0] if header_content else "unknown_accession"
        )
        identifier = identifier or accession

    return accession, identifier, header_content


def process_fasta_stream(
    fasta_stream: TextIO, input_gene_symbol: str, max_levenshtein_distance: int = 4
) -> Tuple[List[ProcessedProtein], Dict[str, Any]]:
    log.info(f"Processing FASTA for gene '{input_gene_symbol}'...")
    stats: Dict[str, Any] = {
        "headers_encountered": 0,
        "parsing_skipped_or_incomplete": 0,
        "duplicates_accession_skipped": 0,
        "initial_unique_sequences_parsed": 0,
        "removed_identical_sequences": 0,
        "removed_near_identical_sequences": 0,
        "removed_fragment_sequences": 0,
        "final_sequences_kept": 0,
    }

    parsed_proteins_map: Dict[str, ProcessedProtein] = {}
    seen_accessions_in_file: Set[str] = set()

    for record in SeqIO.parse(fasta_stream, "fasta"):
        stats["headers_encountered"] += 1
        raw_header_content = record.description

        accession, identifier_from_header, full_header_for_storage = (
            _parse_fasta_header(raw_header_content)
        )

        if not accession or not identifier_from_header:
            log.warning(
                f"Gene {input_gene_symbol}: Skipping header due to parsing issue: '{raw_header_content}'"
            )
            stats["parsing_skipped_or_incomplete"] += 1
            continue

        sequence_str = str(record.seq).upper()
        if not sequence_str:
            log.warning(
                f"Gene {input_gene_symbol}: Skipping accession '{accession}' due to empty sequence."
            )
            stats["parsing_skipped_or_incomplete"] += 1
            continue

        if accession in seen_accessions_in_file:
            log.debug(
                f"Gene {input_gene_symbol}: Duplicate accession '{accession}' in this file. Keeping first seen."
            )
            stats["duplicates_accession_skipped"] += 1
            continue

        seen_accessions_in_file.add(accession)
        parsed_proteins_map[accession] = ProcessedProtein(
            accession=accession,
            gene_symbol_input=input_gene_symbol,
            identifier_from_header=identifier_from_header,
            sequence=sequence_str,
            full_header=full_header_for_storage,
        )

    stats["initial_unique_sequences_parsed"] = len(parsed_proteins_map)
    log.info(
        f"Gene {input_gene_symbol}: Initial parsing yielded {stats['initial_unique_sequences_parsed']} unique accessions from {stats['headers_encountered']} input records."
    )

    if not parsed_proteins_map:
        log.info(f"Gene {input_gene_symbol}: No valid sequences parsed from input.")
        return [], stats

    current_proteins: List[ProcessedProtein] = list(parsed_proteins_map.values())

    log.debug(
        f"Gene {input_gene_symbol}: [Filter 1/3] Identifying identical sequences among {len(current_proteins)} entries..."
    )
    start_time_filter1 = time.time()
    sequence_to_proteins_map: Dict[str, List[ProcessedProtein]] = defaultdict(list)
    for protein in current_proteins:
        sequence_to_proteins_map[protein.sequence].append(protein)

    proteins_after_identity_check: List[ProcessedProtein] = []
    for _, protein_group in sequence_to_proteins_map.items():
        if len(protein_group) > 1:
            protein_group.sort(key=lambda p: p.accession)
            kept_protein = protein_group[0]
            proteins_after_identity_check.append(kept_protein)
            stats["removed_identical_sequences"] += len(protein_group) - 1
            log.debug(
                f"Gene {input_gene_symbol}: Kept '{kept_protein.accession}', removed {len(protein_group)-1} identical sequence proteins: {[p.accession for p in protein_group[1:]]}"
            )
        else:
            proteins_after_identity_check.append(protein_group[0])
    current_proteins = proteins_after_identity_check
    log.debug(
        f"Gene {input_gene_symbol}:  - Removed {stats['removed_identical_sequences']} identical sequences. Remaining: {len(current_proteins)}. Time: {time.time() - start_time_filter1:.2f}s"
    )

    if max_levenshtein_distance > 0:
        log.debug(
            f"Gene {input_gene_symbol}: [Filter 2/3] Identifying near-identical (Levenshtein <= {max_levenshtein_distance}) among {len(current_proteins)}..."
        )
        start_time_filter2 = time.time()
        current_proteins.sort(key=lambda p: (-len(p.sequence), p.accession))

        representatives: List[ProcessedProtein] = []
        proteins_to_remove_near_identity: Set[str] = set()

        for protein in current_proteins:
            if protein.accession in proteins_to_remove_near_identity:
                continue
            is_near_duplicate = False
            for rep_protein in representatives:
                if (
                    abs(len(protein.sequence) - len(rep_protein.sequence))
                    > max_levenshtein_distance
                ):
                    continue

                distance = Levenshtein.distance(
                    protein.sequence,
                    rep_protein.sequence,
                    score_cutoff=max_levenshtein_distance,
                )
                if distance <= max_levenshtein_distance:
                    proteins_to_remove_near_identity.add(protein.accession)
                    is_near_duplicate = True
                    log.debug(
                        f"Gene {input_gene_symbol}: Near-Identical: Removing '{protein.accession}' (dist {distance} to '{rep_protein.accession}')"
                    )
                    break
            if not is_near_duplicate:
                representatives.append(protein)

        num_removed = len(proteins_to_remove_near_identity)
        stats["removed_near_identical_sequences"] = num_removed
        current_proteins = [
            p
            for p in current_proteins
            if p.accession not in proteins_to_remove_near_identity
        ]
        log.debug(
            f"Gene {input_gene_symbol}:  - Removed {num_removed} near-identical. Remaining: {len(current_proteins)}. Time: {time.time() - start_time_filter2:.2f}s"
        )
    else:
        log.debug(
            f"Gene {input_gene_symbol}: [Filter 2/3] Skipping near-identical check (max_dist=0)."
        )

    log.debug(
        f"Gene {input_gene_symbol}: [Filter 3/3] Identifying fragment sequences among {len(current_proteins)}..."
    )
    start_time_filter3 = time.time()
    current_proteins.sort(key=lambda p: (-len(p.sequence), p.accession))

    proteins_to_remove_fragments: Set[str] = set()
    n_candidates = len(current_proteins)

    for i in range(n_candidates):
        protein1 = current_proteins[i]
        if protein1.accession in proteins_to_remove_fragments:
            continue

        for j in range(n_candidates):
            if i == j:
                continue
            protein2 = current_proteins[j]
            if protein2.accession in proteins_to_remove_fragments:
                continue

            if (
                len(protein2.sequence) < len(protein1.sequence)
                and protein2.sequence in protein1.sequence
            ):
                proteins_to_remove_fragments.add(protein2.accession)
                log.debug(
                    f"Gene {input_gene_symbol}: Fragment: Removing '{protein2.accession}' (substring of '{protein1.accession}')"
                )

    num_removed_frags = len(proteins_to_remove_fragments)
    stats["removed_fragment_sequences"] = num_removed_frags
    final_proteins_kept: List[ProcessedProtein] = [
        p for p in current_proteins if p.accession not in proteins_to_remove_fragments
    ]
    log.debug(
        f"Gene {input_gene_symbol}:  - Removed {num_removed_frags} fragments. Remaining: {len(final_proteins_kept)}. Time: {time.time() - start_time_filter3:.2f}s"
    )

    final_proteins_kept.sort(key=lambda p: p.accession)
    stats["final_sequences_kept"] = len(final_proteins_kept)

    log.info(
        f"Gene {input_gene_symbol}: Processing complete. Started with {stats['headers_encountered']} records (from keyword-filtered input), kept {stats['final_sequences_kept']} sequences after all filters."
    )
    return final_proteins_kept, stats


def write_processed_proteins_to_fasta(
    proteins: List[ProcessedProtein], output_file_path: str, use_full_header: bool
):
    try:
        with open(output_file_path, "w") as outfile:
            for protein in proteins:
                header_line = (
                    protein.get_full_header_line()
                    if use_full_header
                    else protein.get_short_header_line()
                )
                outfile.write(f"{header_line}\n")
                seq = protein.sequence
                for i in range(0, len(seq), 60):
                    outfile.write(seq[i : i + 60] + "\n")
        log.debug(
            f"FASTA ({'full' if use_full_header else 'short'} header) saved to: {output_file_path}"
        )
    except IOError as e:
        log.error(f"Error writing FASTA file {output_file_path}: {e}")
        raise


def write_processed_proteins_to_csv(
    proteins: List[ProcessedProtein], output_file_path: str
):
    try:
        with open(output_file_path, "w", newline="") as outfile_csv:
            csv_writer = csv.writer(outfile_csv)
            csv_writer.writerow(["identifier", "gene", "species"])
            for protein in proteins:
                csv_writer.writerow(protein.get_csv_row())
        log.debug(f"Metadata CSV saved to: {output_file_path}")
    except IOError as e:
        log.error(f"Error writing CSV file {output_file_path}: {e}")
        raise


def deduplicate_processed_proteins(
    proteins: List[ProcessedProtein], dedup_key: str = "accession"
) -> List[ProcessedProtein]:
    seen_keys: Set[str] = set()
    deduplicated_list: List[ProcessedProtein] = []
    for protein in proteins:
        key_val: str = ""
        if dedup_key == "accession":
            key_val = protein.accession
        elif dedup_key == "full_header":
            key_val = protein.full_header
        else:
            raise ValueError(
                f"Invalid dedup_key: {dedup_key}. Must be 'accession' or 'full_header'."
            )

        if key_val not in seen_keys:
            deduplicated_list.append(protein)
            seen_keys.add(key_val)
        else:
            log.debug(
                f"Deduplication: Removing protein with {dedup_key} '{key_val}' as it was already seen."
            )
    return deduplicated_list
