# protfetch/main.py
import argparse
import csv
import logging
import os
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from . import __version__
from .cdhit_manager import (
    get_cdhit_threshold,
    is_cdhit_installed,
    run_cdhit_on_proteins,
)
from .fetcher import (
    configure_entrez,
    fetch_protein_fasta_for_gene,
    filter_fasta_by_keyword,
)
from .processor import (
    ProcessedProtein,
    deduplicate_processed_proteins,
    process_fasta_stream,
    write_processed_proteins_to_csv,
    write_processed_proteins_to_fasta,
)
from .utils import (
    COMBINED_CSV_SUFFIX,
    COMBINED_FASTA_FULL_SUFFIX,
    COMBINED_FASTA_SHORT_SUFFIX,
    DEFAULT_ENTREZ_EMAIL,
    DEFAULT_MAX_DIST,
    DEFAULT_MAX_WORKERS,
    DEFAULT_REQUEST_RETRIES,
    DEFAULT_REQUEST_TIMEOUT,
    OUTPUT_SUBDIR_INDIVIDUAL,
    GeneInput,
    ensure_output_dir,
    log,
    parse_gene_list_file,
    setup_logging,
)


def apply_cdhit_processing(
    proteins_to_filter: List[ProcessedProtein],
    gene_symbol: str,
    enable_cdhit: bool,
    fixed_cdhit_threshold: Optional[float],
    output_dir_individual: Optional[Path],
) -> List[ProcessedProtein]:
    if not enable_cdhit:
        return proteins_to_filter

    if not proteins_to_filter:
        log.info(
            f"CD-HIT: No proteins to process with CD-HIT for gene group '{gene_symbol}'."
        )
        return []

    if output_dir_individual:
        try:
            pre_cdhit_fasta_path = (
                output_dir_individual / f"{gene_symbol}_1.5_pre_cdhit.fasta"
            )
            seq_records = [
                SeqRecord(Seq(p.sequence), id=p.accession, description=p.full_header)
                for p in proteins_to_filter
            ]
            SeqIO.write(seq_records, pre_cdhit_fasta_path, "fasta")
            log.debug(
                f"Saved pre-CD-HIT FASTA for {gene_symbol} ({len(proteins_to_filter)} seqs) to {pre_cdhit_fasta_path}"
            )
        except Exception as e_save_pre_cdhit:
            log.error(
                f"Error saving pre-CD-HIT FASTA for {gene_symbol}: {e_save_pre_cdhit}"
            )

    threshold_to_use: float
    if fixed_cdhit_threshold is not None:
        threshold_to_use = fixed_cdhit_threshold
        log.info(
            f"CD-HIT: Using fixed threshold {threshold_to_use:.2f} for gene group '{gene_symbol}'."
        )
    else:
        count = len(proteins_to_filter)
        threshold_to_use = get_cdhit_threshold(count)
        log.info(
            f"CD-HIT: Applying CD-HIT for gene group '{gene_symbol}' (initial count: {count}, auto-threshold: {threshold_to_use:.2f})."
        )

    cdhit_filtered_proteins = run_cdhit_on_proteins(
        proteins_to_filter, threshold_to_use, gene_symbol
    )
    log.info(
        f"CD-HIT: After CD-HIT for gene group '{gene_symbol}', {len(cdhit_filtered_proteins)} proteins remain."
    )
    return cdhit_filtered_proteins


def process_single_gene_task(
    gene_input: GeneInput,
    output_dir_individual: Union[Path, None],
    max_dist: int,
    entrez_timeout: int,
    entrez_retries: int,
    skip_keyword_filter: bool,
    enable_cdhit: bool,
    fixed_cdhit_threshold: Optional[float],
) -> Union[Tuple[str, List[ProcessedProtein], Dict[str, Any]], None]:
    gene_symbol = gene_input.gene_symbol

    log.info(
        f"Starting processing for gene: {gene_symbol} (Keyword: '{gene_input.query_keyword}')"
    )

    stats: Dict[str, Any] = {"gene_symbol": gene_symbol}

    raw_fasta_content = fetch_protein_fasta_for_gene(
        gene_input, entrez_timeout, entrez_retries
    )
    if not raw_fasta_content:
        log.error(f"Failed to fetch FASTA data for gene {gene_symbol}.")
        stats["error"] = "Fetch failed"
        stats["status"] = "Fetch failed"
        return gene_symbol, [], stats

    if output_dir_individual and raw_fasta_content:
        try:
            raw_ncbi_path = output_dir_individual / f"{gene_symbol}_0_raw_ncbi.fasta"
            with open(raw_ncbi_path, "w") as f_raw:
                f_raw.write(raw_fasta_content)
            log.debug(f"Saved raw NCBI FASTA for {gene_symbol} to {raw_ncbi_path}")
        except Exception as e_save_raw:
            log.error(f"Error saving raw NCBI FASTA for {gene_symbol}: {e_save_raw}")

    keyword_filtered_fasta_content = raw_fasta_content
    if not skip_keyword_filter:
        log.info(
            f"Gene '{gene_symbol}': Applying keyword filter with keyword '{gene_input.query_keyword}'..."
        )
        keyword_filtered_fasta_content = filter_fasta_by_keyword(
            raw_fasta_content, gene_input.query_keyword, gene_symbol_for_log=gene_symbol
        )
        if not keyword_filtered_fasta_content and raw_fasta_content:
            log.warning(
                f"Gene '{gene_symbol}': Keyword filtering resulted in zero sequences."
            )

        if output_dir_individual and keyword_filtered_fasta_content:
            try:
                kw_filtered_path = (
                    output_dir_individual / f"{gene_symbol}_1_keyword_filtered.fasta"
                )
                with open(kw_filtered_path, "w") as f_kw:
                    f_kw.write(keyword_filtered_fasta_content)
                log.debug(
                    f"Saved keyword-filtered FASTA for {gene_symbol} to {kw_filtered_path}"
                )
            except Exception as e_save_kw:
                log.error(
                    f"Error saving keyword-filtered FASTA for {gene_symbol}: {e_save_kw}"
                )
    else:
        log.info(f"Gene '{gene_symbol}': Skipping keyword filtering.")
        if output_dir_individual and keyword_filtered_fasta_content:
            try:
                kw_filtered_path = (
                    output_dir_individual / f"{gene_symbol}_1_keyword_filtered.fasta"
                )
                with open(kw_filtered_path, "w") as f_kw:
                    f_kw.write(keyword_filtered_fasta_content)
                log.debug(
                    f"Saved (unfiltered by keyword) FASTA for {gene_symbol} as {kw_filtered_path}"
                )
            except Exception as e_save_kw_skipped:
                log.error(
                    f"Error saving (unfiltered by keyword) FASTA for {gene_symbol}: {e_save_kw_skipped}"
                )

    if not keyword_filtered_fasta_content.strip():
        log.info(
            f"No FASTA content to process for gene {gene_symbol} after optional keyword filter."
        )
        stats["status"] = "No content post-keyword-filter"
        stats["headers_encountered_in_processor"] = 0
        stats["final_sequences_kept_by_processor"] = 0
        stats["final_sequences_kept"] = 0
        return gene_symbol, [], stats

    fasta_stream = StringIO(keyword_filtered_fasta_content)
    processed_proteins_before_cdhit, processor_stats = process_fasta_stream(
        fasta_stream, gene_input.gene_symbol, max_dist
    )
    stats.update(processor_stats)

    final_processed_proteins = apply_cdhit_processing(
        processed_proteins_before_cdhit,
        gene_symbol,
        enable_cdhit,
        fixed_cdhit_threshold,
        output_dir_individual,
    )

    stats["sequences_before_cdhit"] = len(processed_proteins_before_cdhit)
    stats["sequences_after_cdhit"] = len(final_processed_proteins)
    stats["cdhit_applied"] = enable_cdhit
    if enable_cdhit:
        stats["final_sequences_kept"] = len(final_processed_proteins)
    else:
        stats["final_sequences_kept"] = stats.get(
            "final_sequences_kept", len(processed_proteins_before_cdhit)
        )

    stats["status"] = "Processed by processor"
    if enable_cdhit:
        stats["status"] += " and CD-HIT"

    if not final_processed_proteins:
        log.info(
            f"No proteins kept after all processing (including CD-HIT if applied) for gene {gene_symbol}."
        )
    else:
        log.info(
            f"Gene {gene_symbol}: All processing complete. Kept {len(final_processed_proteins)} proteins."
        )

        if output_dir_individual:
            try:
                final_fasta_short = (
                    output_dir_individual / f"{gene_symbol}_2_final_short.fasta"
                )
                final_fasta_full = (
                    output_dir_individual / f"{gene_symbol}_2_final_full.fasta"
                )
                final_csv = output_dir_individual / f"{gene_symbol}_2_final_meta.csv"

                write_processed_proteins_to_fasta(
                    final_processed_proteins,
                    str(final_fasta_short),
                    use_full_header=False,
                )
                write_processed_proteins_to_fasta(
                    final_processed_proteins,
                    str(final_fasta_full),
                    use_full_header=True,
                )
                write_processed_proteins_to_csv(
                    final_processed_proteins, str(final_csv)
                )
                log.debug(
                    f"Final individual files for {gene_symbol} saved to {output_dir_individual}"
                )
            except Exception as e:
                log.error(
                    f"Error writing final individual files for {gene_symbol}: {e}"
                )

    return gene_symbol, final_processed_proteins, stats


def solo_redundancy_filtering_workflow(args: argparse.Namespace):
    log.info("Starting solo redundancy filtering workflow...")
    if not args.input_fasta_srf or not args.input_csv_srf:
        log.error(
            "--input-fasta-srf and --input-csv-srf are required for --solo-redundancy-filtering mode."
        )
        return 1

    input_fasta_path = Path(args.input_fasta_srf)
    input_csv_path = Path(args.input_csv_srf)
    output_main_dir = ensure_output_dir(args.output_dir)
    output_dir_individual_srf = ensure_output_dir(
        output_main_dir / OUTPUT_SUBDIR_INDIVIDUAL
    )

    if not input_fasta_path.is_file():
        log.error(f"Input FASTA file not found: {input_fasta_path}")
        return 1
    if not input_csv_path.is_file():
        log.error(f"Input CSV file not found: {input_csv_path}")
        return 1

    cdhit_is_enabled_by_args_srf = (
        args.enable_cdhit or args.cdhit_fixed_threshold is not None
    )
    if not cdhit_is_enabled_by_args_srf:
        log.error(
            "For --solo-redundancy-filtering, CD-HIT must be enabled via --enable-cdhit or --cdhit-fixed-threshold."
        )
        return 1

    if not is_cdhit_installed():
        log.error(
            "CD-HIT command not found. Please install CD-HIT and ensure it's in your PATH to use -srf with CD-HIT."
        )
        return 1

    all_proteins_from_fasta: Dict[str, ProcessedProtein] = {}
    try:
        for record in SeqIO.parse(input_fasta_path, "fasta"):
            all_proteins_from_fasta[record.id] = ProcessedProtein(
                accession=record.id,
                gene_symbol_input="unknown_gene_srf",
                identifier_from_header=(
                    record.description.split(None, 1)[0]
                    if record.description
                    else record.id
                ),
                sequence=str(record.seq),
                full_header=record.description,
            )
        log.info(
            f"Read {len(all_proteins_from_fasta)} sequences from {input_fasta_path}"
        )
    except Exception as e:
        log.error(f"Error reading input FASTA {input_fasta_path}: {e}")
        return 1

    proteins_by_gene_group: Dict[str, List[ProcessedProtein]] = defaultdict(list)
    try:
        with open(input_csv_path, "r", newline="") as csvfile:
            reader = csv.DictReader(csvfile)
            if "identifier" not in reader.fieldnames or "gene" not in reader.fieldnames:
                log.error(
                    "Input CSV must contain 'identifier' and 'gene' columns for -srf mode."
                )
                return 1
            for row in reader:
                accession = row["identifier"]
                gene_from_csv = row["gene"]
                gene_group_key = gene_from_csv.lower()
                species_info = row.get("species", accession)

                if accession in all_proteins_from_fasta:
                    protein_obj = all_proteins_from_fasta[accession]
                    protein_obj.gene_symbol_input = gene_group_key
                    protein_obj.identifier_from_header = species_info
                    proteins_by_gene_group[gene_group_key].append(protein_obj)
        log.info(
            f"Processed CSV {input_csv_path}, found {len(proteins_by_gene_group)} gene groups (case-insensitive)."
        )
    except Exception as e:
        log.error(f"Error reading or processing input CSV {input_csv_path}: {e}")
        return 1

    final_proteins_to_combine: List[ProcessedProtein] = []

    num_srf_workers = args.cdhit_group_workers
    log.info(
        f"Processing {len(proteins_by_gene_group)} gene groups for CD-HIT using up to {num_srf_workers} worker(s) in SRF mode..."
    )

    with ThreadPoolExecutor(max_workers=num_srf_workers) as executor:
        future_to_gene_group: Dict[Any, str] = {}
        for gene_group_key, proteins_in_group in proteins_by_gene_group.items():
            if not proteins_in_group:
                continue

            srf_gene_group_output_dir = (
                output_dir_individual_srf / gene_group_key
                if args.save_individual_files
                else None
            )
            if srf_gene_group_output_dir:
                srf_gene_group_output_dir.mkdir(parents=True, exist_ok=True)

            future = executor.submit(
                apply_cdhit_processing,
                proteins_in_group,
                gene_group_key,
                True,  # CD-HIT is enabled for SRF
                args.cdhit_fixed_threshold,
                srf_gene_group_output_dir,
            )
            future_to_gene_group[future] = gene_group_key

        for future in tqdm(
            as_completed(future_to_gene_group),
            total=len(future_to_gene_group),
            desc="SRF CD-HIT Processing",
        ):
            gene_group_key_processed = future_to_gene_group[future]
            try:
                cdhit_filtered_group = future.result()
                final_proteins_to_combine.extend(cdhit_filtered_group)

                # Save individual final files for this gene group if requested
                if args.save_individual_files and cdhit_filtered_group:
                    srf_gene_group_output_dir = (
                        output_dir_individual_srf / gene_group_key_processed
                    )  # Reconstruct path
                    if srf_gene_group_output_dir:  # Ensure it was created
                        srf_final_short = (
                            srf_gene_group_output_dir
                            / f"{gene_group_key_processed}_srf_final_short.fasta"
                        )
                        srf_final_full = (
                            srf_gene_group_output_dir
                            / f"{gene_group_key_processed}_srf_final_full.fasta"
                        )
                        srf_final_csv = (
                            srf_gene_group_output_dir
                            / f"{gene_group_key_processed}_srf_final_meta.csv"
                        )
                        write_processed_proteins_to_fasta(
                            cdhit_filtered_group,
                            str(srf_final_short),
                            use_full_header=False,
                        )
                        write_processed_proteins_to_fasta(
                            cdhit_filtered_group,
                            str(srf_final_full),
                            use_full_header=True,
                        )
                        write_processed_proteins_to_csv(
                            cdhit_filtered_group, str(srf_final_csv)
                        )
            except Exception as exc:
                log.error(
                    f"Gene group '{gene_group_key_processed}' generated an exception during parallel CD-HIT processing: {exc}"
                )

    if not final_proteins_to_combine:
        log.info("No proteins remained after CD-HIT processing in -srf mode.")
        return 0

    dedup_by_accession_list = deduplicate_processed_proteins(
        final_proteins_to_combine, dedup_key="accession"
    )
    dedup_by_accession_list.sort(key=lambda p: p.accession)

    dedup_by_full_header_list = deduplicate_processed_proteins(
        final_proteins_to_combine, dedup_key="full_header"
    )
    dedup_by_full_header_list.sort(key=lambda p: (p.full_header, p.accession))

    input_file_stem = input_fasta_path.stem
    srf_suffix = "_srf_cdhit_combined"

    final_fasta_short_path = (
        output_main_dir / f"{input_file_stem}{srf_suffix}_short.fasta"
    )
    final_fasta_full_path = (
        output_main_dir / f"{input_file_stem}{srf_suffix}_full.fasta"
    )
    final_csv_path = output_main_dir / f"{input_file_stem}{srf_suffix}_meta.csv"

    log.info(
        f"Writing SRF CD-HIT combined short FASTA ({len(dedup_by_accession_list)} proteins) to {final_fasta_short_path}"
    )
    write_processed_proteins_to_fasta(
        dedup_by_accession_list, str(final_fasta_short_path), use_full_header=False
    )

    log.info(
        f"Writing SRF CD-HIT combined full header FASTA ({len(dedup_by_full_header_list)} proteins) to {final_fasta_full_path}"
    )
    write_processed_proteins_to_fasta(
        dedup_by_full_header_list, str(final_fasta_full_path), use_full_header=True
    )

    log.info(
        f"Writing SRF CD-HIT combined metadata CSV ({len(dedup_by_accession_list)} proteins) to {final_csv_path}"
    )
    write_processed_proteins_to_csv(dedup_by_accession_list, str(final_csv_path))

    log.info("Solo redundancy filtering workflow completed.")
    return 0


def main_workflow(args: argparse.Namespace):
    start_time = time.time()

    log_level = logging.DEBUG if args.debug else logging.INFO
    setup_logging(level=log_level)
    log.info(
        f"protfetch version {__version__} starting with log level {logging.getLevelName(log_level)}."
    )

    cdhit_is_enabled_by_args = (
        args.enable_cdhit or args.cdhit_fixed_threshold is not None
    )
    if args.solo_redundancy_filtering:
        cdhit_is_enabled_by_args = True
        if not args.enable_cdhit and args.cdhit_fixed_threshold is None:
            log.error(
                "SRF mode selected but CD-HIT not enabled via --enable-cdhit or --cdhit-fixed-threshold."
            )
            return 1

    if cdhit_is_enabled_by_args and not is_cdhit_installed():
        log.error(
            "CD-HIT processing was requested but cd-hit executable was not found in PATH. "
            "Please install CD-HIT or remove CD-HIT related flags."
        )
        return 1

    if args.solo_redundancy_filtering:
        return solo_redundancy_filtering_workflow(args)

    if args.entrez_email == DEFAULT_ENTREZ_EMAIL:
        log.warning(
            f"Using default Entrez email: {DEFAULT_ENTREZ_EMAIL}. "
            "Please provide your own email using --entrez-email for reliable NCBI access."
        )
    configure_entrez(args.entrez_email, args.entrez_api_key)

    output_main_dir = ensure_output_dir(args.output_dir)
    output_dir_individual: Union[Path, None] = None
    if args.save_individual_files:
        output_dir_individual = ensure_output_dir(
            output_main_dir / OUTPUT_SUBDIR_INDIVIDUAL
        )

    if not args.input_gene_list_file:
        log.error(
            "Input gene list file is required when not in --solo-redundancy-filtering mode."
        )
        parser.print_help()
        return 1

    try:
        genes_to_process: List[GeneInput] = parse_gene_list_file(
            args.input_gene_list_file
        )
    except Exception:
        log.error(
            f"Failed to parse input gene list file: {args.input_gene_list_file}",
            exc_info=True,
        )
        return 1

    if not genes_to_process:
        log.info("No valid genes found in the input list. Exiting.")
        return 0

    log.info(
        f"Found {len(genes_to_process)} genes/queries to process from '{args.input_gene_list_file}'."
    )

    all_processed_proteins: List[ProcessedProtein] = []
    all_stats: List[Dict[str, Any]] = []

    num_workers = min(args.max_workers, len(genes_to_process))
    if num_workers < 1:
        num_workers = 1
    log.info(f"Processing genes using up to {num_workers} concurrent worker(s)...")

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        future_to_gene_input: Dict[Any, GeneInput] = {
            executor.submit(
                process_single_gene_task,
                gene_input,
                output_dir_individual,
                args.max_dist,
                args.timeout,
                args.retries,
                args.skip_keyword_filter,
                cdhit_is_enabled_by_args,
                args.cdhit_fixed_threshold,
            ): gene_input
            for gene_input in genes_to_process
        }

        for future in tqdm(
            as_completed(future_to_gene_input),
            total=len(genes_to_process),
            desc="Processing genes",
        ):
            gene_input_obj = future_to_gene_input[future]
            try:
                result = future.result()
                if result:
                    gene_sym, proteins_list, gene_stats_dict = result
                    if proteins_list:
                        all_processed_proteins.extend(proteins_list)
                    if gene_stats_dict:
                        all_stats.append(gene_stats_dict)
                else:
                    log.error(
                        f"No result tuple returned for gene: {gene_input_obj.gene_symbol}"
                    )
                    all_stats.append(
                        {
                            "gene_symbol": gene_input_obj.gene_symbol,
                            "error": "Task failed unexpectedly",
                            "status": "Task error",
                        }
                    )
            except Exception as e:
                log.error(
                    f"Error processing gene {gene_input_obj.gene_symbol} in worker future: {e}",
                    exc_info=True,
                )
                all_stats.append(
                    {
                        "gene_symbol": gene_input_obj.gene_symbol,
                        "error": str(e),
                        "status": "Worker exception",
                    }
                )

    log.info(f"Finished processing all {len(genes_to_process)} gene inputs.")

    dedup_by_accession_list: List[ProcessedProtein] = []
    if not all_processed_proteins:
        log.info(
            "No proteins were successfully processed from any gene or after CD-HIT. No combined files will be created."
        )
    else:
        log.info(
            f"Combining final results from {len(all_processed_proteins)} proteins before final deduplication..."
        )

        dedup_by_accession_list = deduplicate_processed_proteins(
            all_processed_proteins, dedup_key="accession"
        )
        dedup_by_accession_list.sort(key=lambda p: p.accession)

        dedup_by_full_header_list = deduplicate_processed_proteins(
            all_processed_proteins, dedup_key="full_header"
        )
        dedup_by_full_header_list.sort(key=lambda p: (p.full_header, p.accession))

        base_name_for_combined = "protfetch_all_genes"
        if args.input_gene_list_file:
            base_name_for_combined = Path(args.input_gene_list_file).stem

        combined_fasta_short_path = (
            output_main_dir / f"{base_name_for_combined}{COMBINED_FASTA_SHORT_SUFFIX}"
        )
        combined_fasta_full_path = (
            output_main_dir / f"{base_name_for_combined}{COMBINED_FASTA_FULL_SUFFIX}"
        )
        combined_csv_path = (
            output_main_dir / f"{base_name_for_combined}{COMBINED_CSV_SUFFIX}"
        )

        log.info(
            f"Writing combined short FASTA ({len(dedup_by_accession_list)} proteins) to {combined_fasta_short_path}"
        )
        write_processed_proteins_to_fasta(
            dedup_by_accession_list,
            str(combined_fasta_short_path),
            use_full_header=False,
        )

        log.info(
            f"Writing combined full header FASTA ({len(dedup_by_full_header_list)} proteins) to {combined_fasta_full_path}"
        )
        write_processed_proteins_to_fasta(
            dedup_by_full_header_list,
            str(combined_fasta_full_path),
            use_full_header=True,
        )

        log.info(
            f"Writing combined metadata CSV ({len(dedup_by_accession_list)} proteins) to {combined_csv_path}"
        )
        write_processed_proteins_to_csv(dedup_by_accession_list, str(combined_csv_path))

    log.info("--- protfetch Run Summary ---")
    successful_genes = 0
    failed_genes = 0
    genes_with_no_final_proteins = 0

    for gene_stat_summary in all_stats:
        gs = gene_stat_summary.get("gene_symbol", "Unknown Gene")
        status = gene_stat_summary.get("status", "Unknown status")

        if (
            "error" in gene_stat_summary
            or status == "Fetch failed"
            or status == "Worker exception"
        ):
            failed_genes += 1
            err_msg = gene_stat_summary.get("error", "Unknown error")
            log.info(f"  Gene {gs}: Failed ({status} - {err_msg})")
        elif status == "No content post-keyword-filter":
            log.info(f"  Gene {gs}: Processed, but no content after keyword filtering.")
            successful_genes += 1
            genes_with_no_final_proteins += 1
        elif "Processed by processor" in status:
            final_kept_overall = gene_stat_summary.get("final_sequences_kept", 0)
            initial_for_processor = gene_stat_summary.get("headers_encountered", "N/A")
            cdhit_info = ""
            if gene_stat_summary.get("cdhit_applied"):
                before_cdhit = gene_stat_summary.get("sequences_before_cdhit", "N/A")
                after_cdhit = gene_stat_summary.get("sequences_after_cdhit", "N/A")
                cdhit_info = f" (CD-HIT: {before_cdhit} -> {after_cdhit})"

            log.info(
                f"  Gene {gs}: Processed. Processor stage started with {initial_for_processor} records, final kept: {final_kept_overall}{cdhit_info}."
            )
            successful_genes += 1
            if final_kept_overall == 0:
                genes_with_no_final_proteins += 1
        else:
            log.warning(f"  Gene {gs}: Unknown processing status - {status}")
            failed_genes += 1

    log.info(
        f"Total genes attempted: {len(genes_to_process) if genes_to_process else 0}"
    )
    log.info(f"  Successfully initiated processing for: {successful_genes} gene(s)")
    log.info(f"  Failed during fetch or early processing: {failed_genes} gene(s)")
    log.info(
        f"  Genes with zero proteins after all filtering: {genes_with_no_final_proteins} (subset of successfully processed)"
    )
    log.info(
        f"Total unique proteins in combined outputs (dedup by accession): {len(dedup_by_accession_list)}"
    )

    end_time = time.time()
    log.info(f"Total execution time: {end_time - start_time:.2f} seconds.")
    log.info("--- End of Summary ---")

    if failed_genes > 0 and successful_genes == 0:
        return 1
    return 0


parser = argparse.ArgumentParser(
    description="protfetch: Fetch and process protein FASTA/metadata from NCBI. Optionally use CD-HIT for redundancy reduction.",
    formatter_class=argparse.RawTextHelpFormatter,
)


def cli_entry():
    if not any(action.dest == "input_gene_list_file" for action in parser._actions):
        fetching_group = parser.add_argument_group("Standard Fetching Mode (default)")
        fetching_group.add_argument(
            "input_gene_list_file",
            type=str,
            nargs="?",
            default=None,
            help="Path to the input file containing gene list (required if not in -srf mode). \n"
            "Formats supported:\n"
            "1. 'Protein Name | GENE_SYMBOL' (per line)\n"
            "2. 'GENE_SYMBOL' (per line)",
        )

        srf_group = parser.add_argument_group("Solo Redundancy Filtering Mode")
        srf_group.add_argument(
            "--solo-redundancy-filtering",
            "-srf",
            action="store_true",
            help="Enable solo redundancy filtering mode. Operates on existing FASTA and CSV files using CD-HIT. "
            "Requires --input-fasta-srf, --input-csv-srf, and CD-HIT to be enabled via --enable-cdhit or --cdhit-fixed-threshold.",
        )
        srf_group.add_argument(
            "--input-fasta-srf",
            type=str,
            help="Path to input FASTA file for -srf mode.",
        )
        srf_group.add_argument(
            "--input-csv-srf",
            type=str,
            help="Path to input CSV metadata file for -srf mode (must contain 'identifier' and 'gene' columns).",
        )
        srf_group.add_argument(  # New argument for SRF mode workers
            "--cdhit-group-workers",
            type=int,
            default=6,  # Default to 6 as requested
            help="Maximum number of concurrent workers for processing gene groups with CD-HIT in -srf mode. Default: 6.",
        )

        common_group = parser.add_argument_group("Common Options")
        common_group.add_argument(
            "-o",
            "--output-dir",
            type=str,
            default="protfetch_results",
            help="Directory to save output files (default: protfetch_results).",
        )
        common_group.add_argument(
            "--entrez-email",
            type=str,
            default=DEFAULT_ENTREZ_EMAIL,
            help=f"Email address for NCBI Entrez. Default: {DEFAULT_ENTREZ_EMAIL}",
        )
        common_group.add_argument(
            "--entrez-api-key",
            type=str,
            default=None,
            help="NCBI API key for higher Entrez request rates (optional).",
        )
        common_group.add_argument(
            "--max-dist",
            type=int,
            default=DEFAULT_MAX_DIST,
            help=f"Max Levenshtein distance for internal filtering (default: 4; 0 to disable this filter if using CD-HIT primarily).",
        )

        cdhit_options_group = parser.add_argument_group(
            "CD-HIT Options (require CD-HIT installed)"
        )
        cdhit_options_group.add_argument(  # Renamed group to avoid conflict if 'cdhit_group' was used before
            "--enable-cdhit",
            "-ec",
            action="store_true",
            help="Enable CD-HIT based redundancy filtering. Uses automatic thresholds per gene/class unless --cdhit-fixed-threshold is set.",
        )
        cdhit_options_group.add_argument(
            "--cdhit-fixed-threshold",
            type=float,
            default=None,
            help="Specify a fixed CD-HIT identity threshold (e.g., 0.9, 0.95). "
            "If set, this overrides the automatic threshold logic. Implies CD-HIT is enabled.",
        )

        # --max-workers is for fetching, --cdhit-group-workers is for SRF CD-HIT runs
        fetching_options_group = parser.add_argument_group("Fetching Options")
        fetching_options_group.add_argument(
            "--timeout",
            type=int,
            default=DEFAULT_REQUEST_TIMEOUT,
            help=f"Timeout for NCBI Entrez requests in seconds. Default: {DEFAULT_REQUEST_TIMEOUT}.",
        )
        fetching_options_group.add_argument(
            "--retries",
            type=int,
            default=DEFAULT_REQUEST_RETRIES,
            help=f"Number of retries for failed NCBI Entrez requests. Default: {DEFAULT_REQUEST_RETRIES}.",
        )
        fetching_options_group.add_argument(
            "--max-workers",  # This is for the standard fetching mode
            type=int,
            default=DEFAULT_MAX_WORKERS,
            help=f"Maximum number of concurrent workers for fetching data (standard mode). Default: {DEFAULT_MAX_WORKERS}.",
        )

        output_options_group = parser.add_argument_group("Output Options")
        output_options_group.add_argument(
            "--save-individual-files",
            action="store_true",
            help="Save processed files for each gene individually. Also saves intermediate FASTA files (raw NCBI, keyword-filtered, pre-CD-HIT).",
        )
        output_options_group.add_argument(
            "--skip-keyword-filter",
            action="store_true",
            help="Skip filtering fetched FASTA by keyword.",
        )

        general_group = parser.add_argument_group("General Options")
        general_group.add_argument(
            "--debug", action="store_true", help="Enable detailed debug logging."
        )
        general_group.add_argument(
            "-v", "--version", action="version", version=f"%(prog)s {__version__}"
        )

    args = parser.parse_args()

    if args.solo_redundancy_filtering:
        if not args.input_fasta_srf or not args.input_csv_srf:
            parser.error(
                "--input-fasta-srf and --input-csv-srf are required when using --solo-redundancy-filtering."
            )
        if not args.enable_cdhit and args.cdhit_fixed_threshold is None:
            parser.error(
                "CD-HIT must be enabled (via --enable-cdhit or by setting --cdhit-fixed-threshold) when using --solo-redundancy-filtering."
            )
        if args.cdhit_group_workers < 1:
            parser.error("--cdhit-group-workers must be at least 1.")
    elif not args.input_gene_list_file:
        parser.error(
            "input_gene_list_file is required when not using --solo-redundancy-filtering mode."
        )

    if args.cdhit_fixed_threshold is not None:
        if not (0.4 <= args.cdhit_fixed_threshold <= 1.0):
            parser.error("--cdhit-fixed-threshold must be between 0.4 and 1.0.")

    if args.max_workers < 1:  # This refers to the fetching workers
        parser.error("--max-workers must be at least 1.")
    if args.max_dist < 0:
        parser.error("--max-dist cannot be negative.")

    return main_workflow(args)


if __name__ == "__main__":
    cli_entry()
