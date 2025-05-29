# protfetch/main.py
import argparse
import logging
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO
from pathlib import Path
from typing import Dict, List, Tuple, Union

from tqdm import tqdm

from . import __version__, utils
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


def process_single_gene_task(
    gene_input: GeneInput,
    output_dir_individual: Union[Path, None],  # Corrected type hint for Path | None
    max_dist: int,
    entrez_timeout: int,
    entrez_retries: int,
    skip_keyword_filter: bool,
) -> Union[Tuple[str, List[ProcessedProtein], Dict], None]:  # Corrected type hint
    gene_symbol = gene_input.gene_symbol
    log.info(
        f"Starting processing for gene: {gene_symbol} (Keyword: '{gene_input.query_keyword}')"
    )

    raw_fasta_content = fetch_protein_fasta_for_gene(
        gene_input, entrez_timeout, entrez_retries
    )
    if not raw_fasta_content:
        log.error(f"Failed to fetch FASTA data for gene {gene_symbol}.")
        return gene_symbol, [], {"error": "Fetch failed", "status": "Fetch failed"}

    keyword_filtered_fasta_content = raw_fasta_content  # Initialize with raw content
    if not skip_keyword_filter:
        log.info(
            f"Gene '{gene_symbol}': Applying keyword filter with keyword '{gene_input.query_keyword}'..."
        )
        keyword_filtered_fasta_content = filter_fasta_by_keyword(
            raw_fasta_content, gene_input.query_keyword, gene_symbol_for_log=gene_symbol
        )
        if (
            not keyword_filtered_fasta_content and raw_fasta_content
        ):  # Check if filtering resulted in empty but started with content
            log.warning(
                f"Gene '{gene_symbol}': Keyword filtering resulted in zero sequences."
            )
    else:
        log.info(f"Gene '{gene_symbol}': Skipping keyword filtering.")

    if not keyword_filtered_fasta_content.strip():
        log.info(
            f"No FASTA content to process for gene {gene_symbol} after optional keyword filter."
        )
        # Return stats indicating no content for processor step
        return (
            gene_symbol,
            [],
            {
                "status": "No content post-keyword-filter",
                "headers_encountered_in_processor": 0,
                "final_sequences_kept_by_processor": 0,
            },
        )

    fasta_stream = StringIO(keyword_filtered_fasta_content)
    processed_proteins, stats = process_fasta_stream(
        fasta_stream, gene_symbol, max_dist
    )

    stats["gene_symbol"] = gene_symbol
    stats["status"] = "Processed by processor"

    if not processed_proteins:
        log.info(
            f"No proteins kept after sequence processing filters for gene {gene_symbol}."
        )
    else:
        # This log now reflects the outcome of the processor.py filters
        log.info(
            f"Gene {gene_symbol}: Sequence processing filters complete. Kept {len(processed_proteins)} proteins from {stats.get('headers_encountered', 'N/A')} records fed to processor."
        )

        if output_dir_individual:
            try:
                individual_fasta_short = (
                    output_dir_individual / f"{gene_symbol}_filtered_short.fasta"
                )
                individual_fasta_full = (
                    output_dir_individual / f"{gene_symbol}_filtered_full.fasta"
                )
                individual_csv = (
                    output_dir_individual / f"{gene_symbol}_filtered_meta.csv"
                )

                write_processed_proteins_to_fasta(
                    processed_proteins,
                    str(individual_fasta_short),
                    use_full_header=False,
                )
                write_processed_proteins_to_fasta(
                    processed_proteins, str(individual_fasta_full), use_full_header=True
                )
                write_processed_proteins_to_csv(processed_proteins, str(individual_csv))
                log.debug(
                    f"Individual files for {gene_symbol} saved to {output_dir_individual}"
                )
            except Exception as e:
                log.error(f"Error writing individual files for {gene_symbol}: {e}")

    return gene_symbol, processed_proteins, stats


def main_workflow(args):
    start_time = time.time()

    log_level = logging.DEBUG if args.debug else logging.INFO
    setup_logging(level=log_level)
    log.info(
        f"protfetch version {__version__} starting with log level {logging.getLevelName(log_level)}."
    )

    if args.entrez_email == DEFAULT_ENTREZ_EMAIL:
        log.warning(
            f"Using default Entrez email: {DEFAULT_ENTREZ_EMAIL}. "
            "Please provide your own email using --entrez-email for reliable NCBI access."
        )
    configure_entrez(args.entrez_email, args.entrez_api_key)

    output_main_dir = ensure_output_dir(args.output_dir)
    output_dir_individual = None
    if args.save_individual_files:
        output_dir_individual = ensure_output_dir(
            output_main_dir / OUTPUT_SUBDIR_INDIVIDUAL
        )

    try:
        genes_to_process = parse_gene_list_file(args.input_gene_list_file)
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
    all_stats: List[Dict] = []  # Corrected type hint

    num_workers = min(args.max_workers, len(genes_to_process))
    if num_workers < 1:
        num_workers = 1  # Ensure at least 1 worker
    log.info(f"Processing genes using up to {num_workers} concurrent worker(s)...")

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        future_to_gene_input = {
            executor.submit(
                process_single_gene_task,
                gene_input,
                output_dir_individual,
                args.max_dist,
                args.timeout,
                args.retries,
                args.skip_keyword_filter,
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
                if result:  # result is (gene_symbol, proteins_list, stats_dict)
                    gene_sym, proteins_list, gene_stats_dict = result
                    if proteins_list:
                        all_processed_proteins.extend(proteins_list)
                    if gene_stats_dict:
                        all_stats.append(gene_stats_dict)
                else:  # Should not happen if process_single_gene_task always returns a tuple
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
            "No proteins were successfully processed from any gene. No combined files will be created."
        )
    else:
        log.info(
            f"Combining results from {len(all_processed_proteins)} initially collected proteins before final deduplication..."
        )

        dedup_by_accession_list = deduplicate_processed_proteins(
            all_processed_proteins, dedup_key="accession"
        )
        dedup_by_accession_list.sort(key=lambda p: p.accession)

        dedup_by_full_header_list = deduplicate_processed_proteins(
            all_processed_proteins, dedup_key="full_header"
        )
        dedup_by_full_header_list.sort(key=lambda p: (p.full_header, p.accession))

        input_file_stem = Path(args.input_gene_list_file).stem

        combined_fasta_short_path = (
            output_main_dir / f"{input_file_stem}{COMBINED_FASTA_SHORT_SUFFIX}"
        )
        combined_fasta_full_path = (
            output_main_dir / f"{input_file_stem}{COMBINED_FASTA_FULL_SUFFIX}"
        )
        combined_csv_path = output_main_dir / f"{input_file_stem}{COMBINED_CSV_SUFFIX}"

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
        elif status == "Processed by processor":
            final_kept_processor = gene_stat_summary.get("final_sequences_kept", 0)
            initial_for_processor = gene_stat_summary.get("headers_encountered", "N/A")
            log.info(
                f"  Gene {gs}: Processed. Processor stage started with {initial_for_processor} records, kept {final_kept_processor} sequences."
            )
            successful_genes += 1
            if final_kept_processor == 0:
                genes_with_no_final_proteins += 1
        else:  # Should not happen
            log.warning(f"  Gene {gs}: Unknown processing status - {status}")
            failed_genes += 1

    log.info(f"Total genes attempted: {len(genes_to_process)}")
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


def cli_entry():
    """
    Command Line Interface entry point.
    """
    parser = argparse.ArgumentParser(
        description="protfetch: Fetch and process protein FASTA/metadata from NCBI.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "input_gene_list_file",
        type=str,
        help="Path to the input file containing gene list. \n"
        "Formats supported:\n"
        "1. 'Protein Name | GENE_SYMBOL' (per line)\n"
        "2. 'GENE_SYMBOL' (per line)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        default="protfetch_results",
        help="Directory to save output files (default: protfetch_results).",
    )
    parser.add_argument(
        "--entrez-email",
        type=str,
        default=DEFAULT_ENTREZ_EMAIL,
        help=f"Email address for NCBI Entrez (required for reliable access). Default: {DEFAULT_ENTREZ_EMAIL}",
    )
    parser.add_argument(
        "--entrez-api-key",
        type=str,
        default=None,
        help="NCBI API key for higher Entrez request rates (optional).",
    )
    parser.add_argument(
        "--max-dist",
        type=int,
        default=DEFAULT_MAX_DIST,
        help=f"Max Levenshtein distance for filtering near-identical sequences (0 to disable). Default: {DEFAULT_MAX_DIST}.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=DEFAULT_REQUEST_TIMEOUT,
        help=f"Timeout for NCBI Entrez requests in seconds. Default: {DEFAULT_REQUEST_TIMEOUT}.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=DEFAULT_REQUEST_RETRIES,
        help=f"Number of retries for failed NCBI Entrez requests. Default: {DEFAULT_REQUEST_RETRIES}.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=DEFAULT_MAX_WORKERS,
        help=f"Maximum number of concurrent workers for fetching data. Default: {DEFAULT_MAX_WORKERS}.",
    )
    parser.add_argument(
        "--save-individual-files",
        action="store_true",
        help=f"Save processed FASTA and CSV files for each gene individually in a sub-directory ('{OUTPUT_SUBDIR_INDIVIDUAL}').",
    )
    parser.add_argument(
        "--skip-keyword-filter",
        action="store_true",
        help="Skip the step of filtering fetched FASTA sequences by the derived keyword (protein name or gene symbol).",
    )
    parser.add_argument(
        "--debug", action="store_true", help="Enable detailed debug logging."
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"%(prog)s {__version__}"
    )

    args = parser.parse_args()

    # Basic validation
    if args.max_workers < 1:
        parser.error("--max-workers must be at least 1.")
    if args.max_dist < 0:
        parser.error("--max-dist cannot be negative.")

    return main_workflow(args)


if __name__ == "__main__":
    # This allows running the main.py script directly for development/testing
    # e.g. python -m protfetch.main your_gene_list.txt --entrez-email your@email.com
    # However, the proper way to run after installation is `protfetch ...`
    cli_entry()
