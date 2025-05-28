# protfetch/fetcher.py
import time
from Bio import Entrez, SeqIO
from io import StringIO
import requests  # For retryable sessions

from .utils import log, GeneInput
from .processor import (
    ProcessedProtein,
)  # For type hinting if needed, though fetcher primarily returns raw FASTA


# --- Entrez Configuration ---
def configure_entrez(email: str, api_key: str | None = None):
    """Configures Biopython Entrez."""
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    log.info(
        f"Entrez configured with email: {email}" + (" and API key." if api_key else ".")
    )


# --- NCBI Interaction with Retries ---
def _entrez_retry_call(entrez_func, *args, retries=3, delay=5, **kwargs):
    """
    Wrapper for Entrez calls with retry logic for common transient errors.
    """
    for attempt in range(retries):
        try:
            handle = entrez_func(*args, **kwargs)
            return handle
        except (
            requests.exceptions.ConnectionError,
            requests.exceptions.Timeout,
        ) as e:  # HTTP errors often wrapped by Biopython
            log.warning(
                f"Entrez call failed (Attempt {attempt + 1}/{retries}): {e}. Retrying in {delay}s..."
            )
            if attempt + 1 == retries:
                log.error(f"Entrez call failed after {retries} attempts: {e}")
                raise
            time.sleep(delay)
        except Exception as e:  # Catch other Entrez/BioPython errors
            # Check for specific NCBI related HTTP errors if possible, e.g. from urllib.error.HTTPError
            # Biopython might raise generic IOErrors or RuntimeErrors for some NCBI issues.
            if "HTTP Error" in str(e) or "NCBI" in str(e):
                log.warning(
                    f"Entrez call resulted in an error (Attempt {attempt + 1}/{retries}): {e}. Retrying in {delay}s..."
                )
                if attempt + 1 == retries:
                    log.error(
                        f"Entrez call failed after {retries} attempts with NCBI error: {e}"
                    )
                    raise
                time.sleep(delay)
            else:  # Non-retryable error
                log.error(f"Unexpected error during Entrez call: {e}")
                raise
    return None  # Should not be reached if retries exhausted and error raised


# --- Core Fetching Logic ---
def fetch_protein_fasta_for_gene(
    gene_input: GeneInput, timeout: int, retries: int
) -> str | None:
    """
    Fetches protein FASTA sequences from NCBI for a given gene symbol.
    This replaces the esearch, elink, efetch part of the shell script.

    Args:
        gene_input: GeneInput object containing gene symbol and keyword.
        timeout: Timeout for Entrez requests.
        retries: Number of retries for Entrez requests.

    Returns:
        A string containing FASTA sequences, or None if an error occurs or no data found.
    """
    gene_symbol = gene_input.gene_symbol
    log.info(f"Fetching data for gene symbol: '{gene_symbol}'")

    try:
        # 1. Get Gene UID from gene symbol
        #    (Original: esearch -db gene -query "${gene_symbol}[symbol]" | efetch -format uid)
        log.debug(f"Gene '{gene_symbol}': [1/3] Fetching Gene UIDs...")
        search_term = f"{gene_symbol}[Gene Name] AND (animals[Filter] OR human[Filter]) AND alive[Prop]"  # More specific query
        # Consider organism filter if needed, e.g., "Homo sapiens[Organism]"
        # The original script implies fetching orthologs, so a broad search might be intended.
        # Let's try to be a bit more specific to reduce noise, but this might need adjustment.
        # For true orthologs, one might use a different database or tool, but elink from gene to protein is a common proxy.

        handle_search = _entrez_retry_call(
            Entrez.esearch,
            db="gene",
            term=search_term,
            retmax="10",  # Usually one primary gene, but handle multiple if symbol is ambiguous
            retries=retries,
            delay=5,  # Standard delay between Entrez queries
        )
        if not handle_search:
            return None

        record_search = Entrez.read(handle_search)
        handle_search.close()
        gene_ids = record_search["IdList"]

        if not gene_ids:
            log.warning(
                f"Gene '{gene_symbol}': No Gene UIDs found for symbol '{gene_symbol}'."
            )
            return None
        log.debug(f"Gene '{gene_symbol}': Found Gene UIDs: {gene_ids}")

        # 2. Get linked Protein UIDs from Gene UIDs
        #    (Original: elink -db gene -target protein -input "$uid_file")
        log.debug(
            f"Gene '{gene_symbol}': [2/3] Fetching linked Protein UIDs from Gene UIDs: {gene_ids}..."
        )
        time.sleep(1)  # NCBI rate limiting: 3 req/sec without API key, 10 req/sec with.

        protein_ids_all = []
        # Batch elink if many gene_ids (though typically few for a specific symbol)
        # Max 200 IDs per epost/elink URL, but Entrez.elink handles list of IDs directly.
        handle_elink = _entrez_retry_call(
            Entrez.elink,
            dbfrom="gene",
            db="protein",
            id=gene_ids,
            linkname="gene_protein_refseq",  # Link to RefSeq proteins
            retries=retries,
            delay=5,
        )
        if not handle_elink:
            return None

        record_elink = Entrez.read(handle_elink)
        handle_elink.close()

        for linkset in record_elink:
            if linkset.get("LinkSetDb"):
                for link in linkset["LinkSetDb"][0]["Link"]:
                    protein_ids_all.append(link["Id"])

        if not protein_ids_all:
            log.warning(
                f"Gene '{gene_symbol}': No linked Protein UIDs found for Gene UIDs {gene_ids}."
            )
            # Try a broader link if RefSeq yields nothing
            log.debug(
                f"Gene '{gene_symbol}': Retrying elink with broader scope (gene_protein)."
            )
            time.sleep(1)
            handle_elink_broad = _entrez_retry_call(
                Entrez.elink,
                dbfrom="gene",
                db="protein",
                id=gene_ids,
                retries=retries,
                delay=5,
            )
            if not handle_elink_broad:
                return None
            record_elink_broad = Entrez.read(handle_elink_broad)
            handle_elink_broad.close()
            for linkset in record_elink_broad:
                if linkset.get("LinkSetDb"):
                    for link in linkset["LinkSetDb"][0]["Link"]:
                        protein_ids_all.append(link["Id"])
            if not protein_ids_all:
                log.warning(
                    f"Gene '{gene_symbol}': Still no linked Protein UIDs after broader search."
                )
                return None

        log.debug(
            f"Gene '{gene_symbol}': Found {len(protein_ids_all)} linked Protein UIDs (pre-filter). E.g., {protein_ids_all[:5]}"
        )

        # 3. Fetch Protein FASTA sequences for Protein UIDs
        #    (Original: efetch -db protein -format fasta < linked_protein_uids)
        log.debug(
            f"Gene '{gene_symbol}': [3/3] Fetching FASTA for {len(protein_ids_all)} Protein UIDs..."
        )
        time.sleep(1)

        # Efetch can handle multiple IDs. Batch if very large, but for one gene, usually manageable.
        # Max IDs for GET is around 200. POST is better for large lists. Biopython handles this.
        fasta_data_list = []
        # Fetch in batches to avoid overly large requests
        batch_size = 150
        for i in range(0, len(protein_ids_all), batch_size):
            batch_ids = protein_ids_all[i : i + batch_size]
            log.debug(
                f"Gene '{gene_symbol}': Fetching FASTA batch {i//batch_size + 1} for {len(batch_ids)} IDs."
            )
            handle_efetch = _entrez_retry_call(
                Entrez.efetch,
                db="protein",
                id=batch_ids,
                rettype="fasta",
                retmode="text",
                retries=retries,
                delay=5,
            )
            if not handle_efetch:
                continue  # Skip batch on error, or decide to fail all

            fasta_batch_data = handle_efetch.read()
            handle_efetch.close()
            fasta_data_list.append(fasta_batch_data)
            if i + batch_size < len(protein_ids_all):  # If not the last batch
                time.sleep(1)  # Pause between batches

        raw_fasta_content = "".join(fasta_data_list)

        if not raw_fasta_content.strip():
            log.warning(
                f"Gene '{gene_symbol}': No FASTA data returned for Protein UIDs."
            )
            return None

        log.info(
            f"Gene '{gene_symbol}': Successfully fetched raw FASTA data ({len(raw_fasta_content)} bytes)."
        )
        return raw_fasta_content

    except Exception as e:
        log.error(f"Gene '{gene_symbol}': Error during NCBI fetch: {e}", exc_info=True)
        return None


def filter_fasta_by_keyword(
    fasta_content_string: str, keyword: str, gene_symbol_for_log: str = ""
) -> str:
    """
    Filters FASTA records based on a keyword found in the header.
    Replaces the awk keyword filtering from the shell script.
    """
    if not keyword:  # If no keyword, return all
        return fasta_content_string

    log.debug(
        f"Filtering FASTA content (gene: {gene_symbol_for_log}) with keyword: '{keyword}'"
    )
    filtered_records = []
    num_total_records = 0
    try:
        for record in SeqIO.parse(StringIO(fasta_content_string), "fasta"):
            num_total_records += 1
            # record.description is the full header line after '>'
            if keyword.lower() in record.description.lower():
                filtered_records.append(record)

        if not filtered_records:
            log.warning(
                f"Keyword '{keyword}' not found in any headers for gene '{gene_symbol_for_log}'. "
                f"Original FASTA had {num_total_records} records."
            )
            # Return empty string, or original if no filtering should occur?
            # Original AWK script would produce empty file.
            return ""

        output_fasta_io = StringIO()
        SeqIO.write(filtered_records, output_fasta_io, "fasta")
        filtered_fasta_str = output_fasta_io.getvalue()
        log.debug(
            f"Keyword filtering for gene '{gene_symbol_for_log}': {len(filtered_records)}/{num_total_records} records kept."
        )
        return filtered_fasta_str

    except Exception as e:
        log.error(
            f"Error during keyword filtering for gene '{gene_symbol_for_log}': {e}",
            exc_info=True,
        )
        return ""  # Return empty on error to mimic awk behavior of producing empty on failure
