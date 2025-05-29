import http.client
import time
from io import StringIO
from typing import Any, Callable, List, Union

import requests
from Bio import Entrez, SeqIO

from .utils import GeneInput, log


def configure_entrez(email: str, api_key: Union[str, None] = None):
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    log.info(
        f"Entrez configured with email: {email}" + (" and API key." if api_key else ".")
    )


def _entrez_retry_call(
    entrez_func: Callable[..., Any],
    *args: Any,
    retries: int = 3,
    delay: int = 5,
    **kwargs: Any,
) -> Any:
    for attempt in range(retries):  # This loop runs 'retries' times
        try:
            handle = entrez_func(*args, **kwargs)
            return handle
        except (
            requests.exceptions.ConnectionError,
            requests.exceptions.Timeout,
            requests.exceptions.RequestException,
            http.client.IncompleteRead,
        ) as e:  # Added IncompleteRead here
            current_delay = delay
            # Check for HTTP 429 specifically, even if wrapped
            if (
                isinstance(e, requests.exceptions.HTTPError)
                and e.response.status_code == 429
            ):
                log.warning(
                    f"Entrez call received HTTP 429 (Too Many Requests) (Attempt {attempt + 1}/{retries}). Retrying in 60s..."
                )
                current_delay = 60  # Significantly longer delay for 429
            elif isinstance(e, http.client.IncompleteRead):
                log.warning(
                    f"Entrez call resulted in IncompleteRead (Attempt {attempt + 1}/{retries}): {e}. Retrying in {current_delay}s..."
                )
            else:
                log.warning(
                    f"Entrez call network error (Attempt {attempt + 1}/{retries}): {e}. Retrying in {current_delay}s..."
                )

            if attempt + 1 == retries:
                log.error(
                    f"Entrez call failed after {retries} attempts due to network/request issues: {e}"
                )
                raise
            time.sleep(current_delay)
        except Exception as e:  # Catch other Entrez/BioPython errors
            is_http_error_str = "HTTP Error" in str(e) or "NCBI" in str(e)
            is_http_error_type = isinstance(e, (IOError, RuntimeError))

            if is_http_error_str or is_http_error_type:
                current_delay = delay
                status_code = None
                if hasattr(e, "code") and isinstance(
                    e.code, int
                ):  # For urllib.error.HTTPError
                    status_code = e.code
                elif hasattr(e, "url") and "HTTP Error" in str(
                    e
                ):  # Simple string check
                    try:
                        status_code = int(str(e).split("HTTP Error ")[1].split(":")[0])
                    except:
                        pass  # Ignore if parsing fails

                if status_code == 429:
                    log.warning(
                        f"Entrez call resulted in an NCBI/HTTP error (Attempt {attempt + 1}/{retries}) - Specifically HTTP 429: {e}. Retrying in 60s..."
                    )
                    current_delay = 60
                else:
                    log.warning(
                        f"Entrez call resulted in an NCBI/HTTP error (Attempt {attempt + 1}/{retries}): {e}. Retrying in {current_delay}s..."
                    )

                if attempt + 1 == retries:
                    log.error(
                        f"Entrez call failed after {retries} NCBI/HTTP attempts: {e}"
                    )
                    raise
                time.sleep(current_delay)
            else:  # Non-retryable error
                log.error(
                    f"Unexpected, non-retryable error during Entrez call: {e}",
                    exc_info=True,
                )
                raise
    return None


def fetch_protein_fasta_for_gene(
    gene_input: GeneInput, timeout: int, retries: int
) -> Union[str, None]:
    gene_symbol = gene_input.gene_symbol
    log.info(f"Fetching data for gene symbol: '{gene_symbol}'")

    try:
        log.debug(f"Gene '{gene_symbol}': [1/3] Fetching Gene UIDs...")
        search_term = f"{gene_symbol}[Gene Symbol]"

        handle_search = _entrez_retry_call(
            Entrez.esearch, db="gene", term=search_term, retmax="20", retries=retries
        )
        if not handle_search:
            return None

        record_search = Entrez.read(handle_search)
        handle_search.close()
        gene_ids: List[str] = record_search.get("IdList", [])

        if not gene_ids:
            log.warning(
                f"Gene '{gene_symbol}': No Gene UIDs found for symbol '{gene_symbol}' with term '{search_term}'."
            )
            log.info(
                f"Gene '{gene_symbol}': Retrying Gene UID search with broader term '{gene_symbol}[sym]'"
            )
            search_term_broad = f"{gene_symbol}[sym]"
            handle_search_broad = _entrez_retry_call(
                Entrez.esearch,
                db="gene",
                term=search_term_broad,
                retmax="20",
                retries=retries,
            )
            if not handle_search_broad:
                return None
            record_search_broad = Entrez.read(handle_search_broad)
            handle_search_broad.close()
            gene_ids = record_search_broad.get("IdList", [])
            if not gene_ids:
                log.warning(
                    f"Gene '{gene_symbol}': Still no Gene UIDs found with broader term '{search_term_broad}'."
                )
                return None

        log.debug(f"Gene '{gene_symbol}': Found Gene UIDs: {gene_ids}")

        log.debug(
            f"Gene '{gene_symbol}': [2/3] Fetching linked Protein UIDs from Gene UIDs: {gene_ids}..."
        )
        time.sleep(0.34)

        protein_ids_all: List[str] = []
        handle_elink = _entrez_retry_call(
            Entrez.elink, dbfrom="gene", db="protein", id=gene_ids, retries=retries
        )
        if not handle_elink:
            return None

        record_elink_list = Entrez.read(handle_elink)
        handle_elink.close()

        for record_elink_item in record_elink_list:
            if "LinkSetDb" in record_elink_item and record_elink_item["LinkSetDb"]:
                for link_info in record_elink_item["LinkSetDb"][0].get("Link", []):
                    protein_ids_all.append(link_info["Id"])
            elif "IdList" in record_elink_item:
                for protein_id in record_elink_item["IdList"]:
                    protein_ids_all.append(protein_id)

        if not protein_ids_all:
            log.warning(
                f"Gene '{gene_symbol}': No linked Protein UIDs found from direct gene-protein elink."
            )
            return None

        protein_ids_all = sorted(list(set(protein_ids_all)))

        log.debug(
            f"Gene '{gene_symbol}': Found {len(protein_ids_all)} unique linked Protein UIDs. E.g., {protein_ids_all[:5]}"
        )

        log.debug(
            f"Gene '{gene_symbol}': [3/3] Fetching FASTA for {len(protein_ids_all)} Protein UIDs..."
        )
        time.sleep(0.34)

        fasta_data_list: List[str] = []
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
            )
            if not handle_efetch:
                continue

            try:
                fasta_batch_data = handle_efetch.read()
            except http.client.IncompleteRead as e_read:
                log.error(
                    f"Gene '{gene_symbol}': Persistent IncompleteRead error during efetch.read() for batch {batch_ids}: {e_read}. Skipping batch."
                )
                handle_efetch.close()
                continue  # Skip this problematic batch
            finally:
                handle_efetch.close()

            handle_efetch.close()
            fasta_data_list.append(fasta_batch_data)
            if i + batch_size < len(protein_ids_all):
                time.sleep(0.34)

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
    if not keyword:
        return fasta_content_string

    log.debug(
        f"Filtering FASTA content (gene: {gene_symbol_for_log}) with keyword: '{keyword}' using at-least-2-words logic."
    )

    keyword_words = {
        word.lower() for word in keyword.split() if len(word) > 1
    }  # Ignore very short words like 'a', 'of'
    if not keyword_words:  # If keyword was only very short words or empty
        log.debug(
            f"Keyword '{keyword}' resulted in no usable words for filtering. Returning all records."
        )
        return fasta_content_string

    filtered_records: List[Any] = []
    num_total_records = 0

    try:
        for record in SeqIO.parse(StringIO(fasta_content_string), "fasta"):
            num_total_records += 1
            header_lower = record.description.lower()

            # Count matches
            matches = 0
            for kw_word in keyword_words:
                if (
                    kw_word in header_lower
                ):  # Simple substring check for each keyword word
                    matches += 1

            should_keep = False
            if len(keyword_words) == 1:
                if matches >= 1:
                    should_keep = True
            elif (
                len(keyword_words) > 1
            ):  # Covers keywords with 2 or more significant words
                if matches >= 2:
                    should_keep = True

            if should_keep:
                filtered_records.append(record)

        if not filtered_records:
            log.warning(
                f"Keyword '{keyword}' (significant words: {keyword_words}) did not match enough words in any headers for gene '{gene_symbol_for_log}'. "
                f"Original FASTA had {num_total_records} records."
            )
            return ""

        output_fasta_io = StringIO()
        SeqIO.write(filtered_records, output_fasta_io, "fasta")
        filtered_fasta_str = output_fasta_io.getvalue()
        log.debug(
            f"Keyword filtering for gene '{gene_symbol_for_log}': {len(filtered_records)}/{num_total_records} records kept with keyword '{keyword}'."
        )
        return filtered_fasta_str

    except Exception as e:
        log.error(
            f"Error during keyword filtering for gene '{gene_symbol_for_log}': {e}",
            exc_info=True,
        )
        return ""
