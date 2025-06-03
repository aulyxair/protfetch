import http.client
import time
from io import StringIO
from typing import Any, Callable, Dict, List, Union

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
    parser_func: Union[
        Callable[[Any], Any], None
    ] = None,  # Optional parser like Entrez.read
    *args: Any,
    retries: int = 3,
    delay: int = 5,
    **kwargs: Any,
) -> Any:
    last_exception = None
    for attempt in range(retries):
        handle = None
        try:
            handle = entrez_func(*args, **kwargs)
            if parser_func:
                # If a parser is provided, attempt to parse here. Errors will be caught.
                result = parser_func(handle)
                return result  # Return parsed result
            else:
                return handle
        except (
            requests.exceptions.ConnectionError,
            requests.exceptions.Timeout,
            requests.exceptions.RequestException,
            http.client.IncompleteRead,
        ) as e:
            last_exception = e
            current_delay = delay
            if (
                isinstance(e, requests.exceptions.HTTPError)
                and hasattr(e, "response")
                and e.response is not None
                and e.response.status_code == 429
            ):
                log.warning(
                    f"Entrez call received HTTP 429 (Too Many Requests) (Attempt {attempt + 1}/{retries}). Retrying in 60s..."
                )
                current_delay = 60
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
                    f"Entrez call failed after {retries} attempts due to network/request issues: {last_exception}"
                )
                raise last_exception
            time.sleep(current_delay)
        except Exception as e:
            last_exception = e
            is_http_error_str = (
                "HTTP Error" in str(e) or "NCBI" in str(e) or "callMLink" in str(e)
            )  # Added callMLink
            is_http_error_type = isinstance(e, (IOError, RuntimeError))

            if is_http_error_str or is_http_error_type:
                current_delay = delay
                status_code = None
                if hasattr(e, "code") and isinstance(e.code, int):
                    status_code = e.code
                elif hasattr(e, "url") and "HTTP Error" in str(e):
                    try:
                        status_code = int(str(e).split("HTTP Error ")[1].split(":")[0])
                    except:
                        pass
                elif hasattr(e, "response") and hasattr(e.response, "status_code"):
                    status_code = e.response.status_code

                if status_code == 429:
                    log.warning(
                        f"Entrez call resulted in an NCBI/HTTP error (Attempt {attempt + 1}/{retries}) - Specifically HTTP 429: {e}. Retrying in 60s..."
                    )
                    current_delay = 60
                else:
                    log.warning(
                        f"Entrez call resulted in an NCBI/HTTP error or parse error (Attempt {attempt + 1}/{retries}): {e}. Retrying in {current_delay}s..."
                    )

                if attempt + 1 == retries:
                    log.error(
                        f"Entrez call failed after {retries} NCBI/HTTP or parse attempts: {last_exception}"
                    )
                    raise last_exception
                time.sleep(current_delay)
            else:
                log.error(
                    f"Unexpected, non-retryable error during Entrez call or parsing: {e}",
                    exc_info=True,
                )
                raise
        finally:
            if handle and parser_func and hasattr(handle, "close"):
                handle.close()

    if last_exception:
        log.error(f"All Entrez retries failed. Last error: {last_exception}")
        raise last_exception  # Re-raise the last known exception
    return None  # Fallback, though should be an exception


def fetch_protein_fasta_for_gene(
    gene_input: GeneInput, timeout: int, retries: int
) -> Union[str, None]:
    gene_symbol = gene_input.gene_symbol
    log.info(f"Fetching data for gene symbol: '{gene_symbol}'")

    gene_search_webenv: Union[str, None] = None
    gene_search_query_key: Union[str, None] = None
    gene_search_count: int = 0

    try:
        log.debug(f"Gene '{gene_symbol}': [1/3] Searching for Gene UIDs...")
        search_term = f"{gene_symbol}[symbol]"

        search_results = _entrez_retry_call(
            Entrez.esearch,
            parser_func=Entrez.read,
            db="gene",
            term=search_term,
            usehistory="y",
            retries=retries,
        )
        if not search_results:
            return None

        gene_search_count = int(search_results["Count"])
        gene_search_webenv = search_results.get("WebEnv")
        gene_search_query_key = search_results.get("QueryKey")

        log.info(
            f"Gene '{gene_symbol}': Initial esearch with '{search_term}' found {gene_search_count} Gene UIDs. WebEnv: {gene_search_webenv}"
        )

        if gene_search_count == 0:
            log.warning(
                f"Gene '{gene_symbol}': No Gene UIDs found with primary term '{search_term}'. Retrying with '{gene_symbol}[Gene Symbol]'."
            )
            search_term_specific = f"{gene_symbol}[Gene Symbol]"
            search_results_specific = _entrez_retry_call(
                Entrez.esearch,
                parser_func=Entrez.read,
                db="gene",
                term=search_term_specific,
                usehistory="y",
                retries=retries,
            )
            if not search_results_specific:
                return None

            gene_search_count = int(search_results_specific["Count"])
            gene_search_webenv = search_results_specific.get("WebEnv")
            gene_search_query_key = search_results_specific.get("QueryKey")
            if gene_search_count == 0:
                log.warning(
                    f"Gene '{gene_symbol}': Still no Gene UIDs found with fallback term '{search_term_specific}'."
                )
                return None
            log.info(
                f"Gene '{gene_symbol}': Fallback esearch found {gene_search_count} Gene UIDs. WebEnv: {gene_search_webenv}"
            )

        if not gene_search_webenv or not gene_search_query_key:
            log.error(
                f"Gene '{gene_symbol}': Failed to obtain WebEnv/QueryKey from esearch despite finding {gene_search_count} UIDs."
            )
            return None

        log.info(
            f"Gene '{gene_symbol}': Proceeding to elink {gene_search_count} Gene UIDs from history."
        )
        log.debug(
            f"Gene '{gene_symbol}': [2/3] Linking Gene UIDs from history to Protein UIDs..."
        )
        time.sleep(0.34)

        record_elink_list = _entrez_retry_call(
            Entrez.elink,
            parser_func=Entrez.read,
            dbfrom="gene",
            db="protein",
            WebEnv=gene_search_webenv,
            query_key=gene_search_query_key,
            retries=retries,
        )

        if not record_elink_list:
            return None

        protein_ids_all: List[str] = []
        for record_elink_item in record_elink_list:
            if "LinkSetDb" in record_elink_item and record_elink_item["LinkSetDb"]:
                for link_set_db_entry in record_elink_item["LinkSetDb"]:
                    if link_set_db_entry.get("Link"):
                        for link_info in link_set_db_entry.get("Link", []):
                            protein_ids_all.append(link_info["Id"])
            elif "IdList" in record_elink_item:
                for protein_id in record_elink_item["IdList"]:
                    protein_ids_all.append(protein_id)

        if not protein_ids_all:
            log.warning(
                f"Gene '{gene_symbol}': No linked Protein UIDs found from gene-protein elink using history."
            )
            return None

        protein_ids_all = sorted(list(set(protein_ids_all)))

        log.info(
            f"Gene '{gene_symbol}': Found {len(protein_ids_all)} unique linked Protein UIDs. E.g., {protein_ids_all[:5]}"
        )
        if len(protein_ids_all) == 0:
            log.warning(
                f"Gene '{gene_symbol}': Zero unique protein UIDs after elink and deduplication."
            )
            return None

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
                if hasattr(handle_efetch, "close"):
                    handle_efetch.close()
                continue
            finally:
                if (
                    "handle_efetch" in locals()
                    and hasattr(handle_efetch, "close")
                    and not handle_efetch.closed
                ):
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

    lower_keyword = keyword.lower()
    log.debug(
        f"Filtering FASTA content (gene: {gene_symbol_for_log}) by checking if entire keyword phrase '{lower_keyword}' is in header."
    )

    filtered_records: List[Any] = []
    num_total_records = 0

    try:
        for record in SeqIO.parse(StringIO(fasta_content_string), "fasta"):
            num_total_records += 1
            header_lower = record.description.lower()

            if lower_keyword in header_lower:
                filtered_records.append(record)

        if not filtered_records:
            log.warning(
                f"Keyword phrase '{keyword}' not found as a substring in any headers for gene '{gene_symbol_for_log}'. "
                f"Original FASTA had {num_total_records} records."
            )
            return ""

        output_fasta_io = StringIO()
        SeqIO.write(filtered_records, output_fasta_io, "fasta")
        filtered_fasta_str = output_fasta_io.getvalue()
        log.info(
            f"Keyword filtering for gene '{gene_symbol_for_log}': {len(filtered_records)}/{num_total_records} records kept using keyword phrase '{keyword}'."
        )
        return filtered_fasta_str

    except Exception as e:
        log.error(
            f"Error during keyword filtering for gene '{gene_symbol_for_log}': {e}",
            exc_info=True,
        )
        return ""
