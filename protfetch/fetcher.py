# protfetch/fetcher.py
import time
from Bio import Entrez, SeqIO
from io import StringIO
import requests
from typing import Union, List, Any, Callable

from .utils import log, GeneInput

def configure_entrez(email: str, api_key: Union[str, None] = None):
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    log.info(f"Entrez configured with email: {email}" + (" and API key." if api_key else "."))

def _entrez_retry_call(entrez_func: Callable[..., Any], *args: Any, retries: int = 3, delay: int = 5, **kwargs: Any) -> Any:
    for attempt in range(retries):
        try:
            handle = entrez_func(*args, **kwargs)
            return handle
        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout, requests.exceptions.RequestException) as e:
            log.warning(f"Entrez call network error (Attempt {attempt + 1}/{retries}): {e}. Retrying in {delay}s...")
            if attempt + 1 == retries:
                log.error(f"Entrez call failed after {retries} network attempts: {e}")
                raise
            time.sleep(delay)
        except Exception as e:
            if "HTTP Error" in str(e) or "NCBI" in str(e) or isinstance(e, (IOError, RuntimeError)):
                 log.warning(f"Entrez call resulted in an NCBI/HTTP error (Attempt {attempt + 1}/{retries}): {e}. Retrying in {delay}s...")
                 if attempt + 1 == retries:
                    log.error(f"Entrez call failed after {retries} NCBI/HTTP attempts: {e}")
                    raise
                 time.sleep(delay)
            else:
                log.error(f"Unexpected, non-retryable error during Entrez call: {e}")
                raise
    return None

def fetch_protein_fasta_for_gene(
    gene_input: GeneInput,
    timeout: int,
    retries: int
) -> Union[str, None]:
    gene_symbol = gene_input.gene_symbol
    log.info(f"Fetching data for gene symbol: '{gene_symbol}'")

    try:
        log.debug(f"Gene '{gene_symbol}': [1/3] Fetching Gene UIDs...")
        search_term = f"{gene_symbol}[Gene Name] AND (animals[Filter] OR human[Filter]) AND alive[Prop]"
        
        handle_search = _entrez_retry_call(
            Entrez.esearch,
            db="gene",
            term=search_term,
            retmax="10",
            retries=retries
        )
        if not handle_search: return None
        
        record_search = Entrez.read(handle_search)
        handle_search.close()
        gene_ids: List[str] = record_search.get("IdList", [])

        if not gene_ids:
            log.warning(f"Gene '{gene_symbol}': No Gene UIDs found for symbol '{gene_symbol}'.")
            return None
        log.debug(f"Gene '{gene_symbol}': Found Gene UIDs: {gene_ids}")

        log.debug(f"Gene '{gene_symbol}': [2/3] Fetching linked Protein UIDs from Gene UIDs: {gene_ids}...")
        time.sleep(0.34)
        
        protein_ids_all: List[str] = []
        handle_elink = _entrez_retry_call(
            Entrez.elink,
            dbfrom="gene",
            db="protein",
            id=gene_ids,
            linkname="gene_protein_refseq",
            retries=retries
        )
        if not handle_elink: return None

        record_elink_list = Entrez.read(handle_elink)
        handle_elink.close()

        for record_elink in record_elink_list:
            link_set_db = record_elink.get("LinkSetDb")
            if link_set_db and link_set_db[0].get("Link"):
                for link in link_set_db[0]["Link"]:
                    protein_ids_all.append(link["Id"])
        
        if not protein_ids_all:
            log.warning(f"Gene '{gene_symbol}': No RefSeq Protein UIDs found. Trying broader elink.")
            time.sleep(0.34)
            handle_elink_broad = _entrez_retry_call(
                Entrez.elink, dbfrom="gene", db="protein", id=gene_ids, retries=retries
            )
            if not handle_elink_broad: return None
            record_elink_broad_list = Entrez.read(handle_elink_broad)
            handle_elink_broad.close()
            for record_elink_broad in record_elink_broad_list:
                link_set_db_broad = record_elink_broad.get("LinkSetDb")
                if link_set_db_broad and link_set_db_broad[0].get("Link"):
                    for link in link_set_db_broad[0]["Link"]:
                        protein_ids_all.append(link["Id"])
            if not protein_ids_all:
                 log.warning(f"Gene '{gene_symbol}': Still no linked Protein UIDs after broader search.")
                 return None

        log.debug(f"Gene '{gene_symbol}': Found {len(protein_ids_all)} linked Protein UIDs (pre-filter). E.g., {protein_ids_all[:5]}")

        log.debug(f"Gene '{gene_symbol}': [3/3] Fetching FASTA for {len(protein_ids_all)} Protein UIDs...")
        time.sleep(0.34)
        
        fasta_data_list: List[str] = []
        batch_size = 150 
        for i in range(0, len(protein_ids_all), batch_size):
            batch_ids = protein_ids_all[i:i+batch_size]
            log.debug(f"Gene '{gene_symbol}': Fetching FASTA batch {i//batch_size + 1} for {len(batch_ids)} IDs.")
            handle_efetch = _entrez_retry_call(
                Entrez.efetch,
                db="protein",
                id=batch_ids,
                rettype="fasta",
                retmode="text",
                retries=retries
            )
            if not handle_efetch: continue

            fasta_batch_data = handle_efetch.read()
            handle_efetch.close()
            fasta_data_list.append(fasta_batch_data)
            if i + batch_size < len(protein_ids_all):
                time.sleep(0.34)

        raw_fasta_content = "".join(fasta_data_list)

        if not raw_fasta_content.strip():
            log.warning(f"Gene '{gene_symbol}': No FASTA data returned for Protein UIDs.")
            return None

        log.info(f"Gene '{gene_symbol}': Successfully fetched raw FASTA data ({len(raw_fasta_content)} bytes).")
        return raw_fasta_content

    except Exception as e:
        log.error(f"Gene '{gene_symbol}': Error during NCBI fetch: {e}", exc_info=True)
        return None

def filter_fasta_by_keyword(
    fasta_content_string: str,
    keyword: str,
    gene_symbol_for_log: str = ""
) -> str:
    if not keyword:
        return fasta_content_string

    log.debug(f"Filtering FASTA content (gene: {gene_symbol_for_log}) with keyword: '{keyword}'")
    filtered_records: List[Any] = []
    num_total_records = 0
    try:
        for record in SeqIO.parse(StringIO(fasta_content_string), "fasta"):
            num_total_records +=1
            if keyword.lower() in record.description.lower():
                filtered_records.append(record)
        
        if not filtered_records:
            log.warning(f"Keyword '{keyword}' not found in any headers for gene '{gene_symbol_for_log}'. "
                        f"Original FASTA had {num_total_records} records.")
            return "" 

        output_fasta_io = StringIO()
        SeqIO.write(filtered_records, output_fasta_io, "fasta")
        filtered_fasta_str = output_fasta_io.getvalue()
        log.debug(f"Keyword filtering for gene '{gene_symbol_for_log}': {len(filtered_records)}/{num_total_records} records kept.")
        return filtered_fasta_str

    except Exception as e:
        log.error(f"Error during keyword filtering for gene '{gene_symbol_for_log}': {e}", exc_info=True)
        return ""
