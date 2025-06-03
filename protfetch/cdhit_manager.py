# protfetch/cdhit_manager.py
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .processor import ProcessedProtein
from .utils import log


def is_cdhit_installed() -> bool:
    return shutil.which("cd-hit") is not None


def get_cdhit_threshold(sequence_count: int) -> float:
    if sequence_count > 1000:
        return 0.7
    elif sequence_count >= 500:
        return 0.8
    elif sequence_count >= 250:
        return 0.85
    elif sequence_count >= 100:
        return 0.9
    else:
        return 0.95


def run_cdhit_on_proteins(
    proteins: List[ProcessedProtein],
    threshold: float,
    gene_symbol_for_log: str = "group",
) -> List[ProcessedProtein]:
    if not proteins:
        log.debug(
            f"CD-HIT: No proteins provided for {gene_symbol_for_log}, skipping CD-HIT."
        )
        return []

    if not is_cdhit_installed():
        log.error(
            "CD-HIT command not found. This function should not have been called."
        )
        return proteins

    word_size = 5
    if threshold < 0.4:
        log.warning(
            f"CD-HIT: Threshold {threshold} is very low, CD-HIT might be slow or ineffective. Word size will be 2."
        )
        word_size = 2
    elif threshold < 0.5:
        word_size = 2
    elif threshold < 0.6:
        word_size = 3
    elif threshold < 0.7:
        word_size = 4

    with tempfile.NamedTemporaryFile(
        mode="w", delete=False, suffix=".fasta"
    ) as tmp_input_fasta, tempfile.NamedTemporaryFile(
        mode="r", delete=False, suffix=".fasta"
    ) as tmp_output_fasta:

        tmp_input_fasta_path = Path(tmp_input_fasta.name)
        tmp_output_fasta_path = Path(tmp_output_fasta.name)
        tmp_output_clstr_path = Path(str(tmp_output_fasta_path) + ".clstr")

        seq_records_to_write = []
        for protein in proteins:
            seq_object = Seq(protein.sequence)
            seq_rec = SeqRecord(seq_object, id=protein.accession, description="")
            seq_records_to_write.append(seq_rec)

        SeqIO.write(seq_records_to_write, tmp_input_fasta_path, "fasta")
        tmp_input_fasta.close()

        log.info(
            f"CD-HIT: Running CD-HIT for {gene_symbol_for_log} on {len(proteins)} sequences with threshold {threshold:.2f} (word size: {word_size})..."
        )

        cmd = [
            "cd-hit",
            "-i",
            str(tmp_input_fasta_path),
            "-o",
            str(tmp_output_fasta_path),
            "-c",
            str(threshold),
            "-n",
            str(word_size),
            "-d",
            "0",
            "-M",
            "0",
            "-T",
            "0",
        ]

        final_proteins_after_cdhit = proteins
        try:
            log.debug(f"CD-HIT command: {' '.join(cmd)}")
            process = subprocess.run(cmd, capture_output=True, text=True, check=False)

            if process.returncode == 0:
                log.info(
                    f"CD-HIT for {gene_symbol_for_log} completed. stdout (first 200 chars): {process.stdout[:200]}..."
                )

                representative_accessions = set()
                if (
                    tmp_output_fasta_path.exists()
                    and tmp_output_fasta_path.stat().st_size > 0
                ):
                    for record in SeqIO.parse(tmp_output_fasta_path, "fasta"):
                        representative_accessions.add(record.id)

                final_proteins_after_cdhit = [
                    p for p in proteins if p.accession in representative_accessions
                ]
                log.info(
                    f"CD-HIT: Kept {len(final_proteins_after_cdhit)} / {len(proteins)} representative sequences for {gene_symbol_for_log}."
                )

            else:
                log.error(
                    f"CD-HIT failed for {gene_symbol_for_log} with return code {process.returncode}."
                )
                log.error(f"CD-HIT stderr: {process.stderr}")
                log.error(f"CD-HIT stdout: {process.stdout}")

        except FileNotFoundError:
            log.error(
                "CD-HIT command not found during execution. This check should happen earlier."
            )
        except Exception as e:
            log.error(
                f"An error occurred while running CD-HIT for {gene_symbol_for_log}: {e}",
                exc_info=True,
            )
        finally:
            if tmp_input_fasta_path.exists():
                tmp_input_fasta_path.unlink(missing_ok=True)
            if tmp_output_fasta_path.exists():
                tmp_output_fasta_path.unlink(missing_ok=True)
            if tmp_output_clstr_path.exists():
                tmp_output_clstr_path.unlink(missing_ok=True)

        return final_proteins_after_cdhit
