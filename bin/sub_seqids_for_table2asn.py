#!/usr/bin/env python

import logging
from pathlib import Path

import typer
from rich.logging import RichHandler


def init_logging():
    from rich.traceback import install

    install(show_locals=True, width=120, word_wrap=True)
    logging.basicConfig(
        format="%(message)s",
        datefmt="[%Y-%m-%d %X]",
        level=logging.DEBUG,
        handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
    )


def canonicalize_seqid(seqid: str) -> str:
    """
    Normalize sequence identifiers by stripping whitespace and splitting on
    whitespace and pipe characters. This helps ensure that identifiers such as
    those in GISAID headers (which may include spaces or '|' characters) are
    consistently recognized.
    """
    return seqid.strip().split()[0].split("|")[0]


def main(
    input_fasta: Path = typer.Option(..., "--input-fasta", "-f"),
    input_feature_table: Path = typer.Option(..., "--input-feature-table", "-t"),
    outdir: Path = typer.Option(Path("./"), "--outdir", "-o"),
    prefix: str = typer.Option("SAMPLE", "--prefix", "-p"),
):
    init_logging()
    logging.info(f"{input_fasta=}")
    logging.info(f"{input_feature_table=}")
    logging.info(f"{outdir=}")
    logging.info(f"{prefix=}")

    ft_seqids = []  # List of tuples (canonical_seqid, new_seqid)
    ft_out = outdir / f"{prefix}.tbl"
    seqid_count = 0
    with input_feature_table.open() as fh, ft_out.open("w") as fout:
        for line in fh:
            if line.startswith(">Feature "):
                # Remove the leading tag and get the original header
                original = line.replace(">Feature ", "").strip()
                canon = canonicalize_seqid(original)
                seqid_count += 1
                new_seqid = f"SEQUENCE-{seqid_count}"
                ft_seqids.append((canon, new_seqid))
                # Replace only the canonical identifier in this header
                fout.write(line.replace(canon, new_seqid))
            else:
                fout.write(line)
    logging.info(f"table2asn compatible feature table written to '{ft_out}'")

    namesub_path = outdir / f"{prefix}.namesub.txt"
    with namesub_path.open("w") as fout:
        for canon, new_seqid in ft_seqids:
            fout.write(f"{new_seqid}\t{canon}\n")
    logging.info(f"Text file with original and table2asn seqids written to '{namesub_path}'")

    fasta_out = outdir / f"{prefix}.fa"
    seqid_to_new_seqid = dict(ft_seqids)
    logging.info(f"{seqid_to_new_seqid=}")
    with input_fasta.open() as fh, fasta_out.open("w") as fout:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip()
                canon = canonicalize_seqid(header)
                try:
                    fout.write(line.replace(canon, seqid_to_new_seqid[canon]))
                except KeyError:
                    logging.warning(f"{canon=} not in feature table seqids!")
            else:
                fout.write(line)
    logging.info(f"table2asn compatible fasta written to '{fasta_out}'")


if __name__ == "__main__":
    typer.run(main)