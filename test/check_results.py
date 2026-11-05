"""Check FASTA records against merged TSV states (plain or gzipped)."""

import csv
import gzip
import sys
from pathlib import Path

from Bio.SeqIO.FastaIO import SimpleFastaParser


def check_results(directory):
    result = Path(directory) / "30_result"
    table = result / "ancestral_state_result.tsv"
    opener = open
    if not table.exists():
        table = table.with_suffix(".tsv.gz")
        opener = gzip.open
    nodes = {}
    with opener(table, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            nodes.setdefault(row["Node"], []).append((int(row["Site"]), row["State"]))
    expected = [
        (node, "".join(state for _, state in sorted(sites)))
        for node, sites in nodes.items()
    ]
    for suffix, records in (
        ("_with_gap", expected),
        ("", [(node, seq.replace("-", "")) for node, seq in expected]),
    ):
        path = result / f"ancestral_state_result{suffix}.fasta"
        with path.open() as handle:
            actual = list(SimpleFastaParser(handle))
        if actual != records:
            raise AssertionError(f"{path}: FASTA records differ from TSV states")
    print(f"{directory}: tsv={len(nodes)} fasta={len(expected)} with_gap={len(expected)}")


if __name__ == "__main__":
    for directory in sys.argv[1:]:
        check_results(directory)
