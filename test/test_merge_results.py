"""Regression tests requiring no IQ-TREE execution.

Run from the repository root: python -m unittest discover -s test -v
"""

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from ancseq.params import Params

# The CLI parses arguments at import time; isolate that from the test runner.
with patch.object(Params, "set_options", return_value=None):
    from ancseq.ancseq import Ancseq

from check_results import check_results


class MergeResultsTest(unittest.TestCase):
    def test_fasta_matches_all_tsv_nodes(self):
        for mode, states in (
            ("DNA", ("A", "C", "T")),
            ("CODON", ("ATG", "CCA", "TTT")),
            ("AA", ("M", "P", "F")),
        ):
            for names in ((), ("Node1",), ("Node3", "Node2", "Node1")):
                with self.subTest(mode=mode, nodes=names), \
                        tempfile.TemporaryDirectory() as tmp:
                    args = SimpleNamespace(
                        seq="input.fasta", out=tmp, mode=mode, min_gap_prob=0.5
                    )
                    runner = Ancseq(args)
                    Path(runner.out_dir_10).mkdir()
                    Path(runner.out_dir_20).mkdir()
                    header = "Node\tSite\tState\tp_A\tp_C\n"
                    state_lines = ["# IQ-TREE ancestral states\n", header]
                    binary_lines = [
                        "# IQ-TREE binary states\n",
                        "Node\tSite\tState\tp_0\tp_1\n",
                    ]
                    gap = "---" if mode == "CODON" else "-"
                    expected_tsv = [f"Node\tSite\tState\tA\tC\t{gap}\n"]
                    for name in names:
                        # Exercise below, at, and above the strict gap threshold.
                        for site, (state, gap_prob) in enumerate(
                            zip(states, (0.1, 0.5, 0.9)), 1
                        ):
                            state_lines.append(
                                f"{name}\t{site}\t{state}\t0.75\t0.25\n"
                            )
                            binary_lines.append(
                                f"{name}\t{site}\t0\t{1 - gap_prob}\t{gap_prob}\n"
                            )
                            merged_state = gap if gap_prob > 0.5 else state
                            expected_tsv.append(
                                f"{name}\t{site}\t{merged_state}\t"
                                f"0.75000\t0.25000\t{gap_prob:.4f}\n"
                            )
                    Path(runner.out_dir_10, "input.fasta.state").write_text(
                        "".join(state_lines)
                    )
                    Path(runner.out_dir_20, "input.fasta.binary.state").write_text(
                        "".join(binary_lines)
                    )
                    runner.merge_results()
                    result = Path(runner.out_dir_30)
                    self.assertEqual(
                        (result / "ancestral_state_result.tsv").read_text(),
                        "".join(expected_tsv),
                    )
                    check_results(tmp)
                    runner.gzip_table()
                    check_results(tmp)


if __name__ == "__main__":
    unittest.main()
