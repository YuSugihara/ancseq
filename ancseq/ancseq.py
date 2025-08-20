#!/usr/bin/env python3

from ancseq.params import Params


pm = Params("ancseq")
args = pm.set_options()

import os
import sys
import gzip
import shutil
import subprocess as sbp
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from ancseq.utils import time_stamp, clean_cmd, call_log


class ancseq(object):

    def __init__(self, args):
        self.args = args

    def built_tree(self):
        print(time_stamp(), "Building tree by IQ-TREE...", flush=True)
        out_dir_00 = os.path.join(self.args.out, "00_tree")
        seq = os.path.join(out_dir_00, os.path.basename(self.args.seq))
        os.mkdir(out_dir_00)
        shutil.copyfile(self.args.seq, seq)
        cmd = f"iqtree -s {seq} \
                       -st {self.args.mode} \
                       -T {self.args.threads} \
                       -m {self.args.model}"
        if self.args.fast:
            cmd += f" --fast \
                      --alrt {self.args.bootstrap}"
        else:
            cmd += f" -B {self.args.bootstrap}"
        if self.args.outgroup != None:
            cmd += f" -o {self.args.outgroup}"
        cmd += f" 1> {out_dir_00}/00_iqtree.out \
                  2> {out_dir_00}/00_iqtree.err"
        cmd = clean_cmd(cmd)
        try:
            sbp.run(cmd, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        except sbp.CalledProcessError:
            call_log(f"{out_dir_00}/00_iqtree.err", cmd)
            sys.exit(1)
        print(time_stamp(), "IQ-TREE successfully finished.", flush=True)

    def check_best_model(self):
        iqtree_log = os.path.join(
            self.args.out, "00_tree", f"{os.path.basename(self.args.seq)}.log"
        )
        with open(iqtree_log) as log:
            for line in log:
                if line.startswith("Best-fit model:"):
                    self.args.model = line.split()[2]
                    print(
                        time_stamp(),
                        f"The best model in IQ-TREE was {self.args.model}.",
                        flush=True,
                    )

    def reconstruct_ancestral_state(self):
        print(time_stamp(), "Reconstructing ancestral state...", flush=True)
        out_dir_10 = os.path.join(self.args.out, "10_asr")
        seq = os.path.join(out_dir_10, os.path.basename(self.args.seq))
        tree = os.path.join(
            self.args.out, "00_tree", f"{os.path.basename(self.args.seq)}.treefile"
        )
        os.mkdir(out_dir_10)
        shutil.copyfile(self.args.seq, seq)
        cmd = f"iqtree -asr \
                       -s {seq} \
                       -te {tree} \
                       -st {self.args.mode} \
                       -T {self.args.threads} \
                       -m {self.args.model} \
                       -keep_empty_seq"
        if self.args.outgroup != None:
            cmd += f" -o {self.args.outgroup}"
        cmd += f" 1> {out_dir_10}/10_iqtree.out \
                  2> {out_dir_10}/10_iqtree.err"
        cmd = clean_cmd(cmd)
        try:
            sbp.run(cmd, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        except sbp.CalledProcessError:
            call_log(f"{out_dir_10}/10_iqtree.err", cmd)
            sys.exit(1)
        print(
            time_stamp(),
            "Ancestral sequence reconstruction successfully finished.",
            flush=True,
        )

    def alignment_to_binary(self):
        fasta_file = open(self.args.seq)
        self.binary_file_name = os.path.join(
            self.args.out, "20_indels", f"{os.path.basename(self.args.seq)}.binary"
        )
        binary_file = open(self.binary_file_name, "w")
        for header, seq in SimpleFastaParser(fasta_file):
            bin_seq = ""
            dna_missing_chars = ["N", "X", "?", ".", "O", "~", "!"]
            aa_missing_chars = ["X", "?", ".", "~", "!"]
            if self.args.mode == "CODON":
                for i in range(0, len(seq), 3):
                    if seq[i : i + 3] == "---":
                        bin_seq += "1"
                    elif (
                        seq[i] in dna_missing_chars
                        or seq[i + 1] in dna_missing_chars
                        or seq[i + 2] in dna_missing_chars
                    ):
                        bin_seq += "-"
                    else:
                        bin_seq += "0"
            elif self.args.mode == "DNA":
                for i in range(len(seq)):
                    if seq[i] == "-":
                        bin_seq += "1"
                    elif seq[i] in dna_missing_chars:
                        bin_seq += "-"
                    else:
                        bin_seq += "0"
            elif self.args.mode == "AA":
                for i in range(len(seq)):
                    if seq[i] == "-":
                        bin_seq += "1"
                    elif seq[i] in aa_missing_chars:
                        bin_seq += "-"
                    else:
                        bin_seq += "0"
            else:
                print(
                    time_stamp(),
                    f"Invalid mode {self.args.mode} was given.",
                    file=sys.stderr,
                    flush=True,
                )
                sys.exit(1)
            binary_file.write(f">{header}\n")
            binary_file.write(f"{bin_seq}\n")
        fasta_file.close()
        binary_file.close()

    def reconstruct_indels(self):
        print(time_stamp(), "Reconstructing INDELs...", flush=True)
        out_dir_20 = os.path.join(self.args.out, "20_indels")
        tree = os.path.join(
            self.args.out, "00_tree", f"{os.path.basename(self.args.seq)}.treefile"
        )
        os.mkdir(out_dir_20)
        self.alignment_to_binary()
        cmd = f"iqtree -asr \
                       -s {self.binary_file_name} \
                       -te {tree} \
                       -st BIN \
                       -T {self.args.threads} \
                       -blfix \
                       -keep_empty_seq \
                       -m JC2"
        if self.args.outgroup != None:
            cmd += f" -o {self.args.outgroup}"
        cmd += f" 1> {out_dir_20}/20_iqtree.out \
                  2> {out_dir_20}/20_iqtree.err"
        cmd = clean_cmd(cmd)
        try:
            sbp.run(cmd, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        except sbp.CalledProcessError:
            call_log(f"{out_dir_20}/20_iqtree.err", cmd)
            sys.exit(1)
        print(time_stamp(), "INDEL reconstruction successfully finished.", flush=True)

    def merge_results(self):
        out_dir_30 = os.path.join(self.args.out, "30_result")
        os.mkdir(out_dir_30)
        probs_result_10 = open(
            os.path.join(
                self.args.out, "10_asr", f"{os.path.basename(self.args.seq)}.state"
            )
        )
        probs_result_20 = open(
            os.path.join(
                self.args.out,
                "20_indels",
                f"{os.path.basename(self.args.seq)}.binary.state",
            )
        )
        merged = open(f"{out_dir_30}/ancestral_state_result.tsv", "w")
        fasta_with_gap = open(
            f"{out_dir_30}/ancestral_state_result_with_gap.fasta", "w"
        )
        fasta = open(f"{out_dir_30}/ancestral_state_result.fasta", "w")
        if self.args.mode == "CODON":
            gap = "---"
        else:
            gap = "-"
        for line_10, line_20 in zip(probs_result_10, probs_result_20):
            if not line_10.startswith("#"):
                line_10 = line_10.rstrip("\n")
                cols_10 = line_10.split("\t")
                if not line_10.startswith("Node\tSite"):
                    line_20 = line_20.rstrip("\n")
                    cols_20 = line_20.split("\t")
                    node_10 = cols_10[0]
                    node_20 = cols_20[0]
                    site_10 = cols_10[1]
                    site_20 = cols_20[1]
                    state = cols_10[2]
                    probs = [float(p) for p in cols_10[3:]]
                    p_0 = float(cols_20[3])
                    p_1 = float(cols_20[4])
                    if (node_10, site_10) == (node_20, site_20):
                        if previous_node != node_10 and previous_node != None:
                            fasta_with_gap.write(f">{previous_node}\n")
                            fasta_with_gap.write(f"{node_seq}\n")
                            fasta.write(f">{previous_node}\n")
                            fasta.write("{}\n".format(node_seq.replace("-", "")))
                            node_seq = ""
                        if p_1 > self.args.min_gap_prob:
                            merged.write(
                                "{}\t{}\t{}\t{}\t{:.4f}\n".format(
                                    node_10,
                                    site_10,
                                    gap,
                                    "\t".join(
                                        [
                                            format(round(prob, 4), ".5f")
                                            for prob in probs
                                        ]
                                    ),
                                    p_1,
                                )
                            )
                            node_seq = node_seq + gap
                        else:
                            merged.write(
                                "{}\t{}\t{}\t{}\t{:.4f}\n".format(
                                    node_10,
                                    site_10,
                                    state,
                                    "\t".join(
                                        [
                                            format(round(prob, 4), ".5f")
                                            for prob in probs
                                        ]
                                    ),
                                    p_1,
                                )
                            )
                            node_seq = node_seq + state
                        previous_node = node_10
                    else:
                        print(line_10)
                        print(line_20)
                        print(
                            time_stamp(),
                            f"{probs_result_10} and {probs_result_20} did not much.",
                            file=sys.stderr,
                            flush=True,
                        )
                        sys.exit(1)
                else:
                    merged.write(
                        "Node\tSite\tState\t{}\t{}\n".format(
                            "\t".join([p.replace("p_", "") for p in cols_10[3:]]), gap
                        )
                    )
                    previous_node = None
                    node_seq = ""
        probs_result_10.close()
        probs_result_20.close()
        fasta_with_gap.close()
        fasta.close()

    def calculate_codon_prob(self):
        print(time_stamp(), "Calcurating the probability of codon...", flush=True)
        out_dir_30 = os.path.join(self.args.out, "30_result")
        merged = open(f"{out_dir_30}/ancestral_state_result.tsv")
        codon_prob = open(f"{out_dir_30}/ancestral_state_result.codon_prob.tsv", "w")
        for line in merged:
            line = line.rstrip("\n")
            cols = line.split("\t")
            if line.startswith("Node\tSite"):
                codon_prob.write(
                    "Node\tSite\tState\t{}\t---\n".format(
                        "\t".join(
                            [
                                "{}{}{}".format(nuc_1, nuc_2, nuc_3)
                                for nuc_1 in ["A", "C", "G", "T"]
                                for nuc_2 in ["A", "C", "G", "T"]
                                for nuc_3 in ["A", "C", "G", "T"]
                            ]
                        )
                    )
                )
            else:
                site = int(cols[1])
                if site % 3 == 1:
                    nuc_dict_1 = {
                        nuc: float(p) for nuc, p in zip(["A", "C", "G", "T"], cols[3:7])
                    }
                    gap_1 = float(cols[7])
                elif site % 3 == 2:
                    nuc_dict_2 = {
                        nuc: float(p) for nuc, p in zip(["A", "C", "G", "T"], cols[3:7])
                    }
                    gap_2 = float(cols[7])
                else:
                    nuc_dict_3 = {
                        nuc: float(p) for nuc, p in zip(["A", "C", "G", "T"], cols[3:7])
                    }
                    gap_3 = float(cols[7])
                    codon_probs = {
                        "{}{}{}".format(nuc_1, nuc_2, nuc_3): nuc_dict_1[nuc_1]
                        * nuc_dict_2[nuc_2]
                        * nuc_dict_3[nuc_3]
                        for nuc_1 in ["A", "C", "G", "T"]
                        for nuc_2 in ["A", "C", "G", "T"]
                        for nuc_3 in ["A", "C", "G", "T"]
                    }
                    codon_probs["---"] = (gap_1 + gap_2 + gap_3) / 3
                    codon_sorted = sorted(
                        codon_probs.items(), key=lambda x: x[1], reverse=True
                    )
                    if codon_probs["---"] >= self.args.min_gap_prob:
                        codon_prob.write(
                            "{}\t{}\t{}\t{}\n".format(
                                cols[0],
                                int(site / 3),
                                "---",
                                "\t".join(
                                    [
                                        format(round(v, 4), ".5f")
                                        for v in codon_probs.values()
                                    ]
                                ),
                            )
                        )
                    else:
                        codon_prob.write(
                            "{}\t{}\t{}\t{}\n".format(
                                cols[0],
                                int(site / 3),
                                codon_sorted[0][0],
                                "\t".join(
                                    [
                                        format(round(v, 4), ".5f")
                                        for v in codon_probs.values()
                                    ]
                                ),
                            )
                        )
        merged.close()
        codon_prob.close()
        print(time_stamp(), "Calcuration of codon probabilities finished.", flush=True)

    def sort_ambiguous_codon(self):
        print(time_stamp(), "Sorting ambiguous codons...", flush=True)
        if self.args.mode == "DNA":
            out_dir_30 = os.path.join(self.args.out, "30_result")
            ancestral_state = open(
                f"{out_dir_30}/ancestral_state_result.codon_prob.tsv"
            )
        elif self.args.mode == "CODON":
            out_dir_30 = os.path.join(self.args.out, "30_result")
            ancestral_state = open(f"{out_dir_30}/ancestral_state_result.tsv")
        ambiguous_site = open(f"{out_dir_30}/ancestral_state_result.sort.tsv", "w")
        for line in ancestral_state:
            line = line.rstrip("\n")
            cols = line.split("\t")
            if line.startswith("Node\tSite"):
                codons = cols[3:]
            else:
                node = cols[0]
                site = cols[1]
                codon_probs = {
                    codon: float(prob) for codon, prob in zip(codons, cols[3:])
                }
                if codon_probs["---"] > self.args.min_gap_prob:
                    top_codons = ["---"] + (self.args.max_report - 1) * [""]
                    top_codon_probs = [format(round(codon_probs["---"], 4), ".5f")] + (
                        self.args.max_report - 1
                    ) * [""]
                    top_aas = ["-"] + (self.args.max_report - 1) * [""]
                    top_aa_probs = [format(round(codon_probs["---"], 4), ".5f")] + (
                        self.args.max_report - 1
                    ) * [""]
                else:
                    del codon_probs["---"]
                    aa_probs = {}
                    for k, v in codon_probs.items():
                        aa = Seq(k).translate()[0]
                        if aa not in aa_probs:
                            aa_probs[aa] = v
                        else:
                            aa_probs[aa] = aa_probs[aa] + v
                    codon_sorted = sorted(
                        codon_probs.items(), key=lambda x: x[1], reverse=True
                    )
                    aa_sorted = sorted(
                        aa_probs.items(), key=lambda x: x[1], reverse=True
                    )
                    top_codons = [
                        k if v > self.args.min_prob else "" for k, v in codon_sorted
                    ][: self.args.max_report]
                    top_codon_probs = [
                        format(round(v, 4), ".5f") if v > self.args.min_prob else ""
                        for k, v in codon_sorted
                    ][: self.args.max_report]
                    top_aas = [
                        k if v > self.args.min_prob else "" for k, v in aa_sorted
                    ][: self.args.max_report]
                    top_aa_probs = [
                        format(round(v, 4), ".5f") if v > self.args.min_prob else ""
                        for k, v in aa_sorted
                    ][: self.args.max_report]
                ambiguous_site.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        node,
                        site,
                        "\t".join(top_codons),
                        "\t".join(top_codon_probs),
                        "\t".join(top_aas),
                        "\t".join(top_aa_probs),
                    )
                )
        ancestral_state.close()
        ambiguous_site.close()
        print(time_stamp(), "Ambiguous codons were sorted.", flush=True)

    def sort_ambiguous_aa(self):
        print(time_stamp(), "Sorting ambiguous amino acids...", flush=True)
        out_dir_30 = os.path.join(self.args.out, "30_result")
        ancestral_state = open(f"{out_dir_30}/ancestral_state_result.tsv")
        ambiguous_site = open(f"{out_dir_30}/ancestral_state_result.sort.tsv", "w")
        for line in ancestral_state:
            line = line.rstrip("\n")
            cols = line.split("\t")
            if line.startswith("Node\tSite"):
                aas = cols[3:]
            else:
                node = cols[0]
                site = cols[1]
                aa_probs = {aa: float(prob) for aa, prob in zip(aas, cols[3:])}
                if aa_probs["-"] > self.args.min_gap_prob:
                    top_aas = ["-"] + (self.args.max_report - 1) * [""]
                    top_aa_probs = [format(round(aa_probs["-"], 4), ".5f")] + (
                        self.args.max_report - 1
                    ) * [""]
                else:
                    del aa_probs["-"]
                    aa_sorted = sorted(
                        aa_probs.items(), key=lambda x: x[1], reverse=True
                    )
                    top_aas = [
                        k if v > self.args.min_prob else "" for k, v in aa_sorted
                    ][: self.args.max_report]
                    top_aa_probs = [
                        format(round(v, 4), ".5f") if v > self.args.min_prob else ""
                        for k, v in aa_sorted
                    ][: self.args.max_report]
                ambiguous_site.write(
                    "{}\t{}\t{}\t{}\n".format(
                        node, site, "\t".join(top_aas), "\t".join(top_aa_probs)
                    )
                )
        ancestral_state.close()
        ambiguous_site.close()
        print(time_stamp(), "Ambiguous amino acids were sorted.", flush=True)

    def copy_treefile(self):
        out_dir_10 = os.path.join(self.args.out, "10_asr")
        out_dir_30 = os.path.join(self.args.out, "30_result")
        treefile = os.path.join(
            out_dir_10, f"{os.path.basename(self.args.seq)}.treefile"
        )
        shutil.copyfile(treefile, f"{out_dir_30}/ancestral_state_result.treefile")

    def gzip_table(self):
        out_dir_10 = os.path.join(self.args.out, "10_asr")
        out_dir_20 = os.path.join(self.args.out, "20_indels")
        out_dir_30 = os.path.join(self.args.out, "30_result")
        statefile_10 = os.path.join(
            out_dir_10, f"{os.path.basename(self.args.seq)}.state"
        )
        statefile_20 = os.path.join(
            out_dir_20, f"{os.path.basename(self.args.seq)}.binary.state"
        )
        marged_prefix = os.path.join(out_dir_30, "ancestral_state_result")
        if os.path.isfile(statefile_10):
            with open(statefile_10, "rb") as f_in:
                with gzip.open(f"{statefile_10}.gz", "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(statefile_10)
        if os.path.isfile(statefile_20):
            with open(statefile_20, "rb") as f_in:
                with gzip.open(f"{statefile_20}.gz", "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(statefile_20)
        if os.path.isfile(f"{marged_prefix}.tsv"):
            with open(f"{marged_prefix}.tsv", "rb") as f_in:
                with gzip.open(f"{marged_prefix}.tsv.gz", "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(f"{marged_prefix}.tsv")
        if os.path.isfile(f"{marged_prefix}.codon_prob.tsv"):
            with open(f"{marged_prefix}.codon_prob.tsv", "rb") as f_in:
                with gzip.open(f"{marged_prefix}.codon_prob.tsv.gz", "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(f"{marged_prefix}.codon_prob.tsv")

    def run(self):
        print(time_stamp(), "Start to run ancseq.", flush=True)
        os.mkdir(self.args.out)
        if self.args.asr_only:
            print(time_stamp(), "Skipped building tree.", flush=True)
            out_dir_00 = os.path.join(self.args.out, "00_tree")
            seq = os.path.join(out_dir_00, os.path.basename(self.args.seq))
            os.mkdir(out_dir_00)
            shutil.copyfile(self.args.seq, seq)
            shutil.copyfile(f"{self.args.seq}.treefile", f"{seq}.treefile")
            if os.path.isfile(f"{self.args.seq}.log"):
                shutil.copyfile(f"{self.args.seq}.log", f"{seq}.log")
        else:
            self.built_tree()
        if self.args.model == "MFP" and os.path.isfile(
            f"{self.args.out}/00_tree/{self.args.seq}.log"
        ):
            self.check_best_model()
        self.reconstruct_ancestral_state()
        self.reconstruct_indels()
        self.merge_results()
        if self.args.mode == "DNA" and not self.args.stop_codon_prob:
            self.calculate_codon_prob()
            self.sort_ambiguous_codon()
        elif self.args.mode == "CODON":
            self.sort_ambiguous_codon()
        else:
            self.sort_ambiguous_aa()
        self.copy_treefile()
        self.gzip_table()
        print(time_stamp(), "ancseq successfully finished!", flush=True)


def main():
    ancseq(args).run()


if __name__ == "__main__":
    main()
