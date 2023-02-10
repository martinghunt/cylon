from collections import Counter
import logging
import os
import pysam
import shutil
import subprocess
import sys

import pyfastaq
import pysam


def syscall(command, allow_fail=False, cwd=None):
    logging.debug(f"Start running command: {command}")
    completed_process = subprocess.run(
        command,
        shell=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        cwd=cwd,
    )
    if (not allow_fail) and completed_process.returncode != 0:
        print("Error running this command:", command, file=sys.stderr)
        print("Return code:", completed_process.returncode, file=sys.stderr)
        print(
            "\nOutput from stdout:", completed_process.stdout, sep="\n", file=sys.stderr
        )
        print(
            "\nOutput from stderr:", completed_process.stderr, sep="\n", file=sys.stderr
        )
        raise Exception("Error in system call. Cannot continue")

    logging.debug(
        f"Command finished with return code {completed_process.returncode}: {command}"
    )
    return completed_process


def look_for_required_binaries_in_path():
    expected_binaries = ["racon", "minimap2"]
    results = {}
    all_ok = True
    for binary in expected_binaries:
        results[binary] = shutil.which(binary)
        if results[binary] is None:
            logging.warning(f"Required dependency not found in $PATH: {binary}")
            all_ok = False
        else:
            logging.info(f"Found dependency {binary} in $PATH: {results[binary]}")
    return all_ok


def load_single_seq_fasta(infile):
    d = {}
    pyfastaq.tasks.file_to_dict(infile, d)
    if len(d) != 1:
        raise Exception(
            f"Expected exatcly 1 sequence in {infile} but got {len(d)} sequences"
        )
    ref = list(d.values())[0]
    ref.id = ref.id.split()[0]
    return ref


def rm_rf(filename):
    syscall(f"rm -rf {filename}")


def mask_low_coverage(ref_seq, reads_file, outprefix, min_depth=5, debug=False):
    ref_fasta = f"{outprefix}.fa"
    with open(ref_fasta, "w") as f:
        print(">ref", file=f)
        print(ref_seq, file=f)
    unsorted_sam = f"{outprefix}.sam"
    bam = f"{outprefix}.bam"
    syscall(f"minimap2 -x map-ont -t 1 -a {ref_fasta} {reads_file} > {unsorted_sam}")
    pysam.sort("--output-fmt", "BAM", "-o", bam, unsorted_sam)
    pysam.index(bam)

    aln_file = pysam.AlignmentFile(bam, "rb")
    new_seq = ["N"] * len(ref_seq)
    acgt = {"A", "C", "G", "T"}
    for p in aln_file.pileup(
        fastaFile=ref_fasta,
        stepper="samtools",
        compute_baq=False,
        ignore_overlaps=False,
    ):
        counts = Counter(map(str.upper, p.get_query_sequences()))
        total = 0
        max_count = 0
        best = None
        for k, v in counts.items():
            total += v
            if v > max_count:
                max_count = v
                best = k

        if total < min_depth:
            continue

        if best is not None and best in acgt:
            new_seq[p.pos] = best

    if not debug:
        os.unlink(unsorted_sam)
        os.unlink(bam)
        os.unlink(f"{bam}.bai")
        os.unlink(ref_fasta)

    return "".join(new_seq)
