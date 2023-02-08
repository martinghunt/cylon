import os
import pytest


from cylon import consensus_callers, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "racon")


def test_make_consensus():
    fa_to_polish = os.path.join(data_dir, "run_racon.to_polish.fa")
    seq_to_polish = utils.load_single_seq_fasta(fa_to_polish)
    reads = os.path.join(data_dir, "run_racon.reads.fa")
    outdir = "tmp.run_minia"
    utils.rm_rf(outdir)
    got_polished = consensus_callers.make_consensus(
        "minia", seq_to_polish, reads, outdir, debug=True
    )
    assert (
        got_polished
        == "CGTTAATCCTAGGGCAGTTAAAAGCCCCATTTTGTACAGCTTTTTCTAGAACAGTCAGGGCGCGCTCCCAGGAGTTGCTTCGCTTCCAGCTAGAAATGATCATCGAACCTGGGTAAGGGCATAATACGAGAATGCTGCCCTATTGCCAGTGCTTAGAAATGGACTGGTGTTACGTCCACGAAATCTGCAACAAGCCCGGT"
    )
    utils.rm_rf(outdir)

    outdir = "tmp.run_fermilite"
    utils.rm_rf(outdir)
    got_polished = consensus_callers.make_consensus(
        "fermilite", seq_to_polish, reads, outdir, debug=True
    )
    assert (
        got_polished
        == "CGTTAATCCTAGGGCAGTTAAAAGCCCCATTTTGTACAGCTTTTTCTAGAACAGTCAGGGCGCGCTCCCAGGAGTTGCTTCGCTTCCAGCTAGAAATGATCATCGAACCTGGGTAAGGGCTTAATACGAGAATGCTGCCCTATTGCCAGTGCTTAGAAATGGACTGGTGTTACGTCCACGAAATCTGCAACAAGCCCGGT"
    )
    utils.rm_rf(outdir)

    outdir = "tmp.run_racon"
    utils.rm_rf(outdir)
    got_polished = consensus_callers.make_consensus(
        "racon", seq_to_polish, reads, outdir, debug=True
    )
    assert (
        got_polished
        == "CGTTAATCCTAGGGCAGTTAAAAGCCCCATTTTGTACAGCTTTTTCTAGAACAGTCAGGGCGCGCTCCCAGGAGTTGCTTCGCTTCCAGCTAGAAATGATCATCGAACCTGGGTAAGGGCATAATACGAGAATGCTGCCCTATTGCCAGTGCTTAGAAATGGACTGGTGTTACGTCCACGAAATCTGCAACAAGCCCGGT"
    )
    utils.rm_rf(outdir)


def _test_make_consensus_bad_data():
    fa_to_polish = os.path.join(data_dir, "run_racon.to_polish.fa")
    seq_to_polish = utils.load_single_seq_fasta(fa_to_polish)
    reads = os.path.join(data_dir, "run_racon_bad_reads.fa")
    outdir = "tmp.minia_bad_data"
    utils.rm_rf(outdir)
    got_polished = consensus_callers.make_consensus(
        "minia", seq_to_polish, reads, outdir, max_iterations=3, debug=True
    )
    assert got_polished is None
    utils.rm_rf(outdir)

    outdir = "tmp.racon_bad_data"
    got_polished = consensus_callers.make_consensus(
        "racon", seq_to_polish, reads, outdir, max_iterations=3, debug=True
    )
    assert got_polished is None
    utils.rm_rf(outdir)
