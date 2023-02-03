import os

import pyfastaq
import pymummer

from cylon import racon, utils


def get_best_contig(ref_seq, contigs_fa, outprefix):
    ref_fa = f"{outprefix}.ref.fa"
    mummer_out = f"{outprefix}.coords"
    with open(ref_fa, "w") as f:
        if not isinstance(ref_seq, pyfastaq.sequences.Fasta):
            assert isinstance(ref_seq, str)
            print(">ref", file=f)
        print(ref_seq, file=f)

    runner = pymummer.nucmer.Runner(
        ref_fa,
        contigs_fa,
        mummer_out,
        mincluster=5,
        breaklen=500,
        maxmatch=False,
    )
    runner.run()
    best_hit = None
    for hit in pymummer.coords_file.reader(mummer_out):
        if best_hit is None or hit.hit_length_qry > best_hit.hit_length_qry:
            best_hit = hit

    if best_hit is None:
        return None

    all_seqs = {}
    pyfastaq.tasks.file_to_dict(contigs_fa, all_seqs)
    all_seqs = {k.split()[0]: v for k, v in all_seqs.items()}
    if best_hit.qry_name not in all_seqs:
        return None
    contig = all_seqs[best_hit.qry_name]
    if not best_hit.on_same_strand():
        contig.revcomp()
    return contig.seq


def run_minia(
    seq_to_polish,
    reads_filename,
    outdir,
    debug=False,
    minia_opts="-keep-isolated -tip-len-topo-kmult 1.5",
):
    os.mkdir(outdir)
    verbose = "1" if debug else "0"
    command = f"minia {minia_opts} -verbose {verbose} -in {reads_filename} -out out"
    completed_process = utils.syscall(command, allow_fail=True, cwd=outdir)
    if completed_process.returncode != 0:
        return None

    contigs_fa = os.path.join(outdir, "out.contigs.fa")
    if not os.path.exists(contigs_fa):
        return None

    mummer_outprefix = os.path.join(outdir, "mummer")
    return get_best_contig(seq_to_polish, contigs_fa, mummer_outprefix)


def make_consensus(
    method,
    seq_to_polish,
    reads_filename,
    outdir,
    max_iterations=3,
    debug=False,
    minimap_opts=None,
):
    opts = {"debug": debug}
    reads_filename = os.path.abspath(reads_filename)

    if method == "minia":
        f = run_minia
    elif method == "racon":
        f = racon.run_racon_iterations
        opts.update(
            {
                "max_iterations": max_iterations,
                "minimap_opts": minimap_opts,
            }
        )
    else:
        raise NotImplementedError(f"Method not implemented: {method}")

    return f(
        seq_to_polish,
        reads_filename,
        outdir,
        **opts,
    )
