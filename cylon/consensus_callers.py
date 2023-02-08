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
        if best_hit is None or hit.qry_length > best_hit.qry_length:
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


def run_minia_one_kmer(
    seq_to_polish,
    reads_filename,
    outdir,
    kmer,
    debug=False,
    minia_opts="-keep-isolated -tip-len-topo-kmult 1.5",
):
    os.mkdir(outdir)
    verbose = "1" if debug else "0"
    command = f"minia {minia_opts} -kmer-size {kmer} -verbose {verbose} -in {reads_filename} -out out"
    completed_process = utils.syscall(command, allow_fail=True, cwd=outdir)
    if completed_process.returncode != 0:
        return None

    contigs_fa = os.path.join(outdir, "out.contigs.fa")
    if not os.path.exists(contigs_fa):
        return None

    mummer_outprefix = os.path.join(outdir, "mummer")
    return get_best_contig(seq_to_polish, contigs_fa, mummer_outprefix)


def run_minia(
    seq_to_polish,
    reads_filename,
    outdir,
    kmers=[31, 13, 59],
    debug=False,
    minia_opts="-keep-isolated -tip-len-topo-kmult 1.5",
    min_length_prop=0.8,
):
    minia_seqs = []
    os.mkdir(outdir)
    for kmer in kmers:
        got = run_minia_one_kmer(
            seq_to_polish,
            reads_filename,
            os.path.join(outdir, str(kmer)),
            kmer,
            debug=debug,
            minia_opts=minia_opts,
        )
        if got is None:
            continue
        elif len(got) / len(seq_to_polish) >= min_length_prop:
            return got
        else:
            minia_seqs.append(got)

    if len(minia_seqs) > 0:
        minia_seqs.sort(key=len)
        return minia_seqs[-1]
    else:
        return None



def run_velvet(
    seq_to_polish,
    reads_filename,
    outdir,
    kmers=[11,21,31,51],
    debug=False,
    velvetg_opts="-cov_cutoff 20 -exp_cov auto",
    min_length_prop=0.8,
):
    velvet_seqs = []
    os.mkdir(outdir)
    for kmer in kmers:
        out = os.path.join(outdir, str(kmer))
        p = utils.syscall(f"velveth {out} {kmer} -short -fasta {reads_filename}", allow_fail=True)
        if p.returncode != 0:
            continue
        p = utils.syscall(f"velvetg {out} {velvetg_opts}", allow_fail=True)
        if p.returncode != 0:
            continue
        contigs_fa = os.path.join(f"{out}", "contigs.fa")
        if not os.path.exists(contigs_fa):
            return None
        got = get_best_contig(seq_to_polish, contigs_fa, f"{out}.mummer")
        if got is None:
            continue
        elif len(got) / len(seq_to_polish) >= min_length_prop:
            return got
        else:
            velvet_seqs.append(got)

    if len(velvet_seqs) > 0:
        velvet_seqs.sort(key=len)
        return velvet_seqs[-1]
    else:
        return None


def run_fml(
    seq_to_polish,
    reads_filename,
    outdir,
    debug=False,
    fml_opts="-t 1 -e 41 -A",
    min_length_prop=0.8,
):
    # Note on fermilite options: the default error correction kmer length set
    # set with -e introduces indels into the consensus. The default is "auto",
    # not sure what that is but fermikit paper suggests 23. Setting it higher
    # removes these indels.
    # Also, the -c option has a big effect depending on if the reads are all
    # stacked up in the same place or not. Default is ok, except if reads
    # are basically all the same (like in the tests), where -c1,1 works
    c_opts = ["default", "1"]
    fml_seqs = []
    os.mkdir(outdir)
    for c_opt in c_opts:
        out = os.path.join(outdir, c_opt)
        contigs_fq = f"{out}.fq"
        c = "" if c_opt =="default" else f"-c{c_opt},{c_opt}"
        p = utils.syscall(f"fml-asm {c} {fml_opts} {reads_filename} > {contigs_fq}", allow_fail=True)
        if p.returncode != 0:
            continue
        if not os.path.exists(contigs_fq):
            continue
        contigs_fa = f"{out}.fa"
        pyfastaq.tasks.to_fasta(contigs_fq, contigs_fa)
        got = get_best_contig(seq_to_polish, contigs_fa, os.path.join(outdir, "mummer"))
        if got is None:
            continue
        elif len(got) / len(seq_to_polish) >= min_length_prop:
            return got
        else:
            fml_seqs.append(got)

    if len(fml_seqs) > 0:
        fml_seqs.sort(key=len)
        return fml_seqs[-1]
    else:
        return None


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

    if method == "fermilite":
        f = run_fml
    elif method == "minia":
        f = run_minia
    elif method == "velvet":
        f = run_velvet
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
