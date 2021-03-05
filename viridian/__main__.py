#!/usr/bin/env python3

import argparse
import logging
import viridian


def main(args=None):
    parser = argparse.ArgumentParser(
        prog="viridian",
        usage="viridian <command> <options>",
        description="viridian: virus amplicon assembler",
    )

    parser.add_argument("--version", action="version", version=viridian.__version__)
    parser.add_argument(
        "--debug",
        help="More verbose logging, and less file cleaning",
        action="store_true",
    )

    subparsers = parser.add_subparsers(title="Available commands", help="", metavar="")

    # ------------------------ assemble ---------------------------------------
    subparser_assemble = subparsers.add_parser(
        "assemble",
        help="Reference-guided assembly from reads",
        usage="viridian assemble [options] <ref_fasta> <amplicons_bed> <outdir>",
        description="Reference-guided assembly from reads",
        epilog="Required: --bam, or --reads_to_map, or both --reads_to_map and --mates_to_map",
    )
    subparser_assemble.add_argument(
        "ref_fasta",
        help="FASTA file of reference genome",
    )
    subparser_assemble.add_argument(
        "amplicons_bed", help="BED file of amplicon names and positions"
    )
    subparser_assemble.add_argument("outdir", help="Output directory")

    reads_group = subparser_assemble.add_argument_group("Reads options. Must use: --bam; or --reads_to_map; or --reads_to_map and --mates_to_map")
    reads_group.add_argument(
        "--bam",
        help="Input reads in a sorted indexed BAM file",
        metavar="FILENAME",
    )
    reads_group.add_argument(
        "--reads_to_map",
        help="Input reads to be mapped, in FASTA or FASTQ format",
        metavar="FILENAME",
    )
    reads_group.add_argument(
        "--mates_to_map",
        help="Input mate reads to be mapped, in FASTA or FASTQ format. If you have paired reads, use this for second file of reads",
        metavar="FILENAME",
    )
    subparser_assemble.add_argument(
        "--minimap_opts",
        help="Options string to pass to minimap2. Is used for initial mapping if reads provided with --reads_to_map, otherwise is ignored. Do not use -a or -o! This string is not sanity checked - it is up to you to provide valid options. Default is to use 1 thread and the Nanopore preset [%(default)s]",
        metavar="STRING",
        default="-t 1 -x map-ont",
    )
    subparser_assemble.add_argument(
        "--min_mean_coverage",
        help="Minimum mean read depth needed to polish an amplicon. Any amplicon less than this depth is considered failed [%(default)s]",
        default=25,
        type=int,
        metavar="INT",
    )
    subparser_assemble.add_argument(
        "--target_coverage",
        help="Aim for this much coverage when extracting reads for polishing [%(default)s]",
        default=500,
        type=int,
        metavar="INT",
    )
    subparser_assemble.add_argument(
        "--read_end_trim",
        help="Trim this many bases off the end of all reads [%(default)s]",
        default=0,
        type=int,
        metavar="INT",
    )
    subparser_assemble.add_argument(
        "--read_map_tolerance",
        help="Max allowed distance read start or end can be outside amplicon coords [%(default)s]",
        default=20,
        type=int,
        metavar="INT",
    )
    subparser_assemble.add_argument(
        "--min_read_length",
        help="Only use reads at least this long (after trimming with --read_end_trim) [%(default)s]",
        default=200,
        type=int,
        metavar="INT",
    )
    subparser_assemble.add_argument(
        "--racon_iterations",
        help="Run racon up to this many times (stops if no more corrections made) [%(default)s]",
        default=10,
        type=int,
        metavar="INT",
    )
    subparser_assemble.add_argument(
        "--min_depth_for_not_N",
        help="After polishing, each position with read depth less than this value will be masked with an N [%(default)s]",
        default=5,
        type=int,
        metavar="INT",
    )
    subparser_assemble.add_argument(
        "--min_amp_overlap_len",
        help="Minimum perfect overlap length required when overlapping each polished amplicon [%(default)s]",
        default=20,
        type=int,
        metavar="INT",
    )
    subparser_assemble.add_argument(
        "--contig_map_end_allowance",
        help="When mapping contigs ends to the reference to fill in failed amplicons with Ns, allow the mapping to start up to this distance away from the contig end [%(default)s]",
        default=20,
        type=int,
        metavar="INT",
    )
    subparser_assemble.add_argument(
        "--amplicons_to_fail_file",
        help="File of amplicon names to force to count as failed. One name in each line of the file. Names must exactly match those in amplicons_bed file",
        metavar="FILENAME",
    )
    subparser_assemble.add_argument(
        "--force",
        action="store_true",
        help="Overwrite output directory if it already exists",
    )
    subparser_assemble.set_defaults(func=viridian.tasks.assemble.run)

    # ------------------------ amplicon_overlap --------------------------------
    subparser_amp_over = subparsers.add_parser(
        "amplicon_overlap",
        help="Assemble separate amplicons seqs into consensus",
        usage="viridian assemble [options] <ref_fasta> <amplicons_bed> <amplicons_fasta> <outdir>",
        description="Assemble separate amplicons seqs into consensus",
    )
    subparser_amp_over.add_argument(
        "ref_fasta",
        help="FASTA file of reference genome",
    )
    subparser_amp_over.add_argument(
        "amplicons_bed", help="BED file of amplicon names and positions"
    )
    subparser_amp_over.add_argument(
        "amplicons_fasta",
        help="FASTA file of amplicons to be overlapped to make a consensus sequence. Names must exactly match those in BED file. Amplicons in this file must be a subset of amplicons in the BED file.",
    )
    subparser_amp_over.add_argument("outdir", help="Output directory")
    subparser_amp_over.add_argument(
        "--min_amp_overlap_len",
        help="Minimum perfect overlap length required when overlapping each polished amplicon [%(default)s]",
        default=20,
        type=int,
        metavar="INT",
    )
    subparser_amp_over.add_argument(
        "--contig_map_end_allowance",
        help="When mapping contigs ends to the reference to fill in failed amplicons with Ns, allow the mapping to start up to this distance away from the contig end [%(default)s]",
        default=20,
        type=int,
        metavar="INT",
    )
    subparser_amp_over.set_defaults(func=viridian.tasks.amplicon_overlap.run)

    args = parser.parse_args()

    logging.basicConfig(
        format="[%(asctime)s viridian %(levelname)s] %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    log = logging.getLogger()
    if args.debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()