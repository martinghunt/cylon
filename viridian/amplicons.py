from difflib import SequenceMatcher
import json
import logging
import os
import random

import pyfastaq

from viridian import racon, utils


class Amplicon:
    def __init__(self, name, start, end):
        self.name = name
        self.start = start
        self.end = end
        self.polished_seq = None
        self.masked_seq = None
        self.assemble_success = False
        self.polish_data = {
            "Reads matching": 0,
            "Reads matching forward strand": 0,
            "Reads matching reverse strand": 0,
            "Reads used": 0,
            "Coverage for polishing": 0,
            "Polish success": False,
            "Comments": [],
        }

    def __len__(self):
        return self.end - self.start + 1

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def to_dict(self):
        return {
            "name": self.name,
            "start": self.start + 1,
            "end": self.end + 1,
            "polished_seq": self.polished_seq,
            "polished_masked_seq": self.masked_seq,
            "assemble_success": self.assemble_success,
            "polish_data": self.polish_data,
        }

    def clear_seqs_because_overlap_fail(self):
        self.polished_seq = None
        self.masked_seq = None
        self.assemble_success = False
        self.polish_data["Comments"].append("No overlap with adjacent amplicon")

    def force_polish_fail(self):
        self.polish_data["Comments"].append("User chose to fail, so no processing run")

    def has_masked_seq(self):
        return self.masked_seq is not None

    def masked_seq_centre_coord(self):
        if self.has_masked_seq():
            return int(len(self.masked_seq) / 2)
        else:
            return None

    def ref_centre_coord(self):
        return self.start + int((self.end - self.start + 1) / 2)

    def get_reads_for_polishing(
        self,
        ref_name,
        bam,
        reads_outfile,
        min_coverage=25,
        trim_ends=20,
        tolerance=20,
        min_output_length=200,
        target_depth=500,
    ):
        fwd_reads = []
        rev_reads = []
        reads_used = 0

        for read in bam.fetch(ref_name, self.start, self.end + 1):
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
                or read.is_qcfail
            ):
                continue
            if (
                self.start - tolerance
                <= read.reference_start
                < read.reference_end
                <= self.end + tolerance
            ):
                # sometimes the read.query_alignment_sequence is None
                try:
                    ok = (
                        len(read.query_alignment_sequence) - 2 * trim_ends
                        >= min_output_length
                    )
                except:
                    continue
                if ok:
                    # Racon needs to have no duplicate read names.
                    # If we're parsing paired Illumina, then both reads of a
                    # pair probably have the same name. So add /2 if it's the
                    # second of a pair
                    if read.is_paired and read.is_read2:
                        read_name = read.query_name + "/2"
                    else:
                        read_name = read.query_name

                    to_append = rev_reads if read.is_reverse else fwd_reads

                    if trim_ends == 0:
                        to_append.append(
                            (
                                read_name,
                                read.get_forward_sequence(),
                            )
                        )
                    else:
                        to_append.append(
                            (
                                read_name,
                                read.get_forward_sequence()[trim_ends:-trim_ends],
                            )
                        )

        random.seed(42)
        random.shuffle(fwd_reads)
        random.shuffle(rev_reads)
        total_bp_target = target_depth * len(self)
        total_bp = 0

        if len(fwd_reads) == 0 or len(rev_reads) == 0:
            return 0, 0, 0

        with open(reads_outfile, "w") as f_out:
            # for (name, seq) in reads:
            for i in range(0, min(len(fwd_reads), len(rev_reads))):
                print(">" + fwd_reads[i][0], file=f_out)
                print(fwd_reads[i][1], file=f_out)
                print(">" + rev_reads[i][0], file=f_out)
                print(rev_reads[i][1], file=f_out)
                total_bp += len(fwd_reads[i][1]) + len(rev_reads[i][1])
                reads_used += 2
                if total_bp > total_bp_target:
                    break

        coverage = round(total_bp / len(self), 2)
        self.polish_data["Reads matching"] = len(fwd_reads) + len(rev_reads)
        self.polish_data["Reads matching forward strand"] = len(fwd_reads)
        self.polish_data["Reads matching reverse strand"] = len(rev_reads)
        self.polish_data["Reads used"] = reads_used
        self.polish_data["Coverage for polishing"] = coverage
        return len(fwd_reads) + len(rev_reads), reads_used, coverage

    def polish(
        self,
        ref_genome,
        bam,
        outdir,
        min_mean_coverage=25,
        target_coverage=500,
        read_end_trim=20,
        read_map_tolerance=20,
        min_read_length=200,
        racon_iterations=3,
        min_depth_for_not_N=5,
        max_polished_N_prop=0.1,
        debug=False,
    ):
        os.mkdir(outdir)
        reads_file = os.path.join(outdir, "reads.fa")
        total_reads, used_reads, coverage = self.get_reads_for_polishing(
            ref_genome.id,
            bam,
            reads_file,
            min_coverage=min_mean_coverage,
            trim_ends=read_end_trim,
            tolerance=read_map_tolerance,
            min_output_length=min_read_length,
            target_depth=target_coverage,
        )
        logging.debug(
            f"Extracted {total_reads} reads for amplicon {self.name}. Using {used_reads} for polishing, at mean depth of {coverage}"
        )
        if coverage < min_mean_coverage:
            logging.warning(
                f"Mean coverage for amplicon {self.name} is too low: {coverage}. Considering this a failed amplicon"
            )
            self.polish_data["Comments"].append(f"Coverage {coverage} too low")
            return
        amplicon_seq = ref_genome[self.start : self.end + 1]
        racon_dir = os.path.join(outdir, "Racon")
        self.polished_seq = racon.run_racon_iterations(
            amplicon_seq,
            reads_file,
            racon_dir,
            debug=debug,
            max_iterations=racon_iterations,
        )
        if self.polished_seq is None:
            self.polish_data["Comments"].append("No sequenced returned from racon")
            return
        mask_outprefix = os.path.join(outdir, "masked")
        self.masked_seq = utils.mask_low_coverage(
            self.polished_seq,
            reads_file,
            mask_outprefix,
            min_depth=min_depth_for_not_N,
            debug=debug,
        )
        proportion_masked = round(self.masked_seq.count("N") / len(self.masked_seq), 2)
        if proportion_masked > max_polished_N_prop:
            percent_N = 100 * proportion_masked
            self.polish_data["Comments"].append(
                f"Too many Ns ({percent_N}%) after masking polished sequence"
            )
        else:
            self.polish_data["Polish success"] = True
            self.assemble_success = True

        if not debug:
            utils.rm_rf(outdir)

    def masked_overlap(self, other, min_match_length):
        if self.masked_seq is None or other.masked_seq is None:
            return None

        # SequenceMatcher has junk=foo option, which would have been nice
        # to use to get it to not align Ns, but I couldn't get it to work and
        # it kept using them. So replace N with x or y so they can't get aligned.
        # Also use autojunk=False, otherwise the algorithm counts all the ACGTs
        # as junk and doesn't align anything.
        seq_matcher = SequenceMatcher(
            None,
            a=self.masked_seq.replace("N", "x"),
            b=other.masked_seq.replace("N", "y"),
            autojunk=False,
        )
        match = seq_matcher.find_longest_match(
            0, len(self.masked_seq), 0, len(other.masked_seq)
        )
        return match if match.size >= min_match_length else None

    def expected_overlap_length(self, other):
        if self.end < other.start:
            return None
        return self.end - other.start + 1


def load_amplicons_bed_file(infile):
    amplicons = []

    with open(infile) as f:
        for line in f:
            if line.startswith("#"):
                continue
            name, start, end = line.rstrip().split("\t")
            amplicons.append(Amplicon(name, int(start), int(end) - 1))

    return amplicons


def load_amplicons_from_fasta_and_bed(fasta, bed):
    amplicons = load_amplicons_bed_file(bed)
    amp_lookup = {a.name: i for i, a in enumerate(amplicons)}
    reader = pyfastaq.sequences.file_reader(fasta)
    for seq in reader:
        name = seq.id.split()[0]
        if name not in amp_lookup:
            raise Exception(
                f"Amplicon name {name} found in input FASTA file but not in BED file. Cannot continue"
            )
        amplicons[amp_lookup[name]].masked_seq = seq.seq
        amplicons[amp_lookup[name]].assemble_success = True
    return amplicons


def amplicons_to_json(amplicons, outfile):
    data = [a.to_dict() for a in amplicons]
    with open(outfile, "w") as f:
        json.dump(data, f, indent=2, sort_keys=True)