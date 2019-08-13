from collections import defaultdict, namedtuple
import operator
from itertools import tee
import re

from Bio import Alphabet, SeqIO
from Bio.Seq import Seq
import hgvs.edit, hgvs.location, hgvs.posedit, hgvs.sequencevariant
from hgvs.utils.altseq_to_hgvsp import AltSeqToHgvsp
from pyfaidx import Fasta
import vcfpy

from .utils import GenomePosition, GenomeIntervalTree, vcf_alt_affected_range

_VariantTreeTranscriptRecord = namedtuple("VariantTreeTranscriptRecord", [
    "feat", "five_prime_cds_feat", "three_prime_cds_feat",
    "five_prime_utr_feats", "three_prime_utr_feats", "translated_feats"
])

_VariantTreeCodingRecord = namedtuple("VariantTreeCodingRecord",
                                      ["transcript", "feat"])

_LocalizedSeq = namedtuple("LocalizedSeq", ["seq", "pos"])

# Create a mock `RefTranscriptData` type which behaves sufficiently similarly
# for the purposes of this module. Used to avoid the initialization step of the
# actual type.
_MockRefTranscriptData = namedtuple("RefTranscriptData", [
    "transcript_sequence", "aa_sequence", "cds_start", "cds_stop",
    "protein_accession"
])

# Create a mock `AltTranscriptData` type which behaves sufficiently similarly
# for the purposes of this module. Used to avoid the initialization step of the
# actual type.
_MockAltTranscriptData = namedtuple("AltTranscriptData", [
    "transcript_sequence", "aa_sequence", "cds_start", "cds_stop",
    "protein_accession", "is_frameshift", "variant_start_aa",
    "frameshift_start", "is_substitution", "is_ambiguous"
])

ProteinVariantResult = namedtuple(
    "ProteinVariantResult",
    ["query_variant", "predicted_variant", "transcript_feat"])


def _merge_and_sort_intersecting_intervals(intervals):
    current_interval = None

    for interval in sorted(intervals):
        if current_interval is None:
            current_interval = interval
            continue

        i0, i1, c0, c1 = *interval, *current_interval

        if i0 > c1:  # i1 < c0 shouldn't be possible due to sort
            # New interval
            yield current_interval
            current_interval = interval
            continue

        # This interval intersects with the current one; merge them
        current_interval = (min(i0, c0), max(i1, c1))

    if current_interval is not None:
        yield current_interval


class ProteinSequenceVariantMatcher():
    def __init__(self, alt_type, ref_seq, alt_seq):
        self.alt_type = alt_type
        self.ref_seq = ref_seq
        self.alt_seq = alt_seq

        self.ref_len = len(self.ref_seq)
        self.alt_len = len(self.alt_seq)

    # All posedits should be marked `uncertain=True`, as they are ultimately
    # predictions.
    def create_hgvs_posedit(self):
        return self._find_posedit()

    def _find_first_diff_index(self):
        for i in range(min(self.ref_len, self.alt_len)):
            if self.ref_seq[i] != self.alt_seq[i]:
                return i

        # Since the end of the overlapping sequences has been reached, if the
        # sequences differ in length, the first difference must occur at the end
        # of the shorter sequence.
        if self.ref_len != self.alt_len:
            return min(self.ref_len, self.alt_len)

        # The ref and alt sequences are the same.
        return -1

    def _find_posedit(self):
        """Try to find a posedit which heuristically looks like an insertion."""
        first_diff_idx = self._find_first_diff_index()

        seq_len_delta = self.alt_len - self.ref_len
        diff_slice = slice(first_diff_idx, first_diff_idx + abs(seq_len_delta))

        # TODO: Check for extensions.
        # TODO: Consider: "insertions from intron sequences which maintain the normal open reading are described as insertion or Deletion-insertion, intronic insertions which give premature translation termination are described as frame shift"

        diff = self.alt_seq[diff_slice]

        # If there is a difference between the ref and alt past the point of
        # insertion or a stop codon in the insertion, a frameshift likely
        # occurred.
        if (self.ref_seq[diff_slice.stop:] !=
                self.alt_seq[diff_slice.stop:]) or ('*' in diff):
            return self._find_fs_posedit(diff_slice)

        start_ref_aa = self.ref_seq[diff_slice.start - 1]
        end_ref_aa = self.ref_seq[diff_slice.start]

        # 1-indexed
        start = hgvs.location.AAPosition(base=diff_slice.start,
                                         aa=start_ref_aa)
        end = hgvs.location.AAPosition(base=diff_slice.start + 1,
                                       aa=end_ref_aa)
        ival = hgvs.location.Interval(start, end)

        edit = hgvs.edit.AARefAlt(ref=None, alt=str(diff))

        return hgvs.posedit.PosEdit(pos=ival, edit=edit, uncertain=True)

    def _find_fs_posedit(self, diff_slice):
        ref_aa = self.ref_seq[diff_slice.start]
        alt_aa = self.alt_seq[diff_slice.start]

        fs_len = self.alt_seq[diff_slice.start:].find('*')

        # If no termination codon found, default to remaining sequence length.
        if fs_len < 0:
            fs_len = self.alt_len - diff_slice.stop

        # 1-indexed
        pos = hgvs.location.AAPosition(base=diff_slice.start + 1, aa=ref_aa)
        edit = hgvs.edit.AAFs(ref=ref_aa, alt=alt_aa, length=fs_len + 1)

        return hgvs.posedit.PosEdit(pos=pos, edit=edit, uncertain=True)


class ProteinSequenceVariantGenerator():
    def __init__(self, ref_seq, alt_seq):
        self.ref_seq = ref_seq
        self.alt_seq = alt_seq
        # self.alt_offset = alt_offset

        # def get_aa(seq, pos, seq_offset=0):
        #     seq_idx = pos - seq_offset

        #     if seq_idx < 0 or seq_idx >= len(seq):
        #         return None

        #     return seq[seq_idx]

        # aa_seq_len = max(len(self.ref_seq),
        #                  len(self.alt_seq) + self.alt_offset)
        self.aa_seq = [(pos, self.ref_seq[pos], self.alt_seq[pos])
                       for pos in range(len(self.ref_seq) - 1, -1, -1)]

    def generate_hgvs_posedit(self):
        # FIXME: Python 3.8 assignment expressions would be really nice here...

        cterm_ext = self._try_extract_cterm_ext()
        if cterm_ext: return cterm_ext

        nterm_ext = self._try_extract_nterm_ext()
        if nterm_ext: return nterm_ext

    def _try_extract_cterm_ext(self):
        ext_ter_pos = None
        ter_pos = None
        alt_ter_aa = None

        for pos, ref_aa, alt_aa in self.aa_seq:
            if ext_ter_pos and alt_aa != '_' and ref_aa == '*':
                ter_pos = pos
                alt_ter_aa = alt_aa
                break

            if ref_aa == '_' and alt_aa == '*':
                ext_ter_pos = pos
        else:
            # The ref didn't contain a stop codon or one was not found prior to
            # it in the alt.
            return

        ext_len = ext_ter_pos - ter_pos
        pos = hgvs.location.AAPosition(base=ter_pos + 1, aa='*')
        edit = hgvs.edit.AAExt(ref='*', alt=alt_ter_aa, length=ext_len)

        return hgvs.posedit.PosEdit(pos=pos, edit=edit, uncertain=True)

    def _try_extract_nterm_ext(self):
        ext_met_pos = None
        met_pos = None
        alt_met_aa = None

        # rethink this
        for pos, ref_aa, alt_aa in self.aa_seq:
            if ext_met_pos and alt_aa != '_' and ref_aa == 'M':
                met_pos = pos
                alt_met_aa = alt_aa

            if ref_aa == '_' and alt_aa == 'M':
                ext_met_pos = pos
        else:
            # The ref didn't contain a start codon or one was not found prior to
            # it in the alt.
            return

        ext_len = ext_met_pos - met_pos
        pos = hgvs.location.AAPosition(base=met_pos + 1, aa='M')
        edit = hgvs.edit.AAExt(ref='M', alt=alt_met_aa, length=ext_len)

        return hgvs.posedit.PosEdit(pos=pos, edit=edit, uncertain=True)


class AAVariantPredictor():
    # TODO: Now I think it would be nice if this was written to use GFF3 instead
    # of GTF (hierarchic features) and the entire genome transcript rather than
    # fragmented protein-coding transcripts. However, using the entire genome
    # transcript would require a more sophisticated method for querying regions...

    def __init__(self, refgenome_tree, genome_fasta_file):
        # transcript_reader = SeqIO.parse(protein_coding_transcript, "fasta")
        self.genome_fasta = Fasta(genome_fasta_file)

        transcript_feats_dict = defaultdict(lambda: defaultdict(list))
        for feature in refgenome_tree.records:
            if feature.attributes.get("transcript_type") != "protein_coding":
                continue

            feature_type = feature.type
            if feature_type not in [
                    "transcript", "CDS", "UTR", "start_codon", "stop_codon"
            ]:
                continue

            transcript_id = feature.attributes["transcript_id"]
            transcript_feats_dict[transcript_id][feature_type].append(feature)

        transcript_records = {}
        cds_records = []
        for transcript_id, features in transcript_feats_dict.items():
            # There should only be one transcript record...
            transcript_feat = features["transcript"][0]

            # Hopefully redundant check
            if not transcript_feat:
                print(
                    f"Unexpected missing transcript feature for transcript ID {transcript_id}"
                )
                continue

            cds_feats = features["CDS"]

            # Sort the CDS features so that the CDS adjacent to the 5' UTR is
            # first and the CDS adjacent to the 3' UTR is last.
            is_forward_stranded = transcript_feat.is_forward_stranded
            cds_feats.sort(reverse=(not is_forward_stranded),
                           key=lambda feat: feat.pos.start
                           if is_forward_stranded else feat.pos.end)

            five_prime_cds_feat = cds_feats[0]
            three_prime_cds_feat = cds_feats[-1]

            utr_feats = features["UTR"]

            is_before = lambda rh: lambda lh: lh.pos.end <= rh.pos.start
            is_after = lambda rh: lambda lh: is_before(lh)(rh)

            utr_filter_5p = is_before(
                five_prime_cds_feat) if is_forward_stranded else is_after(
                    five_prime_cds_feat)
            utr_filter_3p = is_after(
                three_prime_cds_feat) if is_forward_stranded else is_before(
                    three_prime_cds_feat)

            five_prime_utr_feats = filter(utr_filter_5p, utr_feats)
            three_prime_utr_feats = filter(utr_filter_3p, utr_feats)

            transcript_record = _VariantTreeTranscriptRecord(
                feat=transcript_feat,
                five_prime_cds_feat=five_prime_cds_feat,
                three_prime_cds_feat=three_prime_cds_feat,
                five_prime_utr_feats=five_prime_utr_feats,
                three_prime_utr_feats=three_prime_utr_feats,
                translated_feats=features["start_codon"] + cds_feats +
                features["stop_codon"])

            transcript_records[transcript_id] = transcript_record

            for cds_feat in cds_feats:
                cds_records.append(
                    _VariantTreeCodingRecord(transcript=transcript_record,
                                             feat=cds_feat))

        self.transcript_records = transcript_records
        self.tree = GenomeIntervalTree(lambda record: record.feat.pos,
                                       cds_records)

    # this method could be named more concisely maybe
    def _splice_codon_padding(self, seq, transcript_record, offset):
        transcript_coding_seq = transcript_record.seq[
            transcript_record.cds_slice]
        padded_transcript_coding_seq = (
            'N' * -(transcript_record.five_prime_cds.phase % -3)
        ) + transcript_coding_seq

        if len(padded_transcript_coding_seq) % 3 != 0:
            # FIXME(remove): temp trap
            print("Bad padded tx seq len")

        # Trim the first `offset` characters if `offset` is negative.
        trimmed_start = max(0, -offset)
        trimmed_end = min(
            len(seq),
            # Length of the transcript sequence relative to the sequence
            len(padded_transcript_coding_seq) - offset)

        splice_start = offset + trimmed_start
        splice_end = offset + trimmed_end

        nearest_codon_start = splice_start - (splice_start % 3)
        nearest_codon_end = splice_end - (splice_end % -3)

        # NOTE: For sequences which start before or end after the transcript
        # sequence, the nearest codon start/end indices may be greater than the
        # splice start/end indices, which is fine. In these cases, nothing will
        # be prepended/appended.

        splice_prepend_seq = padded_transcript_coding_seq[nearest_codon_start:
                                                          splice_start]
        splice_append_seq = padded_transcript_coding_seq[splice_end:
                                                         nearest_codon_end]

        return splice_prepend_seq + seq + splice_append_seq

    def _get_translatable_locseq(self, transcript_record):
        transcript_coding_seq = transcript_record.seq[
            transcript_record.cds_slice]

        inverse_phase = -(transcript_record.five_prime_cds.phase % -3)

        if transcript_record.feat.is_reverse_stranded:
            inverse_phase = -inverse_phase

        padded_transcript_coding_seq = (
            'N' * abs(inverse_phase)) + transcript_coding_seq

        pos = transcript_record.five_prime_cds.pos
        pos = pos.shifted_by(-inverse_phase, 0)

        return _LocalizedSeq(seq=padded_transcript_coding_seq, pos=pos)

    def _slice_seq_codons(self, seq, seq_slice):
        nearest_codon_start = seq_slice.start - (seq_slice.start % 3)
        nearest_codon_end = seq_slice.stop - (seq_slice.stop % -3)

        return seq[nearest_codon_start:nearest_codon_end]

    # -> Get coding sequence (including stop codon)
    # -> If 5' CDS has a nonzero phase, pad Ns in new seq
    #   -> What about 3'?
    # -> Replace REF in translatable sequence with ALT
    # -> Translate ref and alt sequences

    # Alternative procedure...
    # -> Get entire transcript sequence
    # -> Make copy of tx seq with REF->ALT
    # -> Splice together ref seq from CDS+stop codon bounds
    # -> Splice together alt seq from CDS+stop codon bounds shifted by len(ALT) - len(REF) after POS
    # -> Copy alt and ref seq and add phase padding to 5' (maybe 3'??)
    # -> Translate alt, ref
    # -> Add gaps to ref to maintain alignment (or not)
    # -> If n' not incomplete & start/stop codon missing, try to take orig. non-padded seq and splice in n' UTR until respective codon (for 5' & 3')
    # -> Re-translate alt, add gaps to ref

    def _splice_seq(self, seq, intervals):
        ordered_intervals = _merge_and_sort_intersecting_intervals(intervals)

        # Create empty seq with same alphabet
        spliced_seq = seq[0:0]

        for interval in ordered_intervals:
            spliced_seq += seq[slice(*interval)]

        return spliced_seq

    def predict_for_vcf_record(self, vcf_record):
        variant_results = []

        record_pos = GenomePosition.from_vcf_record_pos(vcf_record)

        ref = vcf_record.REF
        for alt in vcf_record.ALT:
            # Create a GenomePosition representing the range affected by the ALT
            # sequence.
            affected_pos = record_pos.shifted_by(
                vcf_alt_affected_range(ref, alt))

            # TODO: Would be better in theory if overlaps were used rather than
            # containments, but containments allow a very useful set of
            # assumptions to be made.
            coding_overlaps = self.tree.get_all_overlaps(affected_pos)
            transcript_ids = set(
                record.transcript.feat.attributes["transcript_id"]
                for record in coding_overlaps)

            for tx_id in transcript_ids:
                transcript = self.transcript_records[tx_id]
                tx_pos = transcript.feat.pos
                ref_tx_seq = Seq(self.genome_fasta["chr" + tx_pos.chrom]
                                 [tx_pos.start:tx_pos.end].seq,
                                 alphabet=Alphabet.generic_dna)

                ref_pos = record_pos.shifted_by(0, len(ref) - 1)
                ref_slice = ref_pos.slice_within(tx_pos)

                alt_tx_seq = ref_tx_seq[:ref_slice.start] \
                           + alt.value \
                           + ref_tx_seq[ref_slice.stop:]

                alt_tx_seq_len_delta = len(alt_tx_seq) - len(ref_tx_seq)

                # i think this is right, maybe could use some comments explaining logic
                # kinda hard to explain though

                ref_splice_intervals, alt_splice_intervals = tee(
                    (feat.pos.slice_within(tx_pos).indices(len(tx_pos))[:2]
                     for feat in transcript.translated_feats), 2)

                record_pos_tx_slice = record_pos.slice_within(tx_pos)

                alt_splice_intervals = (
                    (ival[0] + (0 if record_pos_tx_slice.start >= ival[0] else
                                alt_tx_seq_len_delta),
                     ival[1] + (0 if record_pos_tx_slice.start >= ival[1] else
                                alt_tx_seq_len_delta))
                    for ival in alt_splice_intervals)

                ref_splice_intervals = list(ref_splice_intervals)
                alt_splice_intervals = list(alt_splice_intervals)

                print(
                    list(
                        _merge_and_sort_intersecting_intervals(
                            ref_splice_intervals)))
                print(
                    list(
                        _merge_and_sort_intersecting_intervals(
                            alt_splice_intervals)))

                ref_coding_seq = self._splice_seq(
                    ref_tx_seq, ref_splice_intervals).reverse_complement()
                alt_coding_seq = self._splice_seq(
                    alt_tx_seq, alt_splice_intervals).reverse_complement()

                # FIXME: phase/padding + strandedness

                # coding_seq, cds_pos = self._get_translatable_locseq(transcript)

                # ref_pos = record_pos.shifted_by(0, len(ref) - 1)
                # ref_slice = ref_pos.slice_within(cds_pos)

                # alt_seq_padding = -(len(alt.value) % -3)

                # alt_coding_seq = coding_seq[:ref_slice.start] \
                #                + alt.value \
                #                + coding_seq[ref_slice.stop:] \
                #                + ('N' * alt_seq_padding)

                # aa_seq = coding_seq.translate()
                # alt_aa_seq = alt_coding_seq.translate()

                # psvm = ProteinSequenceVariantMatcher(alt.type, aa_seq,
                #                                      alt_aa_seq)
                #
                # posedit = psvm.create_hgvs_posedit()

                # transcript_id = transcript.feat.attributes["transcript_id"]
                # # FIXME: Should use protein product id, not transcript id
                # seqvar = hgvs.sequencevariant.SequenceVariant(ac=transcript_id,
                #                                               type='p',
                #                                               posedit=posedit)

                ref_aa_seq = ref_coding_seq.translate()
                alt_aa_seq = alt_coding_seq.translate()
                protein_id = transcript.feat.attributes["protein_id"]

                print(
                    f"ref {len(ref_tx_seq)}\t->\t{len(ref_coding_seq)}\t->\t{len(ref_aa_seq)}"
                )
                print(
                    f"alt {len(alt_tx_seq)}\t->\t{len(alt_coding_seq)}\t->\t{len(alt_aa_seq)}"
                )

                # Alt transcript data flags

                is_frameshift = (len(ref_coding_seq) -
                                 len(alt_coding_seq)) % 3 != 0
                variant_start_aa = None
                is_substitution = None
                is_ambiguous = alt_aa_seq.count('*') > 1

                diff_count = 0
                for pos, (ref_aa,
                          alt_aa) in enumerate(zip(ref_aa_seq, alt_aa_seq)):
                    if ref_aa != alt_aa:
                        diff_count += 1

                        if variant_start_aa is None:
                            variant_start_aa = pos + 1  # 1-indexed

                is_substitution = (
                    len(ref_aa_seq) == len(alt_aa_seq)) and diff_count == 1

                ref_tx_data = _MockRefTranscriptData(
                    transcript_sequence=str(ref_coding_seq),
                    aa_sequence=str(ref_aa_seq),
                    cds_start=None,  # Used for initialization
                    cds_stop=None,  # Used for initialization
                    protein_accession=protein_id)
                alt_tx_data = _MockAltTranscriptData(
                    transcript_sequence=str(alt_coding_seq),
                    aa_sequence=str(alt_aa_seq),
                    cds_start=None,  # Used for initialization
                    cds_stop=None,  # Used for initialization
                    protein_accession=protein_id,
                    is_frameshift=is_frameshift,
                    variant_start_aa=variant_start_aa,  # 1-indexed
                    frameshift_start=None,  # Seems to be unused currently
                    is_substitution=is_substitution,
                    # The `AltTranscriptData` builder used by hgvs` sets
                    # `is_ambiguous` to true if there are multiple stop codons
                    # in the AA sequence. This doesn't work in this case,
                    # because frameshifts typically create many new stop
                    # codons. It is uncertain why `hgvs` doesn't run into this
                    # issue; perhaps it's because they indiscriminately pad
                    # sequences that can't be cleanly translated.
                    is_ambiguous=(not is_frameshift) and is_ambiguous)

                protein_variant_builder = AltSeqToHgvsp(
                    ref_tx_data, alt_tx_data)
                protein_variant = protein_variant_builder.build_hgvsp()

                print("[[[ === MATCH ===")
                print(f"<<TX::{protein_id}>>")
                print(f"{ref_coding_seq} => {ref_aa_seq}")
                print("->")
                print(f"{alt_coding_seq} => {alt_aa_seq}")
                print()
                print(protein_variant)
                print("]]]")

                variant_results.append(
                    ProteinVariantResult(query_variant=alt,
                                         predicted_variant=protein_variant,
                                         transcript_feat=transcript))

        return variant_results

    enst_id_pattern = re.compile(r"ENST[\d\.]+")

    def parse_enst_id(self, enst_str):
        match = self.enst_id_pattern.match(enst_str)

        if not match:
            return None

        return match[0]

    cds_range_pattern = re.compile(r"CDS:(\d+)-(\d+)")

    def parse_cds_range(self, cds_str):
        match = self.cds_range_pattern

        if not match:
            return None

        return range(int(match[1]) - 1, int(match[2]))

    # Just jotting down some ideas

    # Main tree generation procedure
    #   -> iterate through transcript records
    #       -> try to grab ENST* transcript ID
    #       -> look up in annotated genome
    #       -> grab CDS entries
    #       -> add record to tree once for each CDS position
    #           -> might have to trim transcript for CDS region
    #       -> also keep track of gene name

    # Check for variant procedure
    #   -> look up variant position in main tree
    #   -> iterate through overlap records
    #       -> compose HGVS StructuralVariant(?)
    #   -> report found struct vars