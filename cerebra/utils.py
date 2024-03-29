import re
import numpy as np

from ncls import NCLS

class GenomePosition():
	genome_pos_pattern = re.compile(r"(.+):(\d+)-(\d+)")

	def __init__(self, chrom, start, end):
		self.chrom = chrom
		self.start = start
		self.end = end

	@classmethod
	def from_str(cls, pos_str):
		match = cls.genome_pos_pattern.match(pos_str)

		if not match:
			return None

		return cls(match[1], int(match[2]) - 1, int(match[3]))

	@classmethod
	def from_vcf_record(cls, record):
		CHROM = record.CHROM.replace("chr", "")
		return cls(CHROM, record.affected_start, record.affected_end)

	@classmethod
	def from_gtf_record(cls, record):
		return cls(record[0].replace("chr", ""), int(record[3]) - 1, int(record[4]))

	def contains(self, other):
		return other.chrom == self.chrom and other.start >= self.start and other.end <= self.end

	def __eq__(self, other):
		return self.chrom == other.chrom and self.start == other.start and self.end == other.end

	def __repr__(self):
		return "%s:%d-%d" % (self.chrom, self.start + 1, self.end)

	def __str__(self):
		return "%s:%d-%d" % (self.chrom, self.start + 1, self.end)

	def __len__(self):
		return self.end - self.start


class GenomeIntervalTree():
	def __init__(self, predicate, records):
		self.predicate = predicate
		self.records = []

		working_tree_map = {}

		for idx, record in enumerate(records):
			genome_pos = predicate(record)

			if genome_pos is None:
				continue

			chrom = genome_pos.chrom

			if not chrom in working_tree_map:
				# (starts, ends, ids)
				working_tree_map[chrom] = ([], [], [])

			starts, ends, ids = working_tree_map[chrom]
			starts.append(genome_pos.start)
			ends.append(genome_pos.end)
			ids.append(idx)

			self.records.append(record)

		tree_map = {}

		for chrom, (starts, ends, ids) in working_tree_map.items():
			tree_map[chrom] = NCLS(
				np.array(starts, dtype=np.long),
				np.array(ends, dtype=np.long),
				np.array(ids, dtype=np.long)
			)

		self.tree_map = tree_map

	def __reduce__(self):
		return (self.__class__, (self.predicate, self.records))

	def _make_query_params(self, genome_pos_list):
		starts = np.array([genome_pos.start for genome_pos in genome_pos_list])
		ends = np.array([genome_pos.end for genome_pos in genome_pos_list])
		ids = np.array(list(range(len(genome_pos_list))))

		return (starts, ends, ids)

	def _pick_best_record(self, from_ids=None, for_pos=None):
		if len(from_ids) < 1:
			return None

		if len(from_ids) == 1:
			return self.records[from_ids[0]]

		records = [self.records[record_id] for record_id in from_ids]

		scored_records = [(record, self._compute_jaccard_index(for_pos, self.predicate(record))) for record in records]
		sorted_records = sorted(scored_records, key=lambda tup: tup[1], reverse=True)

		return sorted_records[0][0]

	def _compute_jaccard_index(self, pos_a, pos_b):
		if pos_a.chrom != pos_b.chrom:
			return 0

		range_a = range(pos_a.start, pos_a.end)
		range_b = range(pos_b.start, pos_b.end)

		if range_b.start > range_a.stop or range_a.start > range_b.stop:
			return 0

		intersection = range(max(range_a.start, range_b.start), min(range_a.stop, range_b.stop))

		# The following is equivalent to |A ∩ B| / |A ∪ B|, but avoids computing
		# a union.
		# |A ∩ B| / (|A| + |B| + |A ∩ B|)
		return len(intersection) / (len(range_a) + len(range_b) - len(intersection))

	def has_overlap(self, genome_pos):
		tree = self.tree_map.get(genome_pos.chrom)

		if not tree:
			return False

		qparams = self._make_query_params([genome_pos])

		# has_overlap exists, but it is not vectorized. NCLS's vectorized
		# operations such as first_overlap_both are much faster.
		return bool(len((tree.first_overlap_both(*qparams)[0])))

	def get_first_overlap(self, genome_pos):
		tree = self.tree_map.get(genome_pos.chrom)

		if not tree:
			return None

		qparams = self._make_query_params([genome_pos])
		_, record_ids = tree.first_overlap_both(*qparams)

		if len(record_ids) < 1:
			return None

		return self.records[record_ids[0]]

	def get_best_overlap(self, genome_pos):
		tree = self.tree_map.get(genome_pos.chrom)

		if not tree:
			return None

		qparams = self._make_query_params([genome_pos])
		_, record_ids = tree.all_overlaps_both(*qparams)

		return self._pick_best_record(from_ids=record_ids, for_pos=genome_pos)

	def get_all_overlaps(self, genome_pos):
		tree = self.tree_map.get(genome_pos.chrom)

		if not tree:
			return []

		qparams = self._make_query_params([genome_pos])
		_, record_ids = tree.all_overlaps_both(*qparams)

		return [self.records[record_id] for record_id in record_ids]

	def get_first_containment(self, genome_pos):
		tree = self.tree_map.get(genome_pos.chrom)

		if not tree:
			return None

		qparams = self._make_query_params([genome_pos])
		_, record_ids = tree.all_containments_both(*qparams)

		if len(record_ids) < 1:
			return None

		return self.records[record_ids[0]]

	def get_best_containment(self, genome_pos):
		tree = self.tree_map.get(genome_pos.chrom)

		if not tree:
			return None

		qparams = self._make_query_params([genome_pos])
		_, record_ids = tree.all_containments_both(*qparams)

		return self._pick_best_record(from_ids=record_ids, for_pos=genome_pos)

	def get_all_containments(self, genome_pos):
		tree = self.tree_map.get(genome_pos.chrom)

		if not tree:
			return []

		qparams = self._make_query_params([genome_pos])
		_, record_ids = tree.all_containments_both(*qparams)

		return [self.records[record_id] for record_id in record_ids]
