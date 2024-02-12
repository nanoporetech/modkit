import os
import pysam

test_dir = os.path.dirname(__file__)


seq_lengths = dict()
with pysam.AlignmentFile(os.path.join(test_dir, "resources", "donor_read_sort.bam"), "rb") as in_bam:
    with pysam.AlignmentFile(os.path.join(test_dir, "resources", "donor_read_sort_mn_tag.bam"), "wb", template=in_bam) as out_bam:
        for record in in_bam:
            seq_len = len(record.query_sequence)
            read_id = record.query_name
            assert not record.has_tag("MN")
            seq_lengths[read_id] = seq_len
            record.set_tag("MN", value=seq_len, value_type="i")
            out_bam.write(record)

print("> added MN to original records")

with pysam.AlignmentFile(os.path.join(test_dir, "resources", "trimmed_read_sort.mapped.bam"), "rb") as in_bam:
    with pysam.AlignmentFile(os.path.join(test_dir, "resources", "trimmed_read_sort_mn_tag.mapped.bam"), "wb", template=in_bam) as out_bam:
        for record in in_bam:
            read_id = record.query_name
            # add the original/incorrect MN tag
            seq_len = seq_lengths[read_id]
            assert not record.has_tag("MN")
            record.set_tag("MN", value=seq_len, value_type="i")
            out_bam.write(record)

print("> added MN to trimmed records")

print("> finished")

