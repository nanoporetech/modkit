import pysam
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("in_bam")
    parser.add_argument("out_bam")
    parser.add_argument("--head-trim", type=int)
    parser.add_argument("--tail-trim", type=int)
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    headtrim = args.head_trim if args.head_trim is not None else 0
    with pysam.AlignmentFile(args.in_bam) as bam:
        with pysam.AlignmentFile(args.out_bam, "wb", template=bam) as out_bam:
            for record in bam:
                tailtrim = -args.tail_trim if args.tail_trim is not None else len(record.query_sequence)
                seq = record.query_sequence[headtrim:tailtrim]
                qual = record.query_qualities[headtrim:tailtrim]
                record.query_sequence = seq
                record.query_qualities = qual
                record.cigar = []
                out_bam.write(record)


    ...

