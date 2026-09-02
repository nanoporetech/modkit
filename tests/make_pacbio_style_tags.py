"""Create tests/resources/pacbio_style_tags{,_all_tagged}.bam from
tests/resources/bc_anchored_10_reads.sorted.bam.

The records mimic PacBio HiFi reads processed with Jasmine >= 26.1.3:
`A+a.;C+h?;C+m?;T-a.` where

  * the sparse `C+h?` track comes from a model that is independent of the
    `C+m?` model, so at a position the two probabilities can sum to more than
    one (here the first C has m=255 and h=178),
  * 6mA is called on both strands (`A+a.` and `T-a.`, implicit mode).

In `pacbio_style_tags.bam` the last record has its MM/ML tags removed,
`pacbio_style_tags_no_untagged.bam` is the same file without that record.

Requires pysam. Run from the repository root:
    python tests/make_pacbio_style_tags.py
"""
import array

import pysam

SRC = "tests/resources/bc_anchored_10_reads.sorted.bam"
OUT = "tests/resources/pacbio_style_tags.bam"
OUT_NO_UNTAGGED = "tests/resources/pacbio_style_tags_no_untagged.bam"


def parse_mm(mm, ml):
    subs = [s for s in mm.rstrip(";").split(";") if s]
    parsed = {}
    i = 0
    for s in subs:
        parts = s.split(",")
        head, deltas = parts[0], parts[1:]
        parsed[head] = (deltas, list(ml[i:i + len(deltas)]))
        i += len(deltas)
    assert i == len(ml)
    return parsed


def transform(rec):
    parsed = parse_mm(rec.get_tag("MM"), rec.get_tag("ML"))
    h_deltas, _h_ml = parsed["C+h?"]
    m_deltas, m_ml = parsed["C+m?"]
    assert h_deltas == m_deltas
    # confident 5mC call on the first C ...
    m_ml[0] = 255
    # ... and an independent, confident 5hmC call at the same position
    h_deltas, h_ml = [h_deltas[0]], [178]
    seq = rec.get_forward_sequence()
    tracks = []
    if "A" in seq:
        tracks.append(("A+a.", ["0"], [10]))
    tracks.append(("C+h?", h_deltas, h_ml))
    tracks.append(("C+m?", m_deltas, m_ml))
    if "T" in seq:
        tracks.append(("T-a.", ["0"], [10]))
    new_mm = ";".join(f"{head}," + ",".join(d) for head, d, _ in tracks) + ";"
    new_ml = [x for _, _, ml in tracks for x in ml]
    rec.set_tag("MM", new_mm, "Z")
    rec.set_tag("ML", array.array("B", new_ml))
    return rec


if __name__ == "__main__":
    src = pysam.AlignmentFile(SRC)
    records = [transform(rec) for rec in src]
    n = len(records)
    with pysam.AlignmentFile(OUT_NO_UNTAGGED, "wb", template=src) as fh:
        for rec in records[:-1]:
            fh.write(rec)
    with pysam.AlignmentFile(OUT, "wb", template=src) as fh:
        for i, rec in enumerate(records):
            if i == n - 1:
                rec.set_tag("MM", None)
                rec.set_tag("ML", None)
            fh.write(rec)
    pysam.index(OUT)
    pysam.index(OUT_NO_UNTAGGED)
    print(f"wrote {OUT} ({n} records) and {OUT_NO_UNTAGGED} ({n - 1} records)")
