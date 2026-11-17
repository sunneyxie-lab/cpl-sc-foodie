#!/usr/bin/env python3
"""Select the Perl Trim Galore 0.6.11 preset without considering BGI."""
import gzip
import sys

PROBES = {
    "illumina": "AGATCGGAAGAGC",
    "nextera": "CTGTCTCTTATA",
    "small_rna": "TGGAATTCTCGG",
}
CONTRIBUTING_RECORDS = 999_999


def select_adapter(counts):
    """Insertion order implements Illumina, then Nextera precedence on ties."""
    return max(PROBES, key=lambda preset: counts[preset])


def count_adapters(records):
    counts = dict.fromkeys(PROBES, 0)
    for index, sequence in enumerate(records):
        # Perl consumes the millionth record but does not count it.
        if index == CONTRIBUTING_RECORDS:
            break
        for preset, probe in PROBES.items():
            counts[preset] += probe in sequence
    return counts


def fastq_sequences(handle):
    while True:
        header = handle.readline()
        if not header:
            return
        sequence = handle.readline().rstrip("\r\n")
        separator = handle.readline()
        quality = handle.readline().rstrip("\r\n")
        if not header.startswith("@") or not separator.startswith("+") or len(sequence) != len(quality):
            raise ValueError("Malformed FASTQ record during adapter selection")
        yield sequence


def main():
    with gzip.open(sys.argv[1], "rt", encoding="ascii") as handle:
        counts = count_adapters(fastq_sequences(handle))
    selected = select_adapter(counts)
    print(f"Legacy adapter counts: {counts}; selected --{selected}", file=sys.stderr)
    print(f"--{selected}")


if __name__ == "__main__":
    main()
