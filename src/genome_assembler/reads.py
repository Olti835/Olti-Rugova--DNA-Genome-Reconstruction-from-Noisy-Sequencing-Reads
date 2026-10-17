import gzip
from pathlib import Path
from shlex import join


def read_fastq(file_path):
    reads = []
    full_read = ""

    with gzip.open(file_path, "rt") as file:
        while True:
            identificator = file.readline()

            if not identificator:
                break

            sequence = file.readline()
            separator = file.readline()
            quality = file.readline()

            if not separator.startswith("+"):
                raise ValueError("Invalid Separator")

            if len(sequence) != len(quality):
                raise ValueError("Sequence and quality lengths do not match")

            reads.append(_quality_validator(sequence, quality))

    return reads


def _quality_validator(sequence, quality, min_phred=30):
    filtered = []
    for base, char in zip(sequence, quality):
        phred = ord(char) - 33

        if phred < min_phred:
            filtered.append("N")
        else:
            filtered.append(base)

    return "".join(filtered)

