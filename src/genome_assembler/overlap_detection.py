
#marrja dhe pergaditja e reads
def get_candidate_reads(pair, reads):
    read1_id, read2_id = pair

    read1 = reads[read1_id]
    read2 = reads[read2_id]

    return read1, read2


def reverse_complement(read):
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "N": "N"
    }

    result = ""

    for nucleotide in read:
        result = complement[nucleotide] + result

    return result


#semi-global alignment
def calculate_alignment_score(nucleotide1, nucleotide2):
    if nucleotide1 == "-" or nucleotide2 == "-":
        return -2
    elif nucleotide1 == nucleotide2:
        return 2

    return -1


def semi_global_alignment(read1, read2):
    gap_penalty = 2
    m = len(read1)
    n = len(read2)

    matrix = initialize_alignment_matrix(read1, read2, gap_penalty)

    for i in range(1, m+1):
        for j in range(1, n+1):

            score = calculate_alignment_score(read1[i-1], read2[j-1])

            diagonal = matrix[i - 1][j - 1] + score
            up = matrix[i - 1][j] - gap_penalty
            left = matrix[i][j - 1] - gap_penalty

            matrix[i][j] = max(diagonal, up, left)

    return matrix


def initialize_alignment_matrix(read1, read2, gap_penalty=2):
    m = len(read1)
    n = len(read2)

    matrix = [[0 for j in range(n + 1)] for i in range(m + 1)]

    for i in range(1, m + 1):
        matrix[i][0] = -i * gap_penalty

    for j in range(n + 1):
        matrix[0][j] = 0

    return matrix


def traceback_alignment(read1, read2, matrix, gap_penalty=2):
    i = len(read1)
    j = max(range(len(read2) + 1), key=lambda x: matrix[i][x])

    aligned1 = ""
    aligned2 = ""

    while i > 0 and j > 0:

        score = calculate_alignment_score(read1[i-1], read2[j-1])

        # Diagonal: Match or Mismatch
        if matrix[i][j] == matrix[i-1][j-1] + score:
            aligned1 = read1[i-1] + aligned1
            aligned2 = read2[j-1] + aligned2
            i -= 1
            j -= 1

        # Up: Deletion
        elif matrix[i][j] == matrix[i-1][j] - gap_penalty:
            aligned1 = read1[i-1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1

        # Left: Insertion
        elif matrix[i][j] == matrix[i][j-1] - gap_penalty:
            aligned1 = "-" + aligned1
            aligned2 = read2[j-1] + aligned2
            j -= 1

    return aligned1, aligned2


#vleresimi i overlapit
def calculate_overlap_length(aligned1, aligned2):
    overlap_length = 0

    for i in range(len(aligned1)):
        if aligned1[i] != "-" and aligned2[i] != "-":
            overlap_length += 1

    return overlap_length


def calculate_sequence_identity(aligned1, aligned2):
    matches = 0
    overlap_length = calculate_overlap_length(aligned1, aligned2)

    for i in range(len(aligned1)):
        if (aligned1[i] != "-" and aligned2[i] != "-" and
            aligned1[i] == aligned2[i]):
            matches += 1

    return matches/overlap_length * 100


def calculate_error_rate(aligned1, aligned2):
    mismatches = 0
    gaps = 0
    alignment_len = len(aligned1)

    if alignment_len == 0:
        return 0

    for i in range(alignment_len):
        if aligned1[i] == "-" or aligned2[i] == "-":
            gaps += 1
        elif aligned1[i] != aligned2[i]:
            mismatches += 1

    return (mismatches + gaps) / alignment_len * 100


#pranimi ose refuzimi i overlap
#kontrollon:
#     min_overlap_length
#     min_sequence_identity
#     max_error_rate
#     min_alignment_score

def validate_overlap(
    overlap_length,
    sequence_identity,
    error_rate,
    alignment_score):

    min_overlap_length = 30
    min_sequence_identity = 90
    max_error_rate = 10
    min_alignment_score = 40

    if (overlap_length >= min_overlap_length
        and sequence_identity >= min_sequence_identity
        and error_rate <= max_error_rate
        and alignment_score >= min_alignment_score):

            return True

    return False


#funksioni kryesor

def approximate_overlap_detection(reads, candidate_keys):

    valid_overlaps = []
    rejected_overlaps = []

    for pair in candidate_keys:

        read1_id, read2_id = pair

        read1 = reads[read1_id]
        read2 = reads[read2_id]

        orientations = [
            ("forward", read2),
            ("reverse", reverse_complement(read2))
        ]

        best_overlap = None
        best_score = float("-inf")

        for orientation, oriented_read2 in orientations:

            matrix = semi_global_alignment(read1, oriented_read2)

            i = len(read1)
            j = 0
            alignment_score = matrix[i][0]

            for col in range(1, len(oriented_read2) + 1):
                if matrix[i][col] > alignment_score:
                    alignment_score = matrix[i][col]
                    j = col

            # 5. Traceback
            aligned1, aligned2 = traceback_alignment(
                read1, oriented_read2, matrix
            )

            # 6. Calculate metrics
            overlap_length = calculate_overlap_length(
                aligned1, aligned2
            )

            sequence_identity = calculate_sequence_identity(
                aligned1, aligned2
            )

            error_rate = calculate_error_rate(
                aligned1, aligned2
            )

            is_valid = validate_overlap(
                overlap_length,
                sequence_identity,
                error_rate,
                alignment_score
            )

            if alignment_score > best_score:
                best_score = alignment_score

                best_overlap = {
                    "read1_id": read1_id,
                    "read2_id": read2_id,
                    "orientation": orientation,
                    "aligned1": aligned1,
                    "aligned2": aligned2,
                    "overlap_length": overlap_length,
                    "sequence_identity": sequence_identity,
                    "error_rate": error_rate,
                    "alignment_score": alignment_score,
                    "valid": is_valid
                }

        if best_overlap["valid"]:
            valid_overlaps.append(best_overlap)
        else:
            rejected_overlaps.append(best_overlap)

    return valid_overlaps, rejected_overlaps
