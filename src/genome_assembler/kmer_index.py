
import math

def kmer_generation(read, k=5):

    if k <= 0:
        raise ValueError("k must be positive")

    for i in range(len(read) - k + 1):

        kmer = read[i:i + k]

        if "N" not in kmer:
            yield i, kmer

def build_kmer_data(reads, k=5):

    frequency = {}
    index = {}

    for read_id, read in enumerate(reads):

        for position, kmer in kmer_generation(read, k):

            frequency[kmer] = frequency.get(kmer, 0) + 1

            if kmer not in index:
                index[kmer] = []

            index[kmer].append((read_id, position))

    return frequency, index

def filter_kmer_frequency(frequency, alpha=0.1):

    if not frequency:
        return frequency, {}

    high_freq_kmer = {}

    max_frequency = max(frequency.values())
    threshold = alpha * max_frequency

    for kmer in list(frequency):

        if frequency[kmer] > threshold:
            high_freq_kmer[kmer] = frequency.pop(kmer)

    return frequency, high_freq_kmer

def filter_kmer_index(index, frequency, c=2):

    high_frequency_index = {}

    threshold = kmer_standard_deviation(frequency, c)

    for kmer in list(index):

        if len(index[kmer]) > threshold:
            high_frequency_index[kmer] = index.pop(kmer)

    return index, high_frequency_index

def kmer_standard_deviation(frequency, c=2):

    if not frequency:
        return 0

    values = frequency.values()

    avg_frequency = sum(values) / len(frequency)

    variance = sum(
        (value - avg_frequency) ** 2
        for value in values
    ) / len(frequency)

    std_dev = math.sqrt(variance)

    return avg_frequency + c * std_dev

def build_kmer_index(reads, k=21, alpha=0.1, c=2):

    frequency, index = build_kmer_data(reads, k)

    normal_frequency, high_freq_kmers = (
        filter_kmer_frequency(frequency, alpha)
    )

    normal_index, high_frequency_index = (
        filter_kmer_index(index, normal_frequency, c)
    )

    return (
        normal_index,
        high_frequency_index,
        normal_frequency,
        high_freq_kmers
    )
