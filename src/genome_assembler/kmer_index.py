import math


def kmer_generation(read, k=5):
    kmers = []

    for i in range(len(read) - k + 1):
        kmer = read[i:i + k]

        if "N" not in kmer:
            kmers.append(kmer)

    return kmers


def kmers_frequency(reads, k):
    frequency = {}

    for read in reads:
        kmers = kmer_generation(read, k)

        for kmer in kmers:
            if kmer not in frequency:
                frequency[kmer] = 1
            else:
                frequency[kmer] +=1


    return frequency

def build_kmer_index(reads, k=5):
    index = {}

    for read_id, read in enumerate(reads):
        kmers = kmer_generation(read, k)

        for position, kmer in enumerate(kmers):
            if kmer not in index:
                index[kmer] = []

            index[kmer].append((read_id, position))

    return index


def filter_kmer_frequency(frequency, alpha=0.1):
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

def kmer_standard_deviation(frequency, c):
    unique_kmers_rec = 1 / len(frequency)
    frequency_sum = 0

    for kmer in frequency:
        frequency_sum += frequency[kmer]

    avg_frequency = unique_kmers_rec * frequency_sum

    sdfm_sum = 0  # Squared deviation from the mean

    for kmer in frequency:
        sdfm_sum += math.pow(frequency[kmer] - avg_frequency, 2)

    variance = unique_kmers_rec * sdfm_sum
    std_dev = math.sqrt(variance)

    return avg_frequency + (c * std_dev)


