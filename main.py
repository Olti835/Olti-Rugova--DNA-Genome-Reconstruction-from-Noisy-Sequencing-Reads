from genome_assembler import read_fastq, build_kmer_index, kmers_frequency, filter_kmer_frequency, filter_kmer_index, candidate_key_generation
reads = read_fastq("data/raw/SRR29154704.fastq.gz")

print("Number of reads:", len(reads))
print("First 3 reads:", reads[:3])

index = build_kmer_index(reads[:20])

print("Number of unique k-mers:", len(index))
print("List of indexes:", index.items())

frequency = kmers_frequency(reads, 5)
filtered_frequency = filter_kmer_frequency(frequency, 0.1)

filter_index = filter_kmer_index(index, filtered_frequency, 2)

print("Number of filtered k-mers:", len(filter_index))
print("Filtered list of indexes:" , filter_index)

candidate_keys = candidate_key_generation(filter_index)

print("Candidate keys:", candidate_keys)
