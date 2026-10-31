def candidate_key_generation(index):
    candidate_keys = {}

    for occurrences in index.values():
        for i in range(len(occurrences)):
            for j in range(i + 1, len(occurrences)):
                pair = tuple(sorted((occurrences[i][0], occurrences[j][0])))
                first_pos, second_pos = occurrences[i][1], occurrences[j][1]

                diff = abs(first_pos - second_pos)

                if pair not in candidate_keys:
                    candidate_keys[pair] = {}

                if diff not in candidate_keys[pair]:
                    candidate_keys[pair][diff] = 1
                else:
                    candidate_keys[pair][diff] += 1

    return candidate_keys

def candidate_key_filter(candidate_keys, min_shared_kmers=5, min_position_consistency=0.7):
    filtered_candidates = {}
    for pair, distances_frequency in candidate_keys.items():
        freq_sum, pos_consistency = calculate_position_consistency(distances_frequency)

        if freq_sum >= min_shared_kmers and pos_consistency >= min_position_consistency:
            continue

        filtered_candidates[pair] = candidate_keys.pop(pair)

    return candidate_keys, filtered_candidates

def calculate_position_consistency(distances_frequency):
    try:
        frequency_sum = sum(distances_frequency.values())
        max_frequency = max(distances_frequency.values())

    except ValueError:
        print("Dictionary is empty")

    return frequency_sum, float(max_frequency/frequency_sum)










