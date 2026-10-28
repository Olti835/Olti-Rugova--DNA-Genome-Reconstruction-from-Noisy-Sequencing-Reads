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









