import random

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Perform global alignment using Needleman-Wunsch algorithm.

    Parameters:
        seq1: str - First sequence
        seq2: str - Second sequence
        match: int - Score for a match (default: 1)
        mismatch: int - Score for a mismatch (default: -1)
        gap: int - Penalty for a gap (default: -2)

    Returns:
        aligned_seq1: str - Aligned version of seq1
        aligned_seq2: str - Aligned version of seq2
        alignment_score: int - Final alignment score
        traceback_path: list - Traceback path used for alignment
        scoring_matrix: list - Scoring matrix used during alignment
    """
    n = len(seq1) + 1
    m = len(seq2) + 1

    # Full scoring matrix for accurate traceback
    scoring_matrix = [[0] * m for _ in range(n)]
    traceback_matrix = [[None] * m for _ in range(n)]

    # Initialize first row and column with gap penalties
    for i in range(1, n):
        scoring_matrix[i][0] = i * gap
        traceback_matrix[i][0] = 'U'
    for j in range(1, m):
        scoring_matrix[0][j] = j * gap
        traceback_matrix[0][j] = 'L'

    # Fill scoring and traceback matrices
    for i in range(1, n):
        for j in range(1, m):
            match_mismatch = scoring_matrix[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            insert_gap_seq1 = scoring_matrix[i - 1][j] + gap
            insert_gap_seq2 = scoring_matrix[i][j - 1] + gap

            scoring_matrix[i][j] = max(match_mismatch, insert_gap_seq1, insert_gap_seq2)

            if scoring_matrix[i][j] == match_mismatch:
                traceback_matrix[i][j] = 'D'
            elif scoring_matrix[i][j] == insert_gap_seq1:
                traceback_matrix[i][j] = 'U'
            else:
                traceback_matrix[i][j] = 'L'

    # Traceback to construct alignment
    aligned_seq1 = []
    aligned_seq2 = []
    traceback_path = []
    i, j = len(seq1), len(seq2)

    while i > 0 or j > 0:
        traceback_path.append((i, j))
        if traceback_matrix[i][j] == 'D':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 'U':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        elif traceback_matrix[i][j] == 'L':
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    aligned_seq1 = ''.join(reversed(aligned_seq1))
    aligned_seq2 = ''.join(reversed(aligned_seq2))
    traceback_path = traceback_path[::-1]

    alignment_score = scoring_matrix[len(seq1)][len(seq2)]

    return aligned_seq1, aligned_seq2, alignment_score, traceback_path, scoring_matrix

def calculate_alignment_statistics(aligned_seq1, aligned_seq2):
    """
    Calculate statistics for the alignment.

    Parameters:
        aligned_seq1: str - Aligned first sequence
        aligned_seq2: str - Aligned second sequence

    Returns:
        dict - Statistics including matches, mismatches, gaps, and percentages.
    """
    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
    mismatches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a != b and a != '-' and b != '-')
    gaps = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == '-' or b == '-')
    total = len(aligned_seq1)
    return {
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "match_percentage": (matches / total) * 100,
        "gap_percentage": (gaps / total) * 100,
    }

def generate_random_sequence(length):
    """
    Generate a random DNA sequence of given length.

    Parameters:
        length: int - Length of the DNA sequence to generate

    Returns:
        str - Random DNA sequence
    """
    return ''.join(random.choice("ACGT") for _ in range(length))

def display_scoring_matrix(scoring_matrix):
    """
    Display the scoring matrix in a readable format.

    Parameters:
        scoring_matrix: list - The scoring matrix to display
    """
    print("\nScoring Matrix:")
    for row in scoring_matrix:
        print("\t".join(map(str, row)))

# Example Usage
if __name__ == "__main__":
    choice = input("Would you like to input sequences (1) or generate random sequences (2)? ").strip()

    if choice == "1":
        seq1 = input("Enter the first HLA-DQ gene sequence: ").strip().upper()
        seq2 = input("Enter the second HLA-DQ gene sequence: ").strip().upper()
    elif choice == "2":
        length1 = int(input("Enter the length of the first random sequence: "))
        length2 = int(input("Enter the length of the second random sequence: "))
        seq1 = generate_random_sequence(length1)
        seq2 = generate_random_sequence(length2)
        print(f"Generated Sequences:\nSequence 1: {seq1}\nSequence 2: {seq2}")
    else:
        print("Invalid choice. Exiting.")
        exit()

    aligned_seq1, aligned_seq2, score, traceback_path, scoring_matrix = needleman_wunsch(seq1, seq2)

    print("\nAlignment:")
    print(aligned_seq1)
    print("".join(["|" if a == b else " " for a, b in zip(aligned_seq1, aligned_seq2)]))
    print(aligned_seq2)
    print(f"\nAlignment Score: {score}")

    stats = calculate_alignment_statistics(aligned_seq1, aligned_seq2)
    print("\nAlignment Statistics:")
    for key, value in stats.items():
        print(f"{key.capitalize()}: {value}")

    print("\nTraceback Path:")
    print(traceback_path)

    display_scoring_matrix(scoring_matrix)

    with open("alignment_results.txt", "w") as f:
        f.write(f"Alignment:\n{aligned_seq1}\n")
        f.write("".join(["|" if a == b else " " for a, b in zip(aligned_seq1, aligned_seq2)]) + "\n")
        f.write(f"{aligned_seq2}\n")
        f.write(f"Score: {score}\n")
        f.write("\nAlignment Statistics:\n")
        for key, value in stats.items():
            f.write(f"{key.capitalize()}: {value}\n")
        f.write("\nTraceback Path:\n")
        f.write(str(traceback_path) + "\n")
        f.write("\nScoring Matrix:\n")
        for row in scoring_matrix:
            f.write("\t".join(map(str, row)) + "\n")
    print("\nResults have been exported to 'alignment_results.txt'")
