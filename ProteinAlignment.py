# Sam Gass
# scg0040
# Computational Biology
# Project 1
# python3 coded in Microsoft Visual Studio Code
# Last Modified: February 6, 2020
#
# This python code takes in two FASTA files representing protein
# sequences of amino acids and runs the Needleman-Wunsch algorithm
# for global alignment and the Watermann-Smith algorithm for local
# alignment on the sequences and outputs their alignment score and
# and a visualization of the protein sequence alignment. It uses
# the BLOSUM62 matrix as the objective function.


# A global dictionary of dictionaries that represents the BLOSUM62 matrix
gapPenalty = -4
blosum = {
    '*': {'*': 1, 'A': -4, 'R': -4, 'N': -4, 'D': -4, 'C': -4, 'Q': -4, 'E': -4, 'G': -4, 'H': -4, 'I': -4, 'L': -4, 'K': -4, 'M': -4, 'F': -4, 'P': -4, 'S': -4, 'T': -4, 'W': -4, 'Y': -4, 'V': -4, 'B': -4, 'Z': -4, 'X': -4},
    'A': {'*': -4, 'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0, 'B': -2, 'Z': -1, 'X': 0},
    'R': {'*': -4, 'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3, 'B': -1, 'Z': 0, 'X': -1},
    'N': {'*': -4, 'A': -2, 'R': 0, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -3, 'L': -3, 'K': 0, 'M': -2, 'F': -3, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3, 'B': 3, 'Z': 0, 'X': -1},
    'D': {'*': -4, 'A': -2, 'R': -2, 'N': 1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4, 'Z': 1, 'X': -1},
    'C': {'*': -4, 'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': 9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1, 'B': -3, 'Z': -3, 'X': -2},
    'Q': {'*': -4, 'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 1, 'M': 0, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2, 'B': 0, 'Z': 3, 'X': -1},
    'E': {'*': -4, 'A': -1, 'R': 0, 'N': 0, 'D': 2, 'C': -4, 'Q': 2, 'E': 5, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -2, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1, 'Z': 4, 'X': -1},
    'G': {'*': -4, 'A': 0, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3, 'B': -1, 'Z': -2, 'X': -1},
    'H': {'*': -4, 'A': -2, 'R': 0, 'N': 1, 'D': -1, 'C': -3, 'Q': 0, 'E': 0, 'G': -2, 'H': 8, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y': 2, 'V': -3, 'B': 0, 'Z': 0, 'X': -1},
    'I': {'*': -4, 'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': 4, 'L': 2, 'K': -3, 'M': 1, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V': 3, 'B': -3, 'Z': -3, 'X': -1},
    'L': {'*': -4, 'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 2, 'L': 4, 'K': -2, 'M': 2, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 1, 'B': -4, 'Z': -3, 'X': -1},
    'K': {'*': -4, 'A': -1, 'R': 2, 'N': 0, 'D': -1, 'C': -3, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K': 5, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'Z': 1, 'X': -1},
    'M': {'*': -4, 'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q': 0, 'E': -2, 'G': -3, 'H': -2, 'I': 1, 'L': 2, 'K': -1, 'M': 5, 'F': 0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1, 'B': -3, 'Z': -1, 'X': -1},
    'F': {'*': -4, 'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0, 'L': 0, 'K': -3, 'M': 0, 'F': 6, 'P': -4, 'S': -2, 'T': -2, 'W': 1, 'Y': 3, 'V': -1, 'B': -3, 'Z': -3, 'X': -1},
    'P': {'*': -4, 'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P': 7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2, 'B': -2, 'Z': -1, 'X': -2},
    'S': {'*': -4, 'A': 1, 'R': -1, 'N': 1, 'D': 0, 'C': -1, 'Q': 0, 'E': 0, 'G': 0, 'H': -1, 'I': -2, 'L': -2, 'K': 0, 'M': -1, 'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'Z': 0, 'X': 0},
    'T': {'*': -4, 'A': 0, 'R': -1, 'N': 0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0, 'B': -1, 'Z': -1, 'X': 0},
    'W': {'*': -4, 'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2, 'K': -3, 'M': -1, 'F': 1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y': 2, 'V': -3, 'B': -4, 'Z': -3, 'X': -2},
    'Y': {'*': -4, 'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F': 3, 'P': -3, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1, 'B': -3, 'Z': -2, 'X': -1},
    'V': {'*': -4, 'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': 3, 'L': 1, 'K': -2, 'M': 1, 'F': -1, 'P': -2, 'S': -2, 'T': 0, 'W': -3, 'Y': -1, 'V': 4, 'B': -3, 'Z': -2, 'X': -1},
    'B': {'*': -4, 'A': -2, 'R': -1, 'N': 3, 'D': 4, 'C': -3, 'Q': 0, 'E': 1, 'G': -1, 'H': 0, 'I': -3, 'L': -4, 'K': 0, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4, 'Z': 1, 'X': -1},
    'Z': {'*': -4, 'A': -1, 'R': 0, 'N': 0, 'D': 1, 'C': -3, 'Q': 3, 'E': 4, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1, 'Z': 4, 'X': -1},
    'X': {'*': -4, 'A': 0, 'R': -1, 'N': -1, 'D': -1, 'C': -2, 'Q': -1, 'E': -1, 'G': -1, 'H': -1, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -1, 'P': -2, 'S': 0, 'T': 0, 'W': -2, 'Y': -1, 'V': -1, 'B': -1, 'Z': -1, 'X': -1}
}

# takes in a string representing a FASTA file name/location, skips the
# first line, then puts every char in a list. Returns list of char.
# if file not found a string is returned


def fileToList(filename):
    seqList = ['*']    # list of char in file

    # open file and check for file not found
    try:
        file = open(filename, 'r')
    except FileNotFoundError:
        return "file not found"

    # skip first line, otherwise add all char to list except newline char
    for line in file:
        for c in line:
            if c == '>':
                break
            if c == '\n':
                continue
            seqList.append(c)

    file.close()  # close file

    # print(seqList)
    return seqList


def needlemanWusch(s1, s2):
    # instantiate an empty matrix
    matrix = []
    for i in range(len(s1)):
        empty = []
        matrix.append(empty)
        for j in range(len(s2)):
            matrix[i].append(0)

    # instantiate base cases of matrix
    # matrix[0][0] = blosum[s1[0]][s2[0]] is a gap vs a gap a one or zero?
    for i in range(1, len(s1)):
        matrix[i][0] = blosum[s1[i]][s2[0]] + matrix[i-1][0]
    for j in range(1, len(s2)):
        matrix[0][j] = blosum[s1[0]][s2[j]] + matrix[0][j-1]

    # dynamic programing! Fill the matrix
    for i in range(1, len(s1)):
        for j in range(1, len(s2)):
            diagonal = blosum[s1[i]][s2[j]] + matrix[i-1][j-1]
            up = matrix[i-1][j] + gapPenalty
            left = matrix[i][j-1] + gapPenalty
            # Choose path
            best = max(diagonal, up, left)
            matrix[i][j] = best

    # sequence to be output
    match1 = []
    match2 = []
    matches = []
    # backtracking loop
    i = len(s1) - 1
    j = len(s2) - 1
    score = matrix[i][j]
    while i != 0 and j != 0:
        diagonal = blosum[s1[i]][s2[j]] + matrix[i-1][j-1]
        up = matrix[i-1][j] + gapPenalty
        left = matrix[i][j-1] + gapPenalty
        # Choose path
        best = max(diagonal, up, left)
        if diagonal == best:
            match1.append(s1[i])
            match2.append(s2[j])
            if blosum[s1[i]][s2[j]] > 0:
                matches.append('|')
            else:
                matches.append('*')
            i = i - 1
            j = j - 1
        elif up == best:
            match1.append(s1[i])
            match2.append('-')
            matches.append(' ')
            i = i - 1
        else:
            match1.append('-')
            match2.append(s2[j])
            matches.append(' ')
            j = j - 1
    match1.reverse()
    match2.reverse()
    matches.reverse()

    # output. Prints header, score, then 80 characters at a time
    # of each sequence with a | if its a match, a * if a mismatch
    # space otherwise
    print('\nNeedleman-Wusch Global Alignment:')
    print('__________________________________\n')
    print("Score:", score, '\n')
    start = 0
    end = len(match1)
    count = 0
    while start < end:
        for i in range(start, end):
            print(match1[i], end='')
            count += 1
            if count >= 80:
                count = 0
                break
        print()
        count = 0
        for j in range(start, end):
            print(matches[j], end='')
            count += 1
            if count >= 80:
                count = 0
                break
        print()
        count = 0
        for k in range(start, end):
            print(match2[k], end='')
            count += 1
            if count >= 80:
                count = 0
                break
        start += 80
        count = 0
        print('\n')


def smithWaterman(s1, s2):
    # instantiate an empty matrix
    matrix = []
    for i in range(len(s1)):
        empty = []
        matrix.append(empty)
        for j in range(len(s2)):
            matrix[i].append(0)

    # dynamic programing! Fill the matrix
    for i in range(1, len(s1)):
        for j in range(1, len(s2)):
            diagonal = blosum[s1[i]][s2[j]] + matrix[i-1][j-1]
            up = matrix[i-1][j] + gapPenalty
            left = matrix[i][j-1] + gapPenalty
            # Choose path or 0
            best = max(diagonal, up, left, 0)
            matrix[i][j] = best

    # find max value
    highest = 0
    coordinates = [0, 0]
    for i in range(len(s1) - 1, -1, -1):
        for j in range(len(s2) - 1, -1, -1):
            if matrix[i][j] > highest:
                highest = matrix[i][j]
                coordinates[0] = i
                coordinates[1] = j

    # Trackback
    match1 = []
    match2 = []
    matches = []

    i = coordinates[0]
    j = coordinates[1]

    while matrix[i][j] != 0:
        diagonal = blosum[s1[i]][s2[j]] + matrix[i-1][j-1]
        up = matrix[i-1][j] + gapPenalty
        left = matrix[i][j-1] + gapPenalty
        # Choose path
        best = max(diagonal, up, left)
        if diagonal == best:
            match1.append(s1[i])
            match2.append(s2[j])
            if blosum[s1[i]][s2[j]] > 0:
                matches.append('|')
            else:
                matches.append('*')
            i = i - 1
            j = j - 1
        elif up == best:
            match1.append(s1[i])
            match2.append('-')
            matches.append(' ')
            i = i - 1
        else:
            match1.append('-')
            match2.append(s2[j])
            matches.append(' ')
            j = j - 1
    match1.reverse()
    match2.reverse()
    matches.reverse()

    # output. Prints header, score, then 80 characters at a time
    # of each sequence with a | if its a match, a * if a mismatch
    # space otherwise
    print('\nSmith-Waterman Local Alignment:')
    print('__________________________________\n')
    print("Score:", highest, '\n')
    start = 0
    end = len(match1)
    count = 0
    while start < end:
        for i in range(start, end):
            print(match1[i], end='')
            count += 1
            if count >= 80:
                count = 0
                break
        print()
        count = 0
        for j in range(start, end):
            print(matches[j], end='')
            count += 1
            if count >= 80:
                count = 0
                break
        print()
        count = 0
        for k in range(start, end):
            print(match2[k], end='')
            count += 1
            if count >= 80:
                count = 0
                break
        start += 80
        count = 0
        print('\n')


if __name__ == '__main__':
    # loop that runs the main interface exits when user says so
    while True:
        print("\nTo escape type exit")
        fn1 = input("FASTA file for first protein sequence: ")
        if fn1 == 'exit':
            break
        sequence1 = fileToList(fn1)
        if isinstance(sequence1, str):  # fileToList only returns str if fnf
            print('file not found')
            continue

        fn2 = input("FASTA file for second protein sequence: ")
        if fn2 == 'exit':
            break
        sequence2 = fileToList(fn2)
        if isinstance(sequence2, str):  # fileToList only returns str if fnf
            print('file not found')
            continue
        # Run algorithms on files
        needlemanWusch(sequence1, sequence2)
        smithWaterman(sequence1, sequence2)
    # end of program
