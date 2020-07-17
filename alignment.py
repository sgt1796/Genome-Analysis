# This file involves Needleman-Wunsch algorithm and Smith-Waterman algorithm
def NWalign(seq1, seq2, match=2, mismatch=-1, gap=-1):
    length1 = len(seq1)
    length2 = len(seq2)
    sMatrix, tbMatrix = [], []
    # initialize scoring matrix
    # 1. create matrix with all 0's
    for x in range(length2 + 1):
        sMatrix.append([0] * (length1 + 1))
        tbMatrix.append([0] * (length1 + 1))
    # 2. initialize matrix with initial scores
    for x in range(length1 + 1):
        sMatrix[0][x] = -x
        tbMatrix[0][x] = 'left'
    for x in range(length2 + 1):
        sMatrix[x][0] = -x
        tbMatrix[x][0] = 'up'
    tbMatrix[0][0] = 'done'

    # Now start scoring, and filling the traceback matrix

    for x in range(1, length2 + 1):
        for y in range(1, length1 + 1):
            check = mismatch
            if seq1[y - 1] == seq2[x - 1]:
                check = match
            score = max(sMatrix[x - 1][y - 1] + check, sMatrix[x - 1][y] + gap, sMatrix[x][y - 1] + gap, 0)
            # this block is for generating tbMatrix
            if score == sMatrix[x - 1][y - 1] + check:
                tbMatrix[x][y] = 'diag'
            elif score == sMatrix[x - 1][y] + gap:
                tbMatrix[x][y] = 'up'
            else:
                tbMatrix[x][y] = 'left'

            sMatrix[x][y] = score

    # now trace back the matrix to find the optimal alignment
    # this method will fix the end together
    pos1, pos2 = length2, length1
    hold = tbMatrix[pos1][pos2]
    output1, output2 = '', ''
    while hold != 'done':
        if hold == 'diag':
            output1 = seq1[pos2 - 1] + output1
            output2 = seq2[pos1 - 1] + output2
            pos1 -= 1
            pos2 -= 1
        elif hold == 'left':
            output1 = seq1[pos2 - 1] + output1
            output2 = '-' + output2
            pos2 -= 1
        else:
            output1 = '-' + output1
            output2 = seq2[pos1 - 1] + output2
            pos1 -= 1
        hold = tbMatrix[pos1][pos2]

    print(output1)
    print(output2)

def max_2dmatrix(sequence):
    if not sequence:
        raise ValueError('empty sequence')

    maximum = sequence[0]

    for item in sequence:
        # Compare elements by their weight stored
        # in their second element.
        if max(item) > max(maximum):
            maximum = item

    return maximum

#Smith–Waterman algorithm
def SWalign(seq1, seq2, match=2, mismatch=-1, gap=-1):
    length1 = len(seq1)
    length2 = len(seq2)
    sMatrix, tbMatrix = [], []
    # initialize scoring matrix
    # 1. create matrix with all 0's
    for x in range(length2 + 1):
        sMatrix.append([0] * (length1 + 1))
        tbMatrix.append(['zero'] * (length1 + 1))
    # 2. initialize matrix with initial scores
    tbMatrix[0][0] = 'done'

    # Now start scoring, and filling the traceback matrix

    for x in range(1, length2 + 1):
        for y in range(1, length1 + 1):
            check = mismatch
            if seq1[y - 1] == seq2[x - 1]:
                check = match

            score = max(sMatrix[x - 1][y - 1] + check, sMatrix[x - 1][y] + gap, sMatrix[x][y - 1] + gap, 0)
            # this block is for generating tbMatrix
            if score == sMatrix[x - 1][y - 1] + check:
                tbMatrix[x][y] = 'diag'
            elif score == sMatrix[x - 1][y] + gap:
                tbMatrix[x][y] = 'up'
            elif score == sMatrix[x][y - 1] + gap:
                tbMatrix[x][y] = 'left'
            else:
                tbMatrix[x][y] = 'zero'

            sMatrix[x][y] = score
    for x in range(length2+1):
        print(sMatrix[x])
    # now trace back the matrix to find the optimal alignment
    # this method will start from the highest score
    pos1 = sMatrix.index(max_2dmatrix(sMatrix))
    pos2 = sMatrix[pos1].index(max(map(max, sMatrix)))

    hold = tbMatrix[pos1][pos2]
    output1, output2 = '', ''
    while hold != 'done' and hold != 'zero':
        if hold == 'diag':
            output1 = seq1[pos2 - 1] + output1
            output2 = seq2[pos1 - 1] + output2
            pos1 -= 1
            pos2 -= 1
        elif hold == 'left':
            output1 = seq1[pos2 - 1] + output1
            output2 = '-' + output2
            pos2 -= 1
        elif hold == 'up':
            output1 = '-' + output1
            output2 = seq2[pos1 - 1] + output2
            pos1 -= 1
        hold = tbMatrix[pos1][pos2]

    print(output1)
    print(output2)