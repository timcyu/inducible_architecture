def best_A_content(oligo):
    '''
    Choose the strand with the lowest A content because A's are harder to
    synthesize. Return true if sequence needs to be reverse complemented
    '''
    rc_oligo = reverse_complement(oligo)

    oligo_As = sum([1 for nt in oligo if nt == 'A'])
    rc_As = sum([1 for nt in rc_oligo if nt == 'A'])

    if oligo_As < rc_As:
        return False
    else:
        return True


def reverse_complement(sequence):
    """Return the reverse complement of the supplied sequence string """
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '.': '.', '*': '*'}
    # Reverse, convert to uppercase
    sequence = sequence[::-1].upper()
    # Complement
    letters = list(sequence)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)

# ======================================================================================================================
# up sequences
up_326x = "GGAAAATTTTTTTTCAAAAGTA" # 22 bp
up_136x = "GAAAATATATTTTTCAAAAGTA"
up_69x =  "AGAAAATTATTTTAAATTTCCT"
no_up =   "AGCTCATTCATTAGGCACCCCA" # taken directly from lacI

# -35 sequences
minus35cons = "TTGACA"
minus35_31A_32C = "TTGCAA"
minus35_33T = "TTTACA"
minus35_30C_33T = "TTTACC"
minus35_31G_33T_35C = "CTTAGA"

# -10 sequences
minus10cons = "TATAAT"
minus10_12G = "GATAAT"
minus10_12A = "AATAAT"
minus10_7A = "TATAAA"

ext_UV5 =   "TCG" # 3 bp
ext_min10 = "TGG" # G is conserved at position 13

# Primers for cloning purposes
# >skpp-218-F sk20mer-6672
fwd_primer = "AATTATTCCTGCTGAGGGCTCGAG"

# >skpp-232-R sk20mer-9146
rev_primer = "GCTAGCGTTTTATAGGGTGCTCTG"


five_prime_bg =  "CCCGCGCGTTGGCCGATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCG" # 60 bp
three_prime_bg = "ACACACAGGAAACAGCTATGACCATGATTACGGATTCACT" # 40 bp, changed 5' T to A
lacUV5_core = "AATGTAAGTTAGCTCATTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGG" # 71 bp

minus35s = [minus35cons, minus35_31A_32C, minus35_33T, minus35_30C_33T, minus35_31G_33T_35C]
minus10s = [minus10cons, minus10_12G, minus10_12A, minus10_7A]
up = [up_326x, up_136x, up_69x, no_up]
ext = [ext_UV5, ext_min10]

minus35s_names = ["minus35cons", "minus35_31A_32C", "minus35_33T", "minus35_30C_33T", "minus35_31G_33T_35C"]
minus10s_names = ["minus10cons", "minus10_12G", "minus10_12A", "minus10_7A"]
up_names = ["up_326x", "up_136x", "up_69x", "no_up"]
ext_names = ["ext_UV5", "ext_min10"]

# 4. [UP] Varying UP sequence position relative to -35
# 15 shift * 4 UP * 5 min35 * 2 extmin10 * 4 min10 = 2400
up_seq_lib = []
up_name = []
temp = ""
temp_name = ""
loc_up = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

for shift in loc_up:
    for upseq, upseq_name in zip(up, up_names):
        for min35, min35_name in zip(minus35s, minus35s_names):
            for extmin10, extmin10_name in zip(ext, ext_names):
                for min10, min10_name in zip(minus10s, minus10s_names):
                    temp = rev_primer[::-1] + three_prime_bg[0:40][::-1] + lacUV5_core[::-1][0:6] + min10[::-1] + extmin10[::-1] + lacUV5_core[::-1][15:29] + min35[::-1] + lacUV5_core[::-1][35:35+shift] + upseq[::-1] + lacUV5_core[::-1][57+shift:71] + five_prime_bg[::-1][0:39] + fwd_primer[::-1]
                    temp_name = "UP" + "-" + upseq_name + "-" + "shift" + str(shift) + "-" + min35_name + "-" + extmin10_name + "-" + min10_name
                    up_seq_lib.append(temp[::-1])
                    up_name.append(temp_name)

count = 0
for i in up_seq_lib:
    print len(i)
    count = count + 1

print count

# mark reversed flipped
for seq, index in zip(up_seq_lib, range(len(up_name))):
    reverse = best_A_content(seq)
    if reverse:
        up_name[index] += '-flipped'
        up_seq_lib[index] = reverse_complement(up_seq_lib[index])

# Creates a fasta format file containing all sequences and names
file = open("upseq.txt", "w")
for i in range(len(up_seq_lib)):
    file.write(">" + up_name[i] + "\t" + up_seq_lib[i] + "\n")
file.close()


string = up_seq_lib[0]
print string.count('A')
print string.count('T')