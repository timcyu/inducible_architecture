### author: Timothy Yu

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
no_up =   "AGCTCATTCATTAGGCACCCCA" # taken directly from lacUV5

# AraC sites placed between the -10 and -35
AraCO2_core = "TGTGGACTTTTCTGCCG"
AraCI1_core = "TATGGATAAAAATGCTA"
AraCI2_core = "TCAGGTAGGATCCGCTA"
AraCFP_core = "TATGGATTAATCTGCTG"
AraCscram_core = "CTTTCGATGTGGCTCTA" # scrambled AraCO2

# AraC sites placed downstream or upstream
AraCO2 = "AATGTGGACTTTTCTGCCGTG"
AraCI1 = "CTTATGGATAAAAATGCTATG"
AraCI2 = "CGTCAGGTAGGATCCGCTAAT"
AraCFP = "CTTATGGATTAATCTGCTGTG"
AraCscram = "GACTTTCGATGTGGCTCTATG" # scrambled AraCO2

# -35 sequences
minus35cons = "TTGACA"
minus35_31A_32C = "TTGCAA"
minus35_33T = "TTTACA"
minus35_30C_33T = "TTTACC"

# -10 sequences
minus10cons = "TATAAT"
minus10_12G = "GATAAT"
minus10_12A = "AATAAT"
minus10_7A = "TATAAA"

ext_UV5 =   "TCG" # 3 bp
ext_min10 = "TGG" # G is conserved at position 13

five_prime_bg =  "CCCGCGCGTTGGCCGATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCG" # 60 bp
three_prime_bg = "ACACACAGGAAACAGCTATGACCATGATTACGGATTCACT" # 40 bp, changed 5' T to A

between_up_minus35 = "GC" # 3 bp (3' most T is omitted to preserve the natural location of UP)

# Primers for cloning purposes

# >skpp-168-F sk20mer-5241
fwd_primer = "ACGCCGTGATTGTTATACCTCGAG"
# >skpp-296-R sk20mer-47488
rev_primer = "GCTAGCACCTCTTTCCTCTCTAGG"

# Core sequence
lacUV5_core = "AATGTAAGTTAGCTCATTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGG" # 71 bp

# Store in lists for convenience
minus35s = [minus35cons, minus35_31A_32C, minus35_33T, minus35_30C_33T]
minus10s = [minus10cons, minus10_12G, minus10_12A, minus10_7A]
araCs = [AraCO2, AraCI1, AraCI2, AraCFP, AraCscram]
araCcores = [AraCO2_core, AraCI1_core, AraCI2_core, AraCFP_core, AraCscram_core]
up = [up_326x, up_136x, up_69x, no_up]
ext = [ext_UV5, ext_min10]

minus35s_names = ["minus35cons", "minus35_31A_32C", "minus35_33T", "minus35_30C_33T"]
minus10s_names = ["minus10cons", "minus10_12G", "minus10_12A", "minus10_7A"]
araCs_names = ["AraCO2", "AraCI1", "AraCI2", "AraCFP", "AraCscram"]
araCcores_names = ["AraCO2_core", "AraCI1_core", "AraCI2_core", "AraCFP_core", "AraCscram_core"]
up_names = ["up_326x", "up_136x", "up_69x", "no_up"]
ext_names = ["ext_UV5", "ext_min10"]

# store sequences in here as they are generated
raw_seq_lib =[]
names = []

# ======================================================================================================================
# Case 1: [COMBO] All combinations of araC sites in the original lac operon
temp = ""
temp_name = ""

for up_araC, up_araC_name in zip(araCs, araCs_names):
    for min35, min35_name in zip(minus35s, minus35s_names):
        for min10, min10_name in zip(minus10s, minus10s_names):
            for down_araC, down_araC_name in zip(araCs, araCs_names):
                temp = rev_primer[::-1] + three_prime_bg[0:14][::-1] + down_araC[::-1] + lacUV5_core[::-1][0:6] + min10[::-1] + lacUV5_core[::-1][12:29] + min35[::-1] + lacUV5_core[::-1][35:71] + five_prime_bg[::-1][0:8] + up_araC[::-1] + five_prime_bg[::-1][8:23] + fwd_primer[::-1]
                temp_name = "ARAC-COMBO" + "-" + up_araC_name + "-" + min35_name + "-" + min10_name + "-" + down_araC_name
                raw_seq_lib.append(temp[::-1])
                names.append(temp_name)

# ======================================================================================================================
# 2. [DISTAL] Inducible promoters with distal and core araC sites (with x bp spacer) make into 40-50 bp loop
# 11 shift * 5 araC * 4 -35 * 5 araC * 4 -10 = 4400

loc_distal = [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
loc_five_prime = [45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35]
for shift, five_prime_shift in zip(loc_distal, loc_five_prime):
    for up_araC, up_araC_name in zip(araCs, araCs_names):
        for min35, min35_name in zip(minus35s, minus35s_names):
            for araCcore, araCcore_name in zip(araCcores, araCcores_names):
                for min10, min10_name in zip(minus10s, minus10s_names):
                    temp = rev_primer[::-1] + three_prime_bg[0:29][::-1] + lacUV5_core[::-1][0:6] + min10[::-1] + araCcore[::-1] + min35[::-1] + lacUV5_core[::-1][35:35+shift] + up_araC[::-1] + five_prime_bg[::-1][0:five_prime_shift] + fwd_primer[::-1]
                    temp_name = "ARAC-DISTAL" + "-" + up_araC_name + "-" + "shift" + str(shift) + "-" + min35_name + "-" + araCcore_name + "-" + min10_name
                    raw_seq_lib.append(temp[::-1])
                    names.append(temp_name)

# ======================================================================================================================
# 3. [STERIC] Using steric hindrance and replacing -35 to create repression
# 4 UP * 5 araC * 2 ext-10 * 4 min10 * 5 araC * 2 positions = 1600

# 55 bp loop
for upseq, upseq_name in zip(up, up_names):
    for up_araC, up_araC_name in zip(araCs, araCs_names):
        for extmin10, extmin10_name in zip(ext, ext_names):
            for min10, min10_name in zip(minus10s, minus10s_names):
                for down_araC, down_araC_name in zip(araCs, araCs_names):
                    temp = rev_primer[::-1] + three_prime_bg[19:23][::-1] + down_araC[::-1] + three_prime_bg[0:19][::-1] + lacUV5_core[::-1][0:6] + min10[::-1] + extmin10[::-1] + up_araC[::-1] + between_up_minus35[::-1] + upseq[::-1] + lacUV5_core[::-1][60:71] + five_prime_bg[::-1][0:35] + fwd_primer[::-1]
                    temp_name = "ARAC-STERIC-55bpLOOP" + "-" + upseq_name + "-" + up_araC_name + "-" + extmin10_name + "-" + min10_name + "-" + down_araC_name
                    raw_seq_lib.append(temp[::-1])
                    names.append(temp_name)

# Short loop with araC site in normal O1 position
for upseq, upseq_name in zip(up, up_names):
    for up_araC, up_araC_name in zip(araCs, araCs_names):
        for extmin10, extmin10_name in zip(ext, ext_names):
            for min10, min10_name in zip(minus10s, minus10s_names):
                for down_araC, down_araC_name in zip(araCs, araCs_names):
                    temp = rev_primer[::-1] + three_prime_bg[0:23][::-1] + down_araC[::-1] + lacUV5_core[::-1][0:6] + min10[::-1] + extmin10[::-1] + up_araC[::-1] + between_up_minus35[::-1] + upseq[::-1] + lacUV5_core[::-1][61:71] + five_prime_bg[::-1][0:36]+ fwd_primer[::-1]
                    temp_name = "ARAC-STERIC-O1LOOP" + "-" + upseq_name + "-" + up_araC_name + "-" + extmin10_name + "-" + min10_name + "-" + down_araC_name
                    raw_seq_lib.append(temp[::-1])
                    names.append(temp_name)

# ======================================================================================================================
# 5. [MULTIPLE] Do multiple araC's help recruit araC to the region to increase repression?
# w

for first_up_araC, first_up_araC_name in zip(araCs, araCs_names):
    for second_up_araC, second_up_araC_name in zip(araCs, araCs_names):
        for min35, min35_name in zip(minus35s, minus35s_names):
            for min10, min10_name in zip(minus10s, minus10s_names):
                for down_araC, down_araC_name in zip(araCs, araCs_names):
                    temp = rev_primer[::-1] + three_prime_bg[0:4][::-1] + down_araC[::-1] + lacUV5_core[::-1][0:6] + min10[::-1] + lacUV5_core[::-1][12:29] + min35[::-1] + lacUV5_core[::-1][35:71] + five_prime_bg[::-1][0:8] + first_up_araC[::-1] + second_up_araC[::-1]  + five_prime_bg[::-1][8:12] + fwd_primer[::-1]
                    temp_name = "ARAC-MULTIPLE" + "-" + second_up_araC_name + "-" + first_up_araC_name + "-" + min35_name + "-" + min10_name + "-" + down_araC_name
                    raw_seq_lib.append(temp[::-1])
                    names.append(temp_name)



# ======================================================================================================================
# 6. [HYBRID] Hybrid AraC + LacI dual repressor

lacO1 = "AATTGTGAGCGGATAACAATT" # 21 bp
lacO2 = "AAATTGTAGCGAGTAACAACC"
lacO3 = "GGCAGTGAGCGCAACGCAATT"
lacOsym = "AAATTGTGAGCGCTCACAATT"
lacOscram = "TTAACGGTGTGCATAATAGAA" # scramble O1

lacIs = [lacO1, lacO2, lacO3, lacOsym, lacOscram]
lacIs_names = ["lacO1", "lacO2", "lacO3", "lacOsym", "lacOscram"]

hybrid_35s = [minus35cons, minus35_31A_32C, minus35_33T]
hybrid_10s = [minus10cons, minus10_12G, minus10_12A]

hybrid_35s_names = ["minus35cons", "minus35_31A_32C", "minus35_33T"]
hybrid_10s_names = ["minus10cons", "minus10_12G", "minus10_12A"]

for up_lacI, up_lacI_name in zip(lacIs, lacIs_names):
    for min35, min35_name in zip(hybrid_35s, hybrid_35s_names):
        for core_araC, core_araC_name in zip(araCcores, araCcores_names):
            for min10, min10_name in zip(hybrid_10s, hybrid_10s_names):
                for down_araC, down_araC_name in zip(araCs, araCs_names):
                    for down_lacI, down_lacI_name in zip(lacIs, lacIs_names):
                        temp = rev_primer[::-1] + three_prime_bg[0:5][::-1] + down_araC[::-1] + down_lacI[::-1] + lacUV5_core[::-1][0:6] + min10[::-1] + core_araC[::-1] + min35[::-1] + lacUV5_core[::-1][35:71] + five_prime_bg[::-1][0:5] + up_lacI[::-1] + five_prime_bg[::-1][5:11] + fwd_primer[::-1]
                        temp_name = "ARAC-HYBRID" + "-" + up_lacI_name + "-" + min35_name + "-" + core_araC_name + "-" + min10_name + "-" + down_araC_name + "-" + down_lacI_name
                        raw_seq_lib.append(temp[::-1])
                        names.append(temp_name)

# ======================================================================================================================

# TESTING PURPOSES

count = 0
for i in raw_seq_lib:
    print len(i)
    count = count + 1
print count

#count = 0
#for i in names:
#    print i
#    count = count + 1
#print count

# mark reversed flipped
for seq, index in zip(raw_seq_lib, range(len(names))):
    reverse = best_A_content(seq)
    if reverse:
        names[index] += '-flipped'
        raw_seq_lib[index] = reverse_complement(raw_seq_lib[index])

# Creates a fasta format file containing all sequences and names
file = open("araC_lib.txt", "w")
for i in range(len(raw_seq_lib)):
    file.write(">" + names[i] + "\t" + raw_seq_lib[i] + "\n")
file.close()

string = raw_seq_lib[0]
print string.count('A')
print string.count('T')
