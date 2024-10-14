from math import sqrt

from Bio.Seq import Seq
from numpy import var

def read_data(file_name):
    virus = open(file_name).read()
    name_end = virus.find("\n")
    name = virus[1:name_end]
    name = name.split(".")[0]
    sequence = virus[name_end+1:]
    return [name, sequence]


def normalize_data(seq):
    seq = seq[seq.find("\n")+1:]
    return seq.replace("\n", "")


def find_pairs(virus):
    start_codons = ["ATG"]
    end_codons = ["TGA", "TAG", "TAA"]

    pairs = []

    s = 0

    while s < len(virus) - 2:
        codon = virus[s: s + 3]
        if codon in start_codons:
            for e in range(s+3, len(virus) - 2, 3):
                codon2 = virus[e: e + 3]
                if codon2 in end_codons:
                    if e - s < 100:
                        break
                    pairs.append([s, e])
                    s = e + 2
                    break
        s = s + 1

    return pairs


def pair_to_sequence(virus, pair):
    return virus[pair[0]: pair[1] + 3]


def find_sequences(virus):

    complement_virus = str(Seq(virus).reverse_complement())

    pairs = find_pairs(virus) + find_pairs(complement_virus)

    sequences = []

    for pair in pairs:
        sequences.append(pair_to_sequence(virus, pair))

    translated_sequences = []

    for seq in sequences:
        translated_sequences.append(Seq(seq).translate())

    return translated_sequences


def find_codon_frequencies(virus):

    translated_sequences = find_sequences(virus)


    codon_frequencies = {}
    codon_count = 0

    for seq in translated_sequences:
        for codon in seq:
            if codon_frequencies.get(codon) is None:
                codon_frequencies[codon] = 1
            else:
                codon_frequencies[codon] += 1
            codon_count += 1

    for key in codon_frequencies.keys():
        codon_frequencies[key] = codon_frequencies[key] / codon_count

    return codon_frequencies


def find_dicodon_frequencies(virus):

    translated_sequences = find_sequences(virus)
    dicodon_frequencies = {}

    dicodon_count = 0

    for seq in translated_sequences:
        for i in range(len(seq) - 1):
            dicodone = str(seq[i: i + 2])
            if dicodon_frequencies.get(dicodone) is None:
                dicodon_frequencies[dicodone] = 1
            else:
                dicodon_frequencies[dicodone] += 1
            dicodon_count += 1

    for key in dicodon_frequencies.keys():
        dicodon_frequencies[key] = dicodon_frequencies[key] / dicodon_count

    return dicodon_frequencies


def find_distance(alphabet, frequencies1, frequencies2):

    distance = 0

    for value in alphabet:
        frequency1 = 0
        frequency2 = 0
        if frequencies1.get(value) is not None:
            frequency1 = frequencies1[value]
        if frequencies2.get(value) is not None:
            frequency2 = frequencies2[value]

        dist = (frequency1 - frequency2)
        distance += dist * dist

    return sqrt(distance)


def build_distance_matrix(frequencies, alphabet):
    distance_matrix = []
    for virus1 in frequencies:
        distances = []
        for virus2 in codon_frequencies:
            distances.append(find_distance(alphabet, virus1, virus2))
        distance_matrix.append(distances)

    return distance_matrix

def print_distance_matrix(matrix):
    print(len(matrix))

    for i, distances in enumerate(matrix):
        print(viruses_names[i], end='')
        print(" ", end='')
        for distance in distances:
            print(distance, end='')
            print(" ", end='')
        print()


bac1 = "viruses/data/bacterial1.fasta"
bac2 = "viruses/data/bacterial2.fasta"
bac3 = "viruses/data/bacterial3.fasta"
bac4 = "viruses/data/bacterial4.fasta"

mam1 = "viruses/data/mamalian1.fasta"
mam2 = "viruses/data/mamalian2.fasta"
mam3 = "viruses/data/mamalian3.fasta"
mam4 = "viruses/data/mamalian4.fasta"

viruses_file_names = [bac1, bac2, bac3, bac4, mam1, mam2, mam3, mam4]

viruses_names = []
viruses_data = []

for virus_file_name in viruses_file_names:
    data = read_data(virus_file_name)
    viruses_names.append(data[0])
    viruses_data.append(data[1])

normalized_data = []

for virus in viruses_data:
    normalized_data.append(normalize_data(virus))

codon_frequencies = []

for virus in normalized_data:
    codon_frequencies.append(find_codon_frequencies(virus))

dicodon_frequencies = []

for virus in normalized_data:
    dicodon_frequencies.append(find_dicodon_frequencies(virus))

codon_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F','P', 'S', 'T', 'W', 'Y', 'V']
dicodon_acids = [i+j for i in codon_acids for j in codon_acids]

codon_distance_matrix = build_distance_matrix(codon_frequencies, codon_acids)
dicodon_distance_matrix = build_distance_matrix(dicodon_frequencies, dicodon_acids)

print_distance_matrix(codon_distance_matrix)
print()
print_distance_matrix(dicodon_distance_matrix)
print()

def find_most_variant(frequencies, alphabet):
    codon_variations = []
    for acid in alphabet:
        codon_frequencies_acid = []
        for codon_frequency in frequencies:
            frequency = codon_frequency.get(acid)
            if frequency is not None:
                codon_frequencies_acid.append(frequency)
            else:
                codon_frequencies_acid.append(0)

        codon_variations.append([acid, var(codon_frequencies_acid)])

    return sorted(codon_variations, key=lambda x: x[1], reverse=True)



print(find_most_variant(codon_frequencies, codon_acids))
print(find_most_variant(dicodon_frequencies, dicodon_acids))



