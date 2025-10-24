# Script projet Chapitre 1 from Bioinformatics Algorithms, Phillip Compeau & Pavel Pevzner
# Justin Acoulon, Kenzo Bouillé, Léna VIN 

#### Code Challenge 1A ####
# boucle avec range qui gÃ©nÃ¨re la suite des entiers
def PatternCount(text,pattern):
    count = 0
    for i in range(len(text)-len(pattern) + 1):
        if text[i:i+len(pattern)] == pattern: # Extraction du motif # count n'est pas utilisÃ© ici car cette mÃ©thode ne gÃ¨re pas les overlaps.
            count = count + 1
    return count

#Test
text = "ACAACTATGCATACTATCGGGAACTATCCT"
pattern = "ACTAT"
print(PatternCount(text,pattern)) # 3

#### Code Challenge 1B ####
""" Version 1 du code avec une liste"""
def frequentWords(text, k):
    FrequentPatterns = []
    COUNT = []
    
    # Comptage des occurrences pour chaque k-mer
    for i in range(len(text) - k + 1):
        Pattern = text[i:i+k]
        COUNT.append(pattern_count(text, Pattern))

    maxCount = max(COUNT)

    # Sélection des motifs les plus fréquents
    for i in range(len(text) - k + 1):
        if COUNT[i] == maxCount:
            FrequentPatterns.append(text[i:i+k])

    # Suppression des doublons
    FrequentPatterns = list(set(FrequentPatterns))
    return FrequentPatterns, maxCount

"""Version 2 du code avec un dictionnaire, plus optimiser"""
def FrequentWords(text,k):
    count = {} #dictionnaire vide qui va stocker les kmers avec leur nombre d'occurences

    #compter tous les k-mer dans texte
    for i in range(len(text)-k+1): #nombre de kmers possibles dans la sÃ©quence (text). Donne les positions possibles dans la sequence de chaque kmer en fonction de k
        pattern = text[i:i+k] #parcourir la sÃ©quence lettre aprÃ¨s lettre pour prendre chaque kmer de longueur k
        if pattern in count: #test d'appartenance
            count[pattern] = count[pattern] + 1 #si le kmer est dÃ©jÃ  dans le dictionnaire, on ajoute 1
        else:
            count[pattern] = 1 #si le kmer n'est pas dans le dictionnaire, on l'ajoute aveec la valeur 1
    
    #kmer le plus nombreux
    max_count = max(count.values()) #valeur maximale dans le dictionnaire
    frequent_patterns = [] #liste vide pour stocker les kmers les plus frÃ©quents prÃ©sents dans le dictionnaire (il peut y avoir plusieurs kmer avec le mÃªme nombre maximal d'occurences)
    ##### Ã©tape Ã  supprimer si on utilise set (pour Ã©viter les doublons (afficher 3 fois le kmer qui apparait 3 fois par exemple))

    for pattern in count: #parcours tous les kmers du dictionnaire
        if count[pattern] == max_count: #si ce kmer a la frequence max
            frequent_patterns.append((pattern, max_count)) #le kmer est ajoutÃ© Ã  la liste des rÃ©sultats avec son nombre d'occurence (tuple)
    return frequent_patterns

#### Code Challenge 1C : Reverse complement problem ####
"""Version 1 sans dictionnaire"""
def ReverseStrand(seq):
    # Remplacer chaque nuclÃ©otide par son complÃ©ment
    complement_nuc = ""
    for nucleotide in seq:
        if nucleotide == "A":
            complement_nuc = complement_nuc + "T"
        elif nucleotide == "T":
            complement_nuc = complement_nuc + "A"
        elif nucleotide == "C":
            complement_nuc = complement_nuc + "G"
        elif nucleotide == "G":
            complement_nuc = complement_nuc + "C"
    
    # sens 5'->3'
    reverse_seq_comp = ""
    for i in range(len(complement_nuc) - 1, -1, -1):
        reverse_seq_comp = reverse_seq_comp + complement_nuc[i]
    
    return reverse_seq_comp

"""Version 2 avec dictionnaire"""
def reverse_complementary(text):
    complement ={
        'A':'T',
        'a':'t',
        'T':'A',
        't':'a',
        'C':'G',
        'c':'g',
        'G':'C',
        'g':'c'
    }
    result = ""

    # Pour chaque lettre du motif (en partant de la fin)
    for base in text[::-1]:
        result += complement[base]  # on ajoute la base complémentaire

    return result

text= "AGTCGCATAGT"
print(reverse_complementary(text)) #ACTATGCGACT

#### Code Challenge 1D : Pattern matching problem ####

def PatternOccurence(seq,pattern):
    count = 0 #compteur
    """On avait pas penser à cette partie au départ, mais permet de faire moins d'erreurs"""
    seq = seq.upper()  # Normaliser en majuscules
    pattern = pattern.upper() # Optionnel si pattern déjà en majuscules
    for i in range(len(seq)- len(pattern)+1): #on parcours toutes les positions possibles dans la sÃ©quences
        test_seq = seq[i:i+len(pattern)] #on va tester la sequence 
        if test_seq == pattern:
            count = count + 1
    return count

#test
seq="atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaacctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgaccacggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgacttgtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggattacgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttaggatagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaattgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaagatcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtttccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc"
pattern="atgatcaag"
print(PatternOccurence(seq,pattern)) #3

#### Code Challenge 1E : Find patterns forming clumps in a str ####

def FindClumps(text, k, L, t):

    clumps = []  # Liste des k-mers qui forment des clumps

    # Parcourir chaque position de dÃ©part possible pour une fenÃªtre de taille L
    for i in range(len(text) - L + 1):
        
        # Extraire la fenÃªtre de taille L
        window = text[i:i+L]
        
        # Pour chaque k-mer possible dans cette fenÃªtre
        for j in range(len(window) - k + 1):
            pattern = window[j:j+k]
            
            # Compter combien de fois ce pattern apparaÃ®t dans la fenÃªtre
            count = 0 #compteur
            for position in range(len(window) - k + 1):
                if window[position:position+k] == pattern:
                    count += 1
            
            # Si le pattern apparaÃ®t assez souvent et n'est pas dÃ©jÃ  dans la liste
            if count >= t and pattern not in clumps:
                clumps.append(pattern)
    
    return clumps

# Test avec l'exemple du livre
text = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA"
k = 5
L = 50  
t = 4

print(FindClumps(text, k, L, t)) #['CGACA', 'GAAGA']


# Test avec seq oriC of V. cholerae
text="CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA"
k = 2
L = 50  
t = 4
print(FindClumps(text, k, L, t)) #['CG', 'GA', 'AC', 'CA', 'AG', 'TG', 'AA']


#### Code Challenge 1F : Minimum Skew Problem ####

def MinimumSkew(DNAstring):
    count_SKEW = []
    count = 0
    for i in range(len(DNAstring)):
        if DNAstring[i:i+1] == "C": 
            count = count + 1
            count_SKEW.append(count)
        else: 
            count = count + 1
    return (count_SKEW)

DNAstring = "CCCCGTTTGGGAAACCCCCCTTT"
print(MinimumSkew(DNAstring)) #[1, 2, 3, 4, 15, 16, 17, 18, 19, 20]

#### Code Challenge 1G : Hamming Distance Problem ####
#HammingDistance = count(p!=q)

def HammingDistance(String1,String2):
    assert len(String1)==len(String2)
    count=0
    for i in range(len(String1)):
         if String1[i] == String2[i]: 
            count = count 
         else:
             count += 1
    return(count)

String1 = "AGTAAAGAA"
String2 = "AGTGAAAAA"
print(HammingDistance(String1,String2)) #2 diffÃ©rences

#### Code Challenge 1H : Approximate Pattern Matching Problem ####
#mixer hammingdistance et patterncount
# prendre tous les kmer 1 Ã  1 le long de la sÃ©quence et les comparer au pattern. Si diff < d alors marquer la position

# crÃ©er une liste avec tous les kmer de la mÃªme taille que pattern. Comparer chaque Ã©lÃ©ment de cette liste Ã  pattern et isoler ceux avec d<x
# Code initial: pb return mal placÃ© (liste retournÃ©e trop tÃ´t, dans la premiÃ¨re boucle) + pb de logique: recherche des k-mers uniques alors qu'on cherche les positions + manque HammingDistance + syntaxe append incorrecte
def ApproximatePatternMatching(text,pattern,d):
    list_kmer = []
    for i in range(len(text)-len(pattern)+1):
        kmer=text[i:i+len(pattern)]
        if kmer in list_kmer:
            pass
        else:
            list_kmer.append(kmer)
        return(list_kmer)
    d_treshold = []
    count = 0
    for j in list_kmer:
        if list_kmer[j] == pattern:
            count = count
        else:
            count += 1 and d_treshold.append(list_kmer[j],count)
        return(d_treshold)
    final_list=[]
    for k in d_treshold:
        if k<=d:
            final_list.append(d_treshold[k])
        else:
            pass
    return(final_list.append)
text="actgactactacgacg"
pattern="act"
d=1
print(ApproximatePatternMatching(text,pattern,d)) #act

def ApproximatePatternMatching(text,pattern,d):
    final_list= [] #rÃ©sultat final des positions
    for i in range(len(text)-len(pattern)+1): #on parcourt tous les kmers possibles dans la sÃ©quence
        kmer=text[i:i+len(pattern)] #extraction des kmers Ã  la position i
#distance de Hamming (kmer Vs pattern)
        count = 0
        for j in range(len(kmer)):
            if kmer[j] != pattern[j]:
                count += 1
#si distance <=d on ajoute la position i dans le compteur
        if count <= d:
            final_list.append(i+1) #python compte Ã  partir de 0

    return(final_list)
text="actgactactacgacg"
pattern="act"
d=1
print(ApproximatePatternMatching(text,pattern,d)) # [14,14,14] = pb d'indented
#[1, 5, 8, 11, 14] ok

#### Code Challenge 1J : Frequent Words with Mismatches and Reverse Complements Problem ####
# on prend en compte les mismatches et les reverse complement.
#1 parcourir tous les kmer de la sequence
#2 donner le RC
#3 pour chaque pattern, comparaison avec tous les kmers de la seq
#4 distance de Hamming
#5 si distance <=d (+ ou - strand) on incremente le ciompteur
#6 extraction pattern avec le max

def FrequentWordsWithMismatchesAndReverse(text,k,d):
    count = {}
    def ReverseStrand(seq):
        complement_nuc=""
        for nucleotide in seq:
            if nucleotide == "A":
                complement_nuc = complement_nuc + "T"
            elif nucleotide == "T":
                complement_nuc = complement_nuc + "A"
            elif nucleotide == "C":
                complement_nuc = complement_nuc + "G"
            elif nucleotide == "G":
                complement_nuc = complement_nuc + "C"
        reverse_seq_comp = ""
        for i in range(len(complement_nuc) - 1, -1, -1):
            reverse_seq_comp = reverse_seq_comp + complement_nuc[i]
        return reverse_seq_comp
#tous les kmers possibles
    for i in range(len(text) - k + 1):
        pattern = text[i:i+k]
        revcomp = ReverseStrand(pattern)  # complÃ©ment inverse du kmer
        count[pattern] = 0
        
        for j in range(len(text) - k + 1): #comaparaison
            test_seq = text[j:j+k]
            diff_count_1 = 0 #distance de H
            for x in range(k):
                if pattern[x] != test_seq[x]:
                    diff_count_1 += 1
            diff_count_2 = 0
            for x in range(k):
                if revcomp[x] != test_seq[x]:
                    diff_count_2 += 1
            if diff_count_1<=d or diff_count_2<=d:
                count[pattern] +=1
    # Identifier la frÃ©quence maximale
    max_count = max(count.values())

    # Lister les patterns (et leurs reverse) les plus frÃ©quents
    frequent_patterns = []
    for pattern in count:
        if count[pattern] == max_count:
            frequent_patterns.append((pattern, max_count))

    return frequent_patterns


# Test
text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
k = 4
d = 1

print(FrequentWordsWithMismatchesAndReverse(text, k, d)) #[('GCAT', 8), ('ATGT', 8), ('ATGA', 8), ('ATGC', 8)]     

#### Code Challenge 1K : Computing the Frequency Array of a Str ####
#crÃ©ation tableau de frÃ©quences = combien de fois chaque kmer possible apparaÃ®t dans une sequence
#on va numÃ©roter les nuclÃ©otides A=0, C=1, G=2 et T=3

def FrequencyArray(text, k):
    #liste liste de taille 4^k (tous les k-mers possibles)
    freq_array = [0] * (4**k)

    # Fonction interne pour convertir un k-mer en entier (index du tableau)
    def PatternToNumber(pattern):
        number = 0
        for i in range(len(pattern)):
            # Calculer la "valeur" du nuclÃ©otide
            if pattern[i] == "A":
                value = 0
            elif pattern[i] == "C":
                value = 1
            elif pattern[i] == "G":
                value = 2
            elif pattern[i] == "T":
                value = 3
            # Mise Ã  jour : passage en base 4
            number = 4 * number + value
        return number

    # Parcourir tous les k-mers possibles du texte
    for i in range(len(text) - k + 1):
        pattern = text[i:i+k]
        j = PatternToNumber(pattern)     # convertir k-mer en index
        freq_array[j] = freq_array[j] + 1   # incrÃ©menter la frÃ©quence Ã  cet index

    return freq_array

# Test FrequencyArray
text = "ACGCGGCTCTGAAA"
k = 2
print(FrequencyArray(text, k)) #[2, 1, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 1, 1, 0]

























