import math
from statistics import stdev

from Bio.Data import CodonTable
from collections import defaultdict

from uorf_predictor.instances import UORFCoordinates

def gc_content(five_sequence: str, uorf: UORFCoordinates) -> float:
    """
    Get the GC content of an uORF.

    If the number of bases chosen ends beyond the 5'UTR end limit (overlapping with the mORF), the number of bases taken will be clipped to
    the length between the uORF stop codon and the mORF start codon.

    :returns: the GC content as a float in range [0, 1].
    """
    total = uorf.uorf.end - uorf.uorf.start

    if uorf.uorf.end > len(five_sequence):
        raise ValueError("uORF overlaps with the mORF")

    if total == 0:
        return 0
    else:
        region = five_sequence[uorf.uorf.start: uorf.uorf.end]
        g = region.count("G")
        c = region.count("C")
                
        return (g+c) / total


def intercistronic_distance(five_sequence: str, uorf: UORFCoordinates) -> int:
    """
    Calculate the intercistronic distance, defined as the distance from the uORF stop codon to the mORF start codon.

    See here: Silva, J., Fernandes, R., Romão, L. (2019). Translational Regulation by Upstream Open Reading Frames and Human Diseases. 
    In: Romão, L. (eds) The mRNA Metabolism in Human Disease. Advances in Experimental Medicine and Biology, vol 1157. 
    Springer, Cham. https://doi.org/10.1007/978-3-030-19966-1_5
    """
    if uorf.uorf.end > len(five_sequence):
        raise ValueError("uORF overlaps with the mORF")

    return len(five_sequence) - uorf.uorf.end


def cap_five_to_uorf_distance(uorf: UORFCoordinates) -> int:
    """
    Calculate the number of bases between the 5' cap and the uORF start codon.
    """
    return uorf.uorf.start


def kozak_sequence_strength(five_sequence: str, uorf: UORFCoordinates) -> int:
    """
    Indicate the difference between a given Kozak sequence and the consensus sequence, based on two residues.

    The Kozak consensus sequence is defined as CCRCCAUGG, with a purine base in position -3 and a guanine base in position +4 
    as most important for initiation. Initiation sequence contexts are frequently classified as strong (both critical residues match the consensus sequence),
    as adequate/intermediate (either residue -3 or +4 matches) or as weak (neither residue matches).

    :returns: integer being: 2 (strong), 1 (adequate) or 0 (weak).

    See here: Wethmar K, Smink JJ, Leutz A. Upstream open reading frames: molecular switches in (patho)physiology. 
    Bioessays. 2010 Oct;32(10):885-93. doi: 10.1002/bies.201000037. Epub 2010 Aug 19. PMID: 20726009; PMCID: PMC3045505.
    """
    purines = ["A", "G"]
    minus_three_residue = five_sequence[uorf.uorf.start - 3]
    plus_four_residue = five_sequence[uorf.uorf.start + 3]

    if minus_three_residue in purines and plus_four_residue == "G":
        return 2
    elif minus_three_residue in purines and plus_four_residue != "G":
        return 1
    elif minus_three_residue not in purines and plus_four_residue == "G":
        return 1
    else:
        return 0
    
class Codon_features:

    def __init__(
        self, 
        uorf_sequence: str,
    ):
        self._uorf_sequence = uorf_sequence
        self._weights = self._obtain_relative_weights()

    def _obtain_relative_weights(self) -> list:
            """
            Retrieve the relative weights of the codons.

            Krishnamurthy Subramanian, Bryan Payne, Felix Feyertag, David Alvarez-Ponce, 
            The Codon Statistics Database: A Database of Codon Usage Bias, Molecular Biology and Evolution, Volume 39, Issue 8, August 2022, msac157, 
            https://doi.org/10.1093/molbev/msac157
            """
            codon_frequency = {
                "AAA": 25.1901, "AAG": 31.4481, "AAT": 17.0444, "AAC": 18.3501,
                "ACA": 15.3184, "ACC": 18.4242, "ACG": 5.9236, "ACT": 13.5182,
                "AGA": 12.3146, "AGC": 19.9075, "AGG": 12.1033, "AGT": 12.6984,
                "ATA": 7.5636, "ATC": 19.3821, "ATT": 15.7146, "ATG": 21.0300,

                "CAA": 12.8031, "CAG": 34.5800, "CAT": 11.1677, "CAC": 15.0600,
                "CCA": 17.5929, "CCC": 20.4863, "CCG": 7.4462, "CCT": 18.2084,
                "CGA": 6.0584, "CGC": 10.3363, "CGG": 11.3963, "CGT": 4.4912,
                "CTA": 7.1547, "CTC": 19.0061, "CTG": 38.7259, "CTT": 13.3863,

            
                "GAA": 30.4087, "GAG": 40.0424, "GAT": 22.0502, "GAC": 24.8107,
                "GCA": 16.1205, "GCC": 27.8499, "GCG": 7.6763, "GCT": 18.3994,
                "GGA": 16.6903, "GGC": 22.2317, "GGG": 16.4639, "GGT": 10.6731,
                "GTA": 7.1527, "GTC": 13.9037, "GTG": 27.0722, "GTT": 10.9735,

                "TAA": 0.4856, "TAG": 0.3822, "TAT": 11.8554, "TAC": 14.3091,
                "TCA": 12.9278, "TCC": 17.7956, "TCG": 4.5676, "TCT": 15.7052, 
                "TGA": 0.8420, "TGC": 12.3082, "TGG": 12.2384, "TGT": 10.7211, 
                "TTA": 7.8952, "TTC": 19.1588, "TTG": 12.9418, "TTT": 16.9475,
            }
            # Normalize the table to relative weigths
            standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
            codon_to_aa = standard_table.forward_table 

            aa_to_codons = defaultdict(list)
            for codon, freq in codon_frequency.items():
                if codon in codon_to_aa:
                    aa = codon_to_aa[codon]
                    aa_to_codons[aa].append((codon, freq))

            cai_weights = {}
            for aa, codon_list in aa_to_codons.items():
                max_freq = max(freq for _, freq in codon_list)
                for codon, freq in codon_list:
                    cai_weights[codon] = freq / max_freq

            codons = [self._uorf_sequence[i:i+3] for i in range(0, len(self._uorf_sequence), 3)]

            weights = []
            for codon in codons:
                weight = cai_weights.get(codon, 0.0)
                weights.append(weight)

            if not weights:
                return 0.0

            return weights
    
    def codon_adaptation_index(self) -> float:
        """
        Calculate the Codon Adaptation Index (CAI) as geometric mean of codon weights.
        """
        log_weights = [math.log(w) for w in self._weights if w > 0]
        if not log_weights:
            return 0.0
        return math.exp(sum(log_weights) / len(log_weights))
    
    def codon_frequency_std(self) -> float:
        """
        Calculate the standard deviation of the codon usage weights.
        """
        return stdev(self._weights)


def codon_count(uorf_sequence: str) -> int:
    """
    Count the number of codon in the uORF.
    Only for non-overlapping uORFs.
    """
    if len(uorf_sequence) % 3 == 0:
        return (len(uorf_sequence) / 3)
    else:
        return 0