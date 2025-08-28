import numpy as np

import RNA

from uorf_predictor.instances import UORFCoordinates


class RNA_folding:
    """
    `RNA_folding` enables the extraction of features related to RNA folding.

    For this class, it has been used the Python's API of the package ViennaRNA.
    See here: Lorenz, Ronny and Bernhart, Stephan H. and Höner zu Siederdissen, Christian and Tafer, Hakim and Flamm, Christoph and Stadler, Peter F. and Hofacker, Ivo L.
    ViennaRNA Package 2.0
    Algorithms for Molecular Biology, 6:1 26, 2011, doi:10.1186/1748-7188-6-26
    """
    def __init__(
        self,
        five_utr_sequence: str,
        uorf: UORFCoordinates,
    ):
        self._sequence = five_utr_sequence
        self._uorf = uorf
        self._fc_sequence = self._fc_sequence_generator()
        self._sequence_structure = self._sequence_structure_generator()
        self._uorf_start_context = self._uorf_start_codon_context()
        self._uorf_sequence = self._uorf_sequence_generator()
        self._fc_uorf_start_context = self._fc_uorf_start()
        self._fc_uorf_sequence = self._fc_uorf_sequence_generator()
        self._uorf_start_structure = self._uorf_start_context_structure_generator()
        self._uorf_sequence_structure = self._uorf_sequence_structure_generator()
        self._uorf_ubox = self._generate_ubox_probs()

    def _fc_sequence_generator(self) -> RNA.fold_compound:
        """"
        Return the `fold_compound` of the 5'UTR sequence. This object is the core data structure of the ViennaRNA package, that holds data related
        to the sequence, stores folding models and energy parameters, matrices and constraints.
        """
        return RNA.fold_compound(self._sequence)
    
    def _uorf_start_codon_context(self) -> str:
        """
        Get the start codon context of the uORF (30 nucleotides up- and downstream).
        """
        return self._sequence[self._uorf._uorf.start - 30: self._uorf._uorf.start - 30]
    
    def _uorf_sequence_generator(self) -> str:
        """
        Get the uORF sequence.
        """
        return self._sequence[self._uorf._uorf.start: self._uorf._uorf.end]
    
    def _fc_uorf_start(self) -> RNA.fold_compound: 
        """
        Return the `fold_compound` of the start codon context.
        """
        return RNA.fold_compound(self._uorf_start_context)
    
    def _fc_uorf_sequence_generator(self) -> RNA.fold_compound: 
        """
        Return the `fold_compound` of the uORF sequence.
        """
        return RNA.fold_compound(self._uorf_sequence)
    
    def _sequence_structure_generator(self) -> str:
        """
        Generate the structure as a dot-bracket representation of the 5'UTR sequence.
        """
        structure, mfe = self._fc_sequence.mfe()
        return structure
    
    def _uorf_start_context_structure_generator(self) -> str:
        """
        Retrieve the structure as a dot-bracket representation of the start codon context.
        """
        structure, mfe = self._fc_uorf_start_context.mfe()
        return structure
    
    def _uorf_sequence_structure_generator(self) -> str:
        """
        Retrieve the structure as a dot-bracket representation of the uORF sequence.
        """
        structure, mfe = self._fc_uorf_sequence.mfe()
        return structure
    
    def five_sequence_mfe(self) -> float:
        """
        Stores the Minimum Free Energy (MFE) ensemble of the 5'UTR sequence.

        Returns: an `float` with the MFE value.
        """
        structure, mfe = self._fc_sequence.mfe()
        return mfe

    def uorf_start_context_mfe(self) -> float:
        """
        Stores the Minimum Free Energy (MFE) ensemble of the uORF start codon context.

        Returns: an `float` with the MFE value.
        """
        structure, mfe = self._fc_uorf_start_context.mfe()
        return mfe
    
    def uorf_sequence_mfe(self) -> float:
        """
        Stores the Minimum Free Energy (MFE) ensemble of the uORF sequence.

        Returns: an `float` with the MFE value.
        """
        structure, mfe = self._fc_uorf_sequence.mfe()
        return mfe

    def start_codon_is_unpaired(self) -> bool:
        """
        Determines if all three bases of the start codon are unpaired.

        Returns: a `bool` with True if all bases unpaired or False if not.
        """
        start_codon_dot = self._sequence_structure[self._uorf._uorf.start:self._uorf._uorf.start + 3]
        return all(dot == '.' for dot in start_codon_dot)
    
    def start_codon_number_unpaired(self) -> int:
        """
        Counts the number of unpaired bases in the start codon.

        Returns: an `int` with the number.
        """
        start_codon_dot = self._sequence_structure[self._uorf._uorf.start:self._uorf._uorf.start + 3]
        return start_codon_dot.count(".")
    
    def uorf_unpaired_bases_percentage(self) -> float:
        """
        Calculate the percentage of unpaired bases in the uORF.

        Returns: a `float` with the percentage.
        """ 
        number_unpaired = self._uorf_sequence_structure.count('.')
        return (number_unpaired * 100) / len(self._uorf_sequence_structure)    

    def uorf_ensemble_diversity(self) -> float:
        """
        Calculate the difference between the Ensemble diversity of the uORF.

        Returns: an `float` with the value.
        """ 
        self._fc_uorf_sequence.pf() 

        return self._fc_uorf_sequence.mean_bp_distance()

    def _generate_ubox_probs(
        self, 
        threshold: float = 1e-5,
    ) -> list:
        """
        Stores the ubox base pairing probabilities of the sequence.

        Args:
            threshold: a `float` to filter base pair probabilities.

        Returns: a `list` containing the ubox base pair probabilities.
        """
        ubox = []
        self._fc_uorf_sequence.pf()
        basepair_probs = self._fc_uorf_sequence.bpp()
        for i in range(1, len(self._uorf_sequence)+1):
            for j in range(i+1, len(self._uorf_sequence)+1):
                p = basepair_probs[i][j]
                if p > threshold:
                    ubox.append({"i": i, "j": j, "score": p})
        return ubox

    def ubox_total_probs_sum(
        self, 
    ) -> float:
        """
        Calculates the sum of all ubox probabilities of the uORF.

        Returns: a `float` with the sum.
        """
        return sum(entry["score"] for entry in self._uorf_ubox)

    def shannon_entropy(
        self,
    ) -> float:
        """
        Calculates the Shannon Entropy of the uORF.

        Returns: a `float` with the entropy value.
        """
        probs = [0.0] * (len(self._uorf_sequence) + 1)

        for entry in self._uorf_ubox:
                i = entry['i']
                j = entry['j']
                p = entry['score']
                probs[i] += p
                probs[j] += p

        total_prob = sum(probs)
        if total_prob == 0:
            return 0.0

        probs = [p / total_prob for p in probs if p > 0]
        entropy = -sum(p * np.log2(p) for p in probs)
        return entropy