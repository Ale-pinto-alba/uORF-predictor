import pytest

from gpsea.model.genome import Region
from uorf_predictor.instances import UORFCoordinates, FiveUTRCoordinates
from uorf_predictor.rna_features import RNA_folding

class TestRNA_folding:

    @pytest.mark.parametrize(
    "region, ouorf, expected",
    [
        (Region(start=16, end=67), False, -298.7),
        (Region(start=302, end=407), False, -298.7),
        (Region(start=510, end=576), False, -298.7),
        ((Region(start=606, end=623)), True, -298.7),
    ]
    )
    def test_mfe_five_prime(
        self,
        hr_five_utr_sequence: str,
        hr_five_utr: FiveUTRCoordinates,
        region: Region, 
        ouorf: bool,
        expected: float,
    ):
        uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf=ouorf)
        rna_instance = RNA_folding(five_utr_sequence=hr_five_utr_sequence, uorf=uorf)
        result = rna_instance.five_sequence_mfe()
        
        assert result == pytest.approx(expected, rel=1e-6)

    @pytest.mark.parametrize(
    "region, ouorf, expected",
    [
        (Region(start=16, end=67), False, -12.6),
        (Region(start=302, end=407), False, -15.1),
        (Region(start=510, end=576), False, -28.7),
        ((Region(start=606, end=623)), True, -18.6),
    ]
    )
    def test_mfe_start_context(
        self,
        hr_five_utr_sequence: str,
        hr_five_utr: FiveUTRCoordinates,
        region: Region, 
        ouorf: bool,
        expected: float,
    ):
        uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf=ouorf)
        rna_instance = RNA_folding(five_utr_sequence=hr_five_utr_sequence, uorf=uorf)
        result = rna_instance.uorf_start_context_mfe()
        
        assert result == pytest.approx(expected, rel=0.1)

    @pytest.mark.parametrize(
    "region, ouorf, expected",
    [
        (Region(start=16, end=67), False, -12.6),
        (Region(start=302, end=407), False, -39),
        (Region(start=510, end=576), False, -28.7),
    ]
    )
    def test_mfe_uorf(
        self,
        hr_five_utr_sequence: str,
        hr_five_utr: FiveUTRCoordinates,
        region: Region, 
        ouorf: bool,
        expected: float,
    ):
        uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf=ouorf)
        rna_instance = RNA_folding(five_utr_sequence=hr_five_utr_sequence, uorf=uorf)
        result = rna_instance.uorf_sequence_mfe()
        
        assert result == pytest.approx(expected, rel=0.1)

    @pytest.mark.parametrize(
    "region, ouorf, expected",
    [
        (Region(start=16, end=67), False, False),
        (Region(start=302, end=407), False, False),
        (Region(start=510, end=576), False, False),
        ((Region(start=606, end=623)), True, False),
    ]
    )
    def test_unpaired_start_codon(
        self,
        hr_five_utr_sequence: str,
        hr_five_utr: FiveUTRCoordinates,
        region: Region, 
        ouorf: bool,
        expected: bool,
    ):
        uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf=ouorf)
        rna_instance = RNA_folding(five_utr_sequence=hr_five_utr_sequence, uorf=uorf)
        result = rna_instance.start_codon_is_unpaired()
        
        assert result == expected

    @pytest.mark.parametrize(
    "region, ouorf, expected",
    [
        (Region(start=16, end=67), False, 1),
        (Region(start=302, end=407), False, 0),
        (Region(start=510, end=576), False, 2),
        ((Region(start=606, end=623)), True, 1),
    ]
    )
    def test_number_unpaired_start_codon(
        self,
        hr_five_utr_sequence: str,
        hr_five_utr: FiveUTRCoordinates,
        region: Region, 
        ouorf: bool,
        expected: int,
    ):
        uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf=ouorf)
        rna_instance = RNA_folding(five_utr_sequence=hr_five_utr_sequence, uorf=uorf)
        result = rna_instance.start_codon_number_unpaired()
        
        assert result == expected

    @pytest.mark.parametrize(
    "region, ouorf, expected",
    [
        (Region(start=16, end=67), False, 56.8),
        (Region(start=302, end=407), False, 40.9),
        (Region(start=510, end=576), False, 45.4),
    ]
    )
    def test_percentage_unpaired(
        self,
        hr_five_utr_sequence: str,
        hr_five_utr: FiveUTRCoordinates,
        region: Region, 
        ouorf: bool,
        expected: float,
    ):
        uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf=ouorf)
        rna_instance = RNA_folding(five_utr_sequence=hr_five_utr_sequence, uorf=uorf)
        result = rna_instance.uorf_unpaired_bases_percentage()
        
        assert result == pytest.approx(expected, rel=0.1)

    @pytest.mark.parametrize(
    "region, ouorf, expected",
    [
        (Region(start=16, end=67), False, 10.4),
        (Region(start=302, end=407), False, 29.5),
        (Region(start=510, end=576), False, 2.5),
    ]
    )
    def test_ensembl_diversity(
        self,
        hr_five_utr_sequence: str,
        hr_five_utr: FiveUTRCoordinates,
        region: Region, 
        ouorf: bool,
        expected: float,
    ):
        uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf=ouorf)
        rna_instance = RNA_folding(five_utr_sequence=hr_five_utr_sequence, uorf=uorf)
        result = rna_instance.uorf_ensemble_diversity()
        
        assert result == pytest.approx(expected, rel=0.1)

    @pytest.mark.parametrize(
    "region, ouorf, expected",
    [
        (Region(start=16, end=67), False, 10.4),
        (Region(start=302, end=407), False, 29.5),
        (Region(start=510, end=576), False, 18.3),
    ]
    )
    def test_ubox_sum(
        self,
        hr_five_utr_sequence: str,
        hr_five_utr: FiveUTRCoordinates,
        region: Region, 
        ouorf: bool,
        expected: float,
    ):
        uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf=ouorf)
        rna_instance = RNA_folding(five_utr_sequence=hr_five_utr_sequence, uorf=uorf)
        result = rna_instance.ubox_total_probs_sum()
        
        assert result == pytest.approx(expected, rel=0.1)

    @pytest.mark.parametrize(
    "region, ouorf, expected",
    [
        (Region(start=16, end=67), False, 5.0),
        (Region(start=302, end=407), False, 6.2),
        (Region(start=510, end=576), False, 5.3),
    ]
    )
    def test_shannon_entropy(
        self,
        hr_five_utr_sequence: str,
        hr_five_utr: FiveUTRCoordinates,
        region: Region, 
        ouorf: bool,
        expected: float,
    ):
        uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf=ouorf)
        rna_instance = RNA_folding(five_utr_sequence=hr_five_utr_sequence, uorf=uorf)
        result = rna_instance.shannon_entropy()
        
        assert result == pytest.approx(expected, rel=0.1)