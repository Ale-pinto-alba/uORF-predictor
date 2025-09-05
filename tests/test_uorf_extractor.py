import pytest

from uorf_predictor.instances import FiveUTRCoordinates
from uorf_predictor.uorf_extractor import fetch_cdna_from_ensembl, obtain_uorf_in_five_utr, check_start_and_stop_codon


@pytest.mark.online
@pytest.mark.parametrize(
     "tx_id, start, end, n_bases",
     [
          ("ENST00000696628", "GGTCGTTCCC", "CTATTTGAAA", 2_412), # tx on the + strand
          ("ENST00000381418", "AGTTGCGCTT", "ATAAGGGTAA", 5_474), # tx on the - strand
     ]
)
def test_fetch_cdna_from_ensembl(
     tx_id: str,
     start: str,
     end: str,
     n_bases: int,
):
     cdna = fetch_cdna_from_ensembl(tx_id)
     assert cdna.startswith(start)
     assert cdna.endswith(end)
     assert len(cdna) == n_bases


def test_uorf_forward_strand(
    ppp_five_utr_sequence: str,
    ppp_five_utr: FiveUTRCoordinates,
):
     uorf_in_five_utr = obtain_uorf_in_five_utr(
          five_utrs=ppp_five_utr, 
          start_uorf=48_872_553, 
          end_uorf=48_872_634,
     )
     
     uorf = "ATGAACGCGCTGGCCTCCCTAACCGTCCGGACCTGTGATCGCTTCTGGCAGACCGAACCGGCGCTCCTGCCCCCGGGGTGA"

     assert ppp_five_utr_sequence[uorf_in_five_utr.uorf.start: uorf_in_five_utr.uorf.end] == uorf
     assert uorf_in_five_utr.uorf.end - uorf_in_five_utr.uorf.start == len(uorf)


def test_uorf_reverse_strand(
    hr_five_utr_sequence: str,
    hr_five_utr: FiveUTRCoordinates,
):
     uorf_in_five_utr = obtain_uorf_in_five_utr(
          five_utrs=hr_five_utr, 
          start_uorf=22_130_604, 
          end_uorf=22_130_708,
     )
     
     uorf = "ATGGCGCAACCTACGGCCTCGGCCCAGAAGCTGGTGCGGCCGATCCGCGCCGTGTGCCGCATCCTGCAGATCCCGGAGTCCGACCCCTCCAACCTGCGGCCCTAG"

     assert hr_five_utr_sequence[uorf_in_five_utr.uorf.start: uorf_in_five_utr.uorf.end] == uorf
     assert uorf_in_five_utr.uorf.end - uorf_in_five_utr.uorf.start == len(uorf)


@pytest.mark.parametrize(
     "uorf_sequence, expected",
     [
          ("ATGCCCTAG", True), 
          ("ATGCCCTCC", False), 
          ("CTGCCCTAG", True), 
          ("ACGCCCTAA", True), 
     ]
)
def test_check_start_and_stop_codon(
     uorf_sequence: str,
     expected: bool,
):
     assert check_start_and_stop_codon(uorf_sequence=uorf_sequence) == expected