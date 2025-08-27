import pytest
import typing

from gpsea.model.genome import Region
from uorf_predictor.instances import FiveUTRCoordinates, UORFCoordinates
from uorf_predictor.uorf_features import gc_content, intercistronic_distance, cap_five_to_uorf_distance, kozak_sequence_strength, codon_adaptation_index


@pytest.mark.parametrize(
    "region, expected",
    [
        (Region(start=16, end=67), ((13+19)/51)),
        (Region(start=302, end=407), ((31+43)/105)),
        (Region(start=510, end=576), ((26+28)/66)),
    ]
)
def test_gc_content(
    hr_five_utr_sequence: str, 
    hr_five_utr: FiveUTRCoordinates, 
    region: Region, 
    expected: float,
):
    uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf= False)
    result = gc_content(five_sequence=hr_five_utr_sequence, uorf=uorf)
    
    assert result == pytest.approx(expected, rel=1e-6)


def test_uorf_ends_out_of_five_prime(
    hr_five_utr_sequence: str,
    hr_five_utr: FiveUTRCoordinates,
):
    overlapping_uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=Region(start=700, end=800), ouorf=True)

    with pytest.raises(ValueError) as e:
        gc_content(five_sequence=hr_five_utr_sequence, uorf=overlapping_uorf)

    assert e.value.args == ("uORF overlaps with the mORF",)


@pytest.mark.parametrize(
        "region, expected",
        [
            ((Region(start=16, end=67)), 556),
            ((Region(start=302, end=407)), 216),
            ((Region(start=510, end=576)), 47),
        ]
)
def test_intercistronic_distances(
    hr_five_utr_sequence: str,
    hr_five_utr: FiveUTRCoordinates,
    region: Region,
    expected: int,
):
    uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf=False)

    assert intercistronic_distance(five_sequence=hr_five_utr_sequence, uorf=uorf) == expected


@pytest.mark.parametrize(
        "region, expected",
        [
            ((Region(start=16, end=67)), 16),
            ((Region(start=302, end=407)), 302),
            ((Region(start=510, end=576)), 510),
        ]
)
def test_five_cap_to_uorf_distance(
    hr_five_utr: FiveUTRCoordinates,
    region: Region,
    expected: int,
):
    uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf=False)

    assert cap_five_to_uorf_distance(uorf=uorf) == expected


@pytest.mark.parametrize(
        "region, ouorf, expected",
        [
            ((Region(start=16, end=67)), False, 2),
            ((Region(start=302, end=407)), False, 1),
            ((Region(start=510, end=576)), False, 1),
            ((Region(start=606, end=623)), True, 0),
        ]
)
def test_kozak_sequence_strength(
    hr_five_utr_sequence: str,
    hr_five_utr: FiveUTRCoordinates,
    region: Region,
    ouorf: bool,
    expected: typing.Optional[int],
):
    uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf=ouorf)

    assert kozak_sequence_strength(five_sequence=hr_five_utr_sequence, uorf=uorf) == expected


@pytest.mark.parametrize(
        "region, expected",
        [
            ((Region(start=16, end=67)), pytest.approx(expected=1.9, rel=0.1)),
            ((Region(start=302, end=407)), pytest.approx(expected=2.2, rel=0.1)),
            ((Region(start=510, end=576)), pytest.approx(expected=2.1, rel=0.1)),
            ((Region(start=606, end=623)), pytest.approx(expected=2.0, rel=0.1)),
        ]
)
def test_codon_adaptation_index(
    hr_five_utr_sequence: str,
    region: Region,
    expected: typing.Optional[int],
):
    uorf = hr_five_utr_sequence[region.start:region.end]

    assert codon_adaptation_index(uorf_sequence= uorf) == expected