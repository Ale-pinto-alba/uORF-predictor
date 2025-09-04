import os
import pytest

from gpsea.model.genome import GenomeBuild, GRCh38, GenomicRegion, Strand
from uorf_predictor.instances import FiveUTRCoordinates

def pytest_addoption(parser):
    parser.addoption(
        "--runonline", action="store_true", default=False, help="run online tests"
    )


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "online: mark test that require internet access to run"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runonline"):
        # --runonline given in cli: do not skip online tests
        return
    skip_online = pytest.mark.skip(reason="need --runonline option to run")
    for item in items:
        if "online" in item.keywords:
            item.add_marker(skip_online)


@pytest.fixture(scope="session")
def fpath_test_dir() -> str:
    return os.path.dirname(__file__)


@pytest.fixture(scope="session")
def fpath_data_dir(fpath_test_dir: str) -> str:
    return os.path.join(fpath_test_dir, "data")


@pytest.fixture(scope="session")
def genome_build() -> GenomeBuild:
    return GRCh38


@pytest.fixture(scope="session")
def hr_five_utr(
    genome_build: GenomeBuild
) -> FiveUTRCoordinates:
    """
    5'UTR Genomic regions corresponding to one of the transcripts of the HR gene (ENSEMBL transcript ID: `ENST00000381418.9`).

    Both Genomic Regions were obtained from the chromosome 8 GTF file.

    see here: https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000168453;r=8:22114419-22133384;t=ENST00000381418
    """
    contig = genome_build.contig_by_name("8")
    assert contig is not None

    return FiveUTRCoordinates(
        regions=(
            GenomicRegion(
                contig=contig,
                start=22_130_427,
                end=22_131_010,
                strand=Strand.POSITIVE
                ).with_strand(other=Strand.NEGATIVE),
            GenomicRegion(
                contig=contig,
                start=22_129_170,
                end=22_129_210,
                strand=Strand.POSITIVE
                ).with_strand(other=Strand.NEGATIVE),
        )
    )


@pytest.fixture(scope="session")
def hr_five_utr_sequence() -> str:
    """
    5'UTR cDNA sequence of the transcript of the HR gene (ENSEMBL transcript ID: `ENST00000381418.9`) taken directly from
    the ENSEMBL website.
    
    see here: https://www.ensembl.org/Homo_sapiens/Transcript/Sequence_cDNA?db=core;g=ENSG00000168453;r=8:22114419-22133384;t=ENST00000381418.
    """
    return "AGTTGCGCTTCTGGCGATGGCGATCAGAGGTCCTGCTGCGCTCTCCGCCG" \
        + "CGCTCTACCTCCATTAGCCGCGCTGCGCGGTGCTGCGCCCTCGCCGGTGC" \
        + "CTCTCTCCTGGGTCCCAGGATCGGCCCCCACCATCCAGGCACGACCCCCT" \
        + "TCCCCGGCCCCTCGGCCTTTCCCCCAACTCGGCCATCTCCGACCCGGGGC" \
        + "GCGTGTTCCCCCCGGCCCGGCGCCTTCTCTCCCTCCGGGGGCACCCGCTC" \
        + "CCTAGCCCCGGCCCGGCCCTCCCCGCGGCGCAGCACGGAGTCTCGGCGTC" \
        + "CCATGGCGCAACCTACGGCCTCGGCCCAGAAGCTGGTGCGGCCGATCCGC" \
        + "GCCGTGTGCCGCATCCTGCAGATCCCGGAGTCCGACCCCTCCAACCTGCG" \
        + "GCCCTAGAGCGCCCCCGCCGCCCCGGGGGAAGGAGAGCGCGAGCGCGCTG" \
        + "AGCAGACAGAGCGGGAGAACGCGTCCTCGCCCGCCGGCCGGGAGGCCCCG" \
        + "GAGCTGGCCCATGGGGAGCAGGCGCCCGGTGCCGGCCACGACGACCGCCA" \
        + "CCGCCCGCGCCGCGACCGGCCGGTGAAGCCCAGGGACCCCCCTCTGGGAG" \
        + "AGCCCCATGAGGGCAGGAGAGTG"


@pytest.fixture(scope="session")
def ppp_five_utr(
    genome_build: GenomeBuild
) -> FiveUTRCoordinates:
    """
    5'UTR Genomic regions corresponding to one of the transcripts of the PPP1R15A-201 gene (ENSEMBL transcript ID: `ENST00000200453.6`).

    see here: https://www.ensembl.org/Homo_sapiens/Transcript/Sequence_cDNA?db=core;g=ENSG00000087074;r=19:48872421-48876058;t=ENST00000200453.
    """
    contig = genome_build.contig_by_name("8")
    assert contig is not None

    return FiveUTRCoordinates(
        regions=(
            GenomicRegion(
                contig=contig,
                start=48_872_421,
                end=48_872_651,
                strand=Strand.POSITIVE,
                ),
            GenomicRegion(
                contig=contig,
                start=48_873_225,
                end=48_873_233,
                strand=Strand.POSITIVE
                ),
        )
    )


@pytest.fixture(scope="session")
def ppp_five_utr_sequence() -> str:
    """
    5'UTR cDNA sequence of the transcript of the PPP1R15A-201 gene (ENSEMBL transcript ID: `ENST00000200453.6`) taken directly from
    the ENSEMBL website.
    
    see here: https://www.ensembl.org/Homo_sapiens/Transcript/Sequence_cDNA?db=core;g=ENSG00000087074;r=19:48872421-48876058;t=ENST00000200453.
    """
    return "GCTCTTATCGGTTCCCATCCCAGTTGTTGATCTTATGCAAGACGCTGCACGACCCCGCGC" \
        + "CCGCTTGTCGCCACGGCACTTGAGGCAGCCGGAGATACTCTGAGTTACTCGGAGCCCGAC" \
        + "GCCTGAGGGTGAGATGAACGCGCTGGCCTCCCTAACCGTCCGGACCTGTGATCGCTTCTG" \
        + "GCAGACCGAACCGGCGCTCCTGCCCCCGGGGTGACGCGCAGCTCCCAGCCGCCCAGACAC"