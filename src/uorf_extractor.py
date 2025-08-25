import typing

import requests

from gpsea.model.genome import Region
from instances import FiveUTRCoordinates, UORFCoordinates


def fetch_cdna_from_ensembl(transcript_id: str, timeout: float = 30.,) -> str:
    """
    Download cDNA sequence (spliced mRNA) for a given transcript from Ensembl's REST API.

    :param transcript_id: Ensembl transcript identifier e.g. `ENST00000381418`
    """
    base_url_cdna = f"https://rest.ensembl.org/sequence/id/{transcript_id}?type=cdna&content-type=text/x-fasta"
    
    response = requests.get(base_url_cdna, timeout=timeout)

    if response.status_code == 200:
        lines = response.text.splitlines()
        return ''.join(lines[1:])
    else:
        response.raise_for_status()


def get_five_prime_sequence(cdna_sequence: str, five_utrs: FiveUTRCoordinates) -> str:
    """
    Return the 5'UTR cDNA sequence of a given transcript nucleotide sequence (spliced mRNA).

    :param transcript_sequence: transcript nucleotide sequence.
    :param five_utrs: 5'UTR Genomic Region(s).
    """
    return cdna_sequence[:len(five_utrs)]


def uorf_extractor(five_utr: FiveUTRCoordinates, five_sequence: str) -> typing.Collection[UORFCoordinates]:
    """
    Take the cDNA nucleotide sequence of a transcript 5'UTR region to extract the uORFs sequences available.

    :param five_utr: list of Genomic Regions corresponding to the 5'UTRs regions.
    :param five_sequence: 5'UTR cDNA sequence.
    """
    uorfs = []
    start_codons = ["ATG", "CTG", "GTG", "ACG", "TTG"]
    start_position = 0

    while start_position < len(five_sequence) - 2:
        codon = five_sequence[start_position:start_position + 3]
        if codon in start_codons:
            found_stop = False
            start_codon = codon

            for i in range(start_position + 3, len(five_sequence) - 2, 3):  
                stop_codon = five_sequence[i:i + 3]
                if stop_codon in ["TAA", "TAG", "TGA"]:
                    stop_index = i + 3
                    found_stop = True
                    break

            if found_stop:
                uorfs.append(UORFCoordinates(
                    five_utr=five_utr,
                    uorf=Region(start=start_position, end=stop_index),
                    ouorf= False,
                    start_codon= start_codon,
                ))
                start_position = stop_index  
            else:
                uorfs.append(UORFCoordinates(
                    five_utr=five_utr,
                    uorf=Region(start=start_position, end=len(five_sequence)),
                    ouorf= True,
                    start_codon= start_codon,
                ))
                start_position = start_position + 1  

    return uorfs