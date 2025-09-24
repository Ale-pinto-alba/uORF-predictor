import requests

from gpsea.model.genome import Region, Strand, GenomicRegion
from uorf_predictor.instances import FiveUTRCoordinates, UORFCoordinates


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


def obtain_uorf_in_five_utr(five_utrs: FiveUTRCoordinates, start_uorf: int, end_uorf: int) -> UORFCoordinates:
    """
    Map the genomic coordinates of the uORF into the 5'UTR of the transcript to obtain the relative position within the sequence.

    :param five_utrs: 5'UTR Genomic regions of the transcript.
    :param start_uorf: Start position of the uORF.
    :param end_uorf: End position of the uORF.
    """
    five_utrs_tuple = [(region.start, region.end) for region in five_utrs.regions]
    gene_strand = five_utrs.regions[0].strand
    contig = five_utrs.regions[0].contig

    uorf_length = end_uorf - start_uorf
    variant_cdna_pos = None

    if gene_strand == Strand.POSITIVE:
        cdna_pos = 0
        for start, end in sorted(five_utrs_tuple):
            five_utr_region_length = end - start
            if start <= start_uorf <= end:
                variant_cdna_pos = cdna_pos + (start_uorf - start + 1)
            cdna_pos += five_utr_region_length
        
        if variant_cdna_pos is not None:
            if end_uorf >= five_utrs_tuple[-1][1]:
                ouorf = True
            else:
                ouorf = False
            return UORFCoordinates(
                        five_utr= five_utrs,
                        uorf= Region(start= variant_cdna_pos, end= variant_cdna_pos + uorf_length),
                        ouorf= ouorf,
                    )
    else:
        five_utrs_tuples_sorted_reversed = sorted(five_utrs_tuple, reverse=True)
        uorf_genomic_region = GenomicRegion(contig=contig, start=start_uorf, end=end_uorf, strand=Strand.POSITIVE)
        cdna_pos = 0
        for start, end in five_utrs_tuples_sorted_reversed:
            five_utr_region_length = end - start
            if start <= uorf_genomic_region.start_on_strand(gene_strand) <= end:
                    variant_cdna_pos = cdna_pos + (end - uorf_genomic_region.start_on_strand(gene_strand))
            cdna_pos += five_utr_region_length

        if variant_cdna_pos is not None:
            relative_position = cdna_pos - variant_cdna_pos
            if end_uorf <= five_utrs_tuple[-1][0]:
                ouorf = True
            else:
                ouorf = False
            return UORFCoordinates(
                five_utr=five_utrs,
                uorf=Region(start= relative_position, end= relative_position + uorf_length + 1),
                ouorf=ouorf,
            )   
            
def check_start_and_stop_codon(uorf_sequence: str) -> bool: 
    """
    Check if the uORF is correctly framed by a start and a stop codon.

    :param uorf_sequence: `str` containing the uORF cDNA sequence.
    """
    start_codons = ["ATG", "CTG", "GTG", "TTG", "ACG"]
    stop_codons = ["TAG", "TAA", "TGA"]

    return (
        any(uorf_sequence.startswith(codon) for codon in start_codons) and
        any(uorf_sequence.endswith(codon) for codon in stop_codons)
    )

import pyBigWig
import numpy as np

def get_mean_phastcons_phylop(bigwig_path, chrom, start, end):
    try:
        with pyBigWig.open(bigwig_path) as bw:
            if chrom not in bw.chroms():
                raise ValueError(f"Chromosome '{chrom}' not in the BigWig file.")

            scores = bw.values(chrom, start, end, numpy=True)
            if scores is None:
                return np.nan
            
            scores_filtered = scores[~np.isnan(scores)]
            return np.mean(scores_filtered) if len(scores_filtered) > 0 else np.nan
    except Exception as e:
        return np.nan