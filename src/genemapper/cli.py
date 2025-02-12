import argparse
import re
import logging
import os
from pathlib import Path
from typing import Dict, List
from functools import partial
from multiprocessing import Pool, cpu_count

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner

logging.basicConfig(
    filename="genome_mapping.log",
    filemode="w",
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="A combined script for GenBank metadata retrieval and CDS mapping."
    )
    subparsers = parser.add_subparsers(
        dest="command", help="Choose a subcommand: 'metadata' or 'map'.", required=True
    )
    meta_parser = subparsers.add_parser(
        "metadata", help="Fetch GenBank metadata for query and reference FASTA."
    )
    meta_parser.add_argument(
        "--query", required=True, help="Path to query FASTA file (containing CDS)."
    )
    meta_parser.add_argument(
        "--reference",
        required=True,
        help="Path to reference FASTA file (containing CDS).",
    )
    meta_parser.add_argument(
        "--email", required=True, help="Email address (required by NCBI E-utilities)."
    )
    map_parser = subparsers.add_parser(
        "map", help="Map query FASTA vs. subject FASTA and produce an Excel file."
    )
    map_parser.add_argument(
        "--query", required=True, help="Path to query FASTA file (containing CDS)."
    )
    map_parser.add_argument(
        "--subject", required=True, help="Path to subject FASTA file (containing CDS)."
    )
    map_parser.add_argument(
        "--output",
        default="queryVsSubjectMapping.xlsx",
        help="Output Excel file name (default: queryVsSubjectMapping.xlsx)",
    )
    map_parser.add_argument(
        "--similarityThreshold",
        type=float,
        default=0.99,
        help="Similarity threshold for matching (default: 0.99).",
    )
    map_parser.add_argument(
        "--numProcesses",
        type=int,
        default=None,
        help="Number of parallel processes to use (default: all CPUs).",
    )
    return parser.parse_args()


def parse_accession_from_fasta(fasta_path: Path) -> str:
    """
    Read the first record in a FASTA, parse out the accession from its ID or description.
    Returns None if unable to find an accession.
    Example FASTA header line:
      >lcl|CP006763.1_cds_AGY74238.1_1 [locus_tag=CAETHG_0001] ...
    We want to extract 'CP006763.1'.
    """
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if not records:
        return None
    first_id = records[0].id
    match = re.match(r"^lcl\|([^_]+)_", first_id)
    if match:
        return match.group(1)
    return None


def fetch_genbank_record(accession: str, email: str):
    """
    Try to fetch a GenBank record from NCBI by accession.
    Return the SeqRecord if successful, or None if not.
    """
    Entrez.email = email
    try:
        print(f"  Fetching GenBank record for accession: {accession} ...")
        with Entrez.efetch(
            db="nucleotide", id=accession, rettype="gb", retmode="text"
        ) as handle:
            record = SeqIO.read(handle, "genbank")
        return record
    except Exception as e:
        print(f"  [!] Could not fetch record for {accession}. Error: {e}")
        return None


def print_local_fasta_info(fasta_path: str):
    """
    Print basic info about the local FASTA file:
    file name, number of sequences, first sequence ID, etc.
    """
    seqs = list(SeqIO.parse(str(fasta_path), "fasta"))
    num_seqs = len(seqs)
    print(f"  File: {fasta_path}")
    print(f"  Number of Sequences: {num_seqs}")
    if seqs:
        print(f"  Example ID: {seqs[0].id}")
        print(f"  Example Description: {seqs[0].description}")
        print(f"  Example Sequence Length: {len(seqs[0].seq)}")


def print_genbank_record_info(record):
    """
    Print some basic info about a GenBank record.
    """
    print(f"    ID:          {record.id}")
    print(f"    Name:        {record.name}")
    print(f"    Description: {record.description}")
    print(f"    Organism:    {record.annotations.get('organism', 'Unknown')}")
    print(f"    Date:        {record.annotations.get('date', 'Unknown')}")


def do_metadata(args):
    """
    Handle the 'metadata' subcommand.
    """
    print("\n[Query FASTA Info]")
    print_local_fasta_info(args.query)
    query_acc = parse_accession_from_fasta(Path(args.query))
    if query_acc:
        print(f"  Accession found: {query_acc}")
        query_record = fetch_genbank_record(query_acc, args.email)
        if query_record:
            print("  Successfully retrieved GenBank record:")
            print_genbank_record_info(query_record)
        else:
            print("  [!] Could not retrieve metadata from GenBank.")
    else:
        print("  [!] No valid accession found in the FASTA header.")
    print("\n[Reference FASTA Info]")
    print_local_fasta_info(args.reference)
    ref_acc = parse_accession_from_fasta(Path(args.reference))
    if ref_acc:
        print(f"  Accession found: {ref_acc}")
        ref_record = fetch_genbank_record(ref_acc, args.email)
        if ref_record:
            print("  Successfully retrieved GenBank record:")
            print_genbank_record_info(ref_record)
        else:
            print("  [!] Could not retrieve metadata from GenBank.")
    else:
        print("  [!] No valid accession found in the FASTA header.")


def parseFastaFile(filePath: str) -> Dict[str, Dict[str, str]]:
    """
    Parse a FASTA file to extract sequences along with their sequence coordinates.
    Returns a dictionary keyed by locus tag, with 'sequence' and 'location' as sub-keys.
    """
    if not os.path.exists(filePath):
        raise FileNotFoundError(f"File not found: {filePath}")
    parsedSequences = {}
    try:
        for record in SeqIO.parse(filePath, "fasta"):
            locusTagMatch = record.description.find("[locus_tag=")
            locationMatch = record.description.find("[location=")
            if locusTagMatch != -1:
                endTag = record.description.find("]", locusTagMatch)
                locusTag = record.description[locusTagMatch + 11 : endTag]
                location = "Unknown"
                if locationMatch != -1:
                    endLocation = record.description.find("]", locationMatch)
                    location = record.description[locationMatch + 10 : endLocation]
                parsedSequences[locusTag] = {
                    "sequence": str(record.seq),
                    "location": location,
                }
    except Exception as e:
        logging.error(f"Error parsing FASTA file {filePath}: {e}")
        return {}
    return parsedSequences


def alignAndScore(
    querySequence: str,
    subjectSequence: str,
    similarityThreshold: float = 0.0,
    aligner: PairwiseAligner = None,
) -> float:
    """
    Perform alignment and calculate a normalized score using PairwiseAligner.
    """
    if aligner is None:
        aligner = PairwiseAligner()
        aligner.mode = "global"
    alignments = aligner.align(querySequence, subjectSequence)
    bestAlignment = alignments[0]
    score = bestAlignment.score / max(len(querySequence), len(subjectSequence))
    if similarityThreshold > 0 and score < similarityThreshold:
        return 0.0
    return score


def findBestMatchesWithLocation(
    queryTag: str,
    querySeq: str,
    queryLoc: str,
    subjectGenome: Dict[str, Dict[str, str]],
    similarityThreshold: float = 0.99,
    aligner: PairwiseAligner = None,
) -> List[Dict]:
    """
    Given a single query gene, find the best match(es) above the similarityThreshold.
    """
    if aligner is None:
        aligner = PairwiseAligner()
        aligner.mode = "global"
    scoreList = []
    for subjTag, subjInfo in subjectGenome.items():
        rawScore = alignAndScore(querySeq, subjInfo["sequence"], 0.0, aligner)
        scoreList.append(
            (subjTag, subjInfo["sequence"], subjInfo["location"], rawScore)
        )
    if not scoreList:
        logging.warning(f"Subject genome is empty. No matches for {queryTag}.")
        return [
            {
                "queryLocusTag": queryTag,
                "querySequence": querySeq,
                "queryLocation": queryLoc,
                "subjectLocusTag": None,
                "subjectSequence": None,
                "subjectLocation": None,
                "similarityScore": 0.0,
                "similarityThresholdPassed": False,
                "tiedForBest": False,
            }
        ]
    bestScore = max(s[3] for s in scoreList)
    topHits = [s for s in scoreList if s[3] == bestScore]
    if bestScore < similarityThreshold:
        logging.warning(
            f"No match above threshold for {queryTag}. Best score={bestScore:.4f}"
        )
        return [
            {
                "queryLocusTag": queryTag,
                "querySequence": querySeq,
                "queryLocation": queryLoc,
                "subjectLocusTag": None,
                "subjectSequence": None,
                "subjectLocation": None,
                "similarityScore": bestScore,
                "similarityThresholdPassed": False,
                "tiedForBest": False,
            }
        ]
    tieForBest = len(topHits) > 1
    if tieForBest:
        logging.warning(
            f"Query gene {queryTag} has a tie for best match. Score={bestScore:.4f}. "
            f"{len(topHits)} subject genes share that top score."
        )
    results = []
    for subjTag, subjSeq, subjLoc, scoreVal in topHits:
        results.append(
            {
                "queryLocusTag": queryTag,
                "querySequence": querySeq,
                "queryLocation": queryLoc,
                "subjectLocusTag": subjTag,
                "subjectSequence": subjSeq,
                "subjectLocation": subjLoc,
                "similarityScore": scoreVal,
                "similarityThresholdPassed": True,
                "tiedForBest": tieForBest,
            }
        )
    return results


def parallelGenomeMappingWithLocation(
    queryFilePath: str,
    subjectFilePath: str,
    outputFile: str = "queryVsSubjectMapping.xlsx",
    numProcesses: int = None,
    similarityThreshold: float = 0.99,
):
    """
    Perform parallel genome sequence mapping with location information. Results are saved to Excel.
    """
    numProcesses = numProcesses or cpu_count()
    logging.info("Parsing FASTA files...")
    queryGenome = parseFastaFile(queryFilePath)
    subjectGenome = parseFastaFile(subjectFilePath)
    logging.info(f"Query genome has {len(queryGenome)} sequences.")
    logging.info(f"Subject genome has {len(subjectGenome)} sequences.")
    logging.info("Performing parallel sequence comparisons...")
    aligner = PairwiseAligner()
    aligner.mode = "global"
    with Pool(processes=numProcesses) as pool:
        partialMatchFunc = partial(
            findBestMatchesWithLocation,
            subjectGenome=subjectGenome,
            similarityThreshold=similarityThreshold,
            aligner=aligner,
        )
        mappingResults = pool.starmap(
            partialMatchFunc,
            [
                (tag, data["sequence"], data["location"])
                for tag, data in queryGenome.items()
            ],
        )
    flatResults = []
    for subList in mappingResults:
        flatResults.extend(subList)
    mappingDf = pd.DataFrame(flatResults)
    mappingDf.to_excel(outputFile, index=False)
    logging.info(f"Mapping completed: {len(mappingDf)} rows generated.")
    logging.info(f"Results saved to {outputFile}.")
    noMatchCount = (mappingDf["similarityThresholdPassed"] == False).sum()
    tieCount = (mappingDf["tiedForBest"] == True).sum()
    totalQueryCount = len(queryGenome)
    logging.info(
        f"Out of {totalQueryCount} query genes, {noMatchCount} had no match >= threshold."
    )
    logging.info(f"Number of rows with a tie for best match: {tieCount}.")
    matchedQueryTags = set(
        mappingDf.loc[mappingDf["similarityThresholdPassed"], "queryLocusTag"].dropna()
    )
    logging.info(
        f"Distinct query genes matched above threshold: {len(matchedQueryTags)}"
    )
    matchedSubjectTags = set(
        mappingDf.loc[
            mappingDf["similarityThresholdPassed"], "subjectLocusTag"
        ].dropna()
    )
    allSubjectTags = set(subjectGenome.keys())
    unmatchedSubjectTags = allSubjectTags - matchedSubjectTags
    logging.info(
        f"Number of subject genes unmatched above threshold: {len(unmatchedSubjectTags)}"
    )


def do_mapping(args):
    """
    Handle the 'map' subcommand.
    """
    print("\nStarting genome mapping ...")
    print(f"  Query file:   {args.query}")
    print(f"  Subject file: {args.subject}")
    print(f"  Output file:  {args.output}")
    print(f"  Similarity threshold: {args.similarityThreshold}")
    if args.numProcesses:
        print(f"  Using {args.numProcesses} processes ...")
    else:
        print(f"  Using all available CPU cores ...")
    parallelGenomeMappingWithLocation(
        queryFilePath=args.query,
        subjectFilePath=args.subject,
        outputFile=args.output,
        numProcesses=args.numProcesses,
        similarityThreshold=args.similarityThreshold,
    )
    print("Mapping complete. See 'genome_mapping.log' for details.")


def main():
    args = parse_arguments()
    if args.command == "metadata":
        do_metadata(args)
    elif args.command == "map":
        do_mapping(args)


if __name__ == "__main__":
    main()
