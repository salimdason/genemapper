GeneMapper
==========

|PyPI version| |Documentation Status| |BioPython| |Pandas| |Python version|

A command-line tool for retrieving GenBank metadata and mapping coding
sequences between query and subject FASTA files to identify their most
similar counterparts.

- **Free software**: MIT license
- **Documentation**: https://genemapper.readthedocs.io

Overview
--------

GeneMapper provides two main functions:

1. **Metadata Retrieval**:
   - Extracts an accession number from a FASTA header and fetches the corresponding GenBank record from NCBI.
   Prints out relevant details such as organism name, record description, and publication date.

2. **Mapping Sequences**:
   - Reads a **query FASTA** and a **subject FASTA**, each containing coding sequences (CDS). 
   It uses pairwise alignment to score each query sequence against possible subject matches, 
   and outputs an Excel file listing best matches, tie scenarios, unmatched sequences e.t.c.

Features
--------

- Automated retrieval of GenBank metadata for query and reference
  sequences files in FASTA format.
- Parallelized mapping of query CDS against subject CDS with scoring
  and threshold-based similarity.
- Detailed log and Excel output of mapping results
- Adjustable similarity thresholds and multi-processor usage.

Installation
------------

.. code:: bash

   pip install genemapper

If installing from source, navigate to the package directory and run:

.. code:: bash

   pip install .

Usage
-----

Below are the two primary subcommands: **metadata** and **map**.

Usage - Metadata Retrieval
--------------------------

Fetch GenBank metadata for query/reference FASTA files. The FASTA headers
must contain an accession in the form ">lcl|ACCESSION_cds_..." for
automatic extraction.

.. code:: bash

    genemapper metadata --query my_query_gene_cds.fasta \
                        --reference my_subject_gene_cds.fasta \
                        --email mohammed.dason@polito.it

What this does:

- **Parses** your query FASTA to extract an accession.
- **Fetches** the corresponding GenBank record from NCBI, printing details like:
  organism name, record date, and a brief description.
- **Repeats** for the reference FASTA.

Usage - Gene Mapping
--------------------

Compare a query FASTA to a subject FASTA, finding the best-matching CDS:

.. code:: bash

    genemapper map --query my_query_gene_cds.fasta \
                   --subject my_subject_gene_cds.fasta \
                   --output results.xlsx \
                   --similarityThreshold 0.99 \
                   --numProcesses 20

What this does:

- **Parses** both FASTA files, extracting locus tags and sequences.
- **Aligns** each query sequence to all subject sequences using a global aligner.
- **Filters** results by the similarity threshold.
- **Saves** everything in an Excel file ("results.xlsx").

How It Works
------------

1. **parse_accession_from_fasta**  
   Looks at the first record in a FASTA, uses a regular expression
   to detect the accession in the header.

2. **fetch_genbank_record**  
   Contacts NCBI using Biopython’s "Entrez.efetch" to download the GenBank
   record. Requires a valid email to respect NCBI's policy.

3. **parseFastaFile**  
   Reads all CDS in a given FASTA, extracting a "locus_tag" (and optional
   "[location=...]") to store each sequence in a python dictionary.

4. **alignAndScore**  
   Performs a global alignment with Biopython’s PairwiseAligner. Normalizes
   the score by the length of the longer sequence, returning a floating-point
   similarity value.

5. **findBestMatchesWithLocation**  
   Loops over all subject sequences, keeps track of the highest-scoring match,
   and flags ties if multiple subjects share the same top score.

6. **parallelGenomeMappingWithLocation**  
   Uses Python's multiprocessing to handle each query locus tag in parallel,
   which speeds up large genome comparisons. Writes results to Excel when finished.

Logging
-------

GeneMapper writes runtime events and warnings to a log file called
"genome_mapping.log". This includes:

- **Metadata retrieval failures** (e.g., no valid accession).
- **Mapping** details, warnings about no matches above threshold, ties, or I/O errors.
- **Summary** of how many query genes matched, remained unmatched, etc.

Check this log file if something goes wrong or for deeper insight into the
script’s operations.

Troubleshooting & Tips
----------------------

1. **No Accession Found**  
   Make sure your FASTA headers follow the format ">lcl|ACCESSION_cds_...".

2. **No Matches Above Threshold**  
   Lower the "--similarityThreshold" or check for sequence divergence. Review
   the log to see actual scores.

3. **Parallel Performance**  
   If you have many CPU cores, increase "--numProcesses". NOTE: by deafault all cores are used.
   If you hit resource limits, please reduce it by setting a suitable value.

4. **Spreadsheet Issues**  
   The script writes an Excel file via **pandas** and **openpyxl**. If needed,
   you can adapt the code to write CSV by modifying "mappingDf.to_csv(...)".

5. **Versions**  
   - Biopython ≥ 1.78 is required for "PairwiseAligner".
   - Python 3.7+ is recommended.
   - Openpyxl 3.1.5 (pinned) is required for the creation of the excel file
   - Pandas > 2.0 is recommended. 

License
-------

**MIT License**: See the licence file for full details.

.. |PyPI version| image:: https://img.shields.io/pypi/v/genemapper.svg
   :target: https://pypi.python.org/pypi/genemapper
   :alt: PyPI version

.. |Documentation Status| image:: https://readthedocs.org/projects/genemapper/badge/?version=latest
   :target: https://genemapper.readthedocs.io/en/latest/?version=latest
   :alt: Documentation Status

.. |BioPython| image:: https://img.shields.io/badge/BioPython-1.78-brightgreen
   :target: https://biopython.org
   :alt: BioPython

.. |Pandas| image:: https://img.shields.io/badge/Pandas-%3E%3D1.0.0-blue
   :target: https://pandas.pydata.org
   :alt: Pandas

.. |Python version| image:: https://img.shields.io/pypi/pyversions/genemapper.svg
   :target: https://www.python.org/
   :alt: Python version
