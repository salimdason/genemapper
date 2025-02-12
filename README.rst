GeneMapper
==========

|PyPI version| |Build Status| |Documentation Status|

A command-line tool for retrieving GenBank metadata and mapping coding
sequences between query and subject FASTA files to identify their most
similar counterparts.

-  **Free software**: MIT license
-  **Documentation**: https://genemapper.readthedocs.io

Features
--------

-  Automated retrieval of GenBank metadata for query and reference
   sequences files in fasta format.
-  Parallelized mapping of query CDS against subject CDS with scoring
   and threshold-based similarity.
-  Detailed log and Excel output of mapping results
-  Adjustable similarity thresholds and multi-processor usage.

Installation
------------

.. code:: bash

   pip install genemapper

Usage - Metadata Retrieval
------------

.. code:: bash

    genemapper metadata --query my_query_gene_cds.fasta --reference my_subject_gene_cds.fasta --email mohammed.dason@polito.it

Usage - Gene Mapping
------------

.. code:: bash

    genemapper map --query my_query_gene_cds.fasta --subject my_subject_gene_cds.fasta --output LAbriniVersusCAETHG_Mapping.xlsx --similarityThreshold 0.99 --numProcesses 20

.. |PyPI version| image:: https://img.shields.io/pypi/v/genemapper.svg
   :target: https://pypi.python.org/pypi/genemapper
.. |Build Status| image:: https://img.shields.io/travis/salimdason/genemapper.svg
   :target: https://travis-ci.com/salimdason/genemapper
.. |Documentation Status| image:: https://readthedocs.org/projects/genemapper/badge/?version=latest
   :target: https://genemapper.readthedocs.io/en/latest/?version=latest
