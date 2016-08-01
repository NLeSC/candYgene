[candYgene](http://software.esciencecenter.nl/project/candygene/)
=========

*SIGA.py* is a command-line tool to generate *Semantically Interoperable Genome Annotations* from
[GFF](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) files according to the Resource Description Framework ([RDF](https://www.w3.org/TR/rdf11-concepts/)) specification.

**Key features**
- process multiple input files in GFF (versions 2 and 3)
- parsed genome annotations stored in a relational schema ([SQLite](https://sqlite.org/)) and RDF graph(s) in one of the supported serializations (plain text formats):
  - [XML](https://www.w3.org/TR/rdf-syntax-grammar/)
  - [N-Triples](https://www.w3.org/TR/n-triples/)
  - [Turtle](https://www.w3.org/TeamSubmission/turtle/)
  - [Notation3](https://www.w3.org/DesignIssues/Notation3.html) (N3)
- currently supported (hierarchy of) features: *gene -> mRNA -> [CDS, exon, intro, five_prime_UTR, three_prime_UTR]*
- typed features and their parent-child relations using [SO(FA)](http://www.sequenceontology.org/) and [FALDO](https://github.com/JervenBolleman/FALDO) ontologies
- ensure referential integrity of data (i.e., parent-child feature relations)

**Installation**

`virtualenv sigaenv`

`source sigaenv/bin/activate`

`pip install -r requirements.txt`

**Example data**

Small GFF files are in `cd examples`

or download tomato genome annotations (ITAG2.4 release) in GFF from the [Sol Genomics Network](https://solgenomics.net) (SGN).

`wget ftp://ftp.solgenomics.net/genomes/Solanum_lycopersicum/annotation/ITAG2.4_release/ITAG2.4_gene_models.gff3`

**Example usage**

Check data integrity (verbose output) and generate RDF triples in Turtle format (default) at base URI

`cd src`

`python SIGA.py -cV ITAG2.4_gene_models.gff3 -b https://solgenomics.net/`

Output files in current working directory:

`ITAG2.4_gene_models.db` # relational database in SQLite

`ITAG2.4_gene_models.ttl` # RDF triples in Turtle

**Import data into Virtuoso RDF Store**

Edit _virtuoso.ini_ config file by adding _/mypath-to-RDF/_ to _DirsAllowed_.

Connect to DB:
`isql 1111 dba`

Register data set:
`ld_dir('/mypath-to-RDF/', 'ITAG2.4_gene_models.ttl', 'https://solgenomics.net#');`

Check if data set is registered:
`SELECT * FROM DB.DBA.load_list;`

Load RDF file:
`rdf_loader_run();`

Count imported triples:
`SPARQL SELECT COUNT(*) FROM <https://solgenomics.net#> { ?s ?p ?o };`
