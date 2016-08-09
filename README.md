[candYgene](http://software.esciencecenter.nl/project/candygene/)
=
*SIGA.py* is a command-line tool to generate *Semantically Interoperable Genome Annotations* from
[GFF](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) files according to the Resource Description Framework ([RDF](https://www.w3.org/TR/rdf11-concepts/)) specification.

**Key features**
- process multiple input files in GFF (versions 2 and 3)
- genome annotations (features) stored in [SQLite](https://sqlite.org/) database and serialized as RDF graph(s) in plain text formats:
  - [XML](https://www.w3.org/TR/rdf-syntax-grammar/)
  - [N-Triples](https://www.w3.org/TR/n-triples/)
  - [Turtle](https://www.w3.org/TeamSubmission/turtle/)
  - [Notation3](https://www.w3.org/DesignIssues/Notation3.html) (N3)
- (hierarchy of) feature types supported: *chromosome -> gene -> mRNA -> [CDS, exon, intro, five_prime_UTR, three_prime_UTR]*
- features and their parent-child relations mapped to [SO(FA)](http://www.sequenceontology.org/) and [FALDO](https://github.com/JervenBolleman/FALDO) ontologies
- parent-child feature relationships checked for referential integrity

**Installation**

`virtualenv sigaenv`

`source sigaenv/bin/activate`

`pip install -r requirements.txt`

**Example data**

`cd examples`

or download the latest tomato genome annotations (ITAG2.4 release) in GFF files from the [Sol Genomics Network](https://solgenomics.net) (SGN).

`wget ftp://ftp.solgenomics.net/genomes/Solanum_lycopersicum/annotation/ITAG2.4_release/ITAG2.4_gene_models.gff3`

**Example usage**

`cd src`

Generate RDF triples in Turtle format (default) from a GFF file (two-steps):

`python SIGA.py db -cV ITAG2.4_gene_models.gff3` # GFF->DB

`python SIGA.py rdf -b https://solgenomics.net ITAG2.4_gene_models.db` # DB->RDF

Summary of I/O files:

`ITAG2.4_gene_models.gff3` # GFF file

`ITAG2.4_gene_models.db`   # SQLite database

`ITAG2.4_gene_models.ttl`  # RDF file in Turtle

**Import RDF graph into Virtuoso**

See the [documentation](http://virtuoso.openlinksw.com/dataspace/doc/dav/wiki/Main/VirtBulkRDFLoader) on bulk data loading.

Edit _virtuoso.ini_ config file by adding _/mydatadir/_ to _DirsAllowed_.

Connect to db server:
`isql 1111 dba`

Delete (old) RDF graph if necessary:
`SPARQL CLEAR GRAPH "https://solgenomics.net#"`

Delete any previously registered data files:
`DELETE FROM DB.DBA.load_list;`

Register data file(s):
`ld_dir('/mydatadir/', 'ITAG2.4_gene_models.ttl', 'https://solgenomics.net#');`

List registered data file(s):
`SELECT * FROM DB.DBA.load_list;`

Bulk data loading:
`rdf_loader_run();`

Note: For loading a single data file one could use the following command:

`SPARQL LOAD "file:///mydatadir/features.ttl" INTO "https://solgenomics.net#";`

However, this approach results in additional triples (generated by Virtuoso) which are not present in the input file.

Count imported triples:
`SPARQL SELECT COUNT(*) FROM <https://solgenomics.net#> { ?s ?p ?o };`

**Import RDF graph into Berkeley DB**

Use Redland RDF processor:

`rdfproc ITAG2.4_gene_models parse ITAG2.4_gene_models.ttl turtle`

Count triples in the input file:
`rapper -i turtle -c ITAG2.4_gene_models.ttl`
