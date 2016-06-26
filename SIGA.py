"""
SIGA.py is a command-line tool to generate Semantically Interoperable Genome Annotations from
GFF files [1] according to the RDF specification [2].

References:
[1] Generic Feature Format specification, http://www.sequenceontology.org/gff3.shtml
[2] Resource Description Framework, https://www.w3.org/TR/rdf11-concepts/

Usage:
  SIGA.py -h | --help
  SIGA.py -v | --version
  SIGA.py [-cV ] [ -d DB_FILE | -e DB_FILE_EXT ] [ -o FORMAT ] GFF_FILE...

Arguments:
  GFF_FILE...      Input file(s) in GFF versions 2 or 3.

Options:
  -h, --help
  -v, --version
  -V, --verbose    Use for debugging.
  -d DB_FILE       Populate a single SQLite database from one or more GFF files.
  -e DB_FILE_EXT   Database file extension. [default: .db].
  -o FORMAT        Select RDF serialization format: turtle (.ttl), n3 (.n3) or xml (.rdf). [default: turtle]
  -c               Check the referential integrity of database(s).

"""

from __future__ import print_function
from docopt import docopt
from rdflib import Graph, URIRef, Literal, BNode
from rdflib.namespace import Namespace, RDF, RDFS, XSD

import os
import gffutils as gff
import sqlite3 as sql

__author__  = 'Arnold Kuzniar'
__version__ = '0.1.3'
__status__  = 'Prototype'
__license__ = 'Apache License, Version 2.0'


def normalize_filext(s):
    dot = '.'
    if s.startswith(dot) is False:
        s = dot + s
    return s

def normalize_feature_id(id):
    return id.split(':')[1] # gene:Solyc00g005000.2 -> Solyc00g005000.2

def triplify(db, format):
    format2fext = dict(turtle = '.ttl', xml = '.rdf', n3 = '.n3')
    assert(format in format2fext), "Unsupported RDF serialization '%s'." % format

    # define additional namespaces
    SO = Namespace('http://purl.obolibrary.org/obo/so.owl#') # Sequence Ontology Feature Annotation
    FALDO = Namespace('http://biohackathon.org/resource/faldo#')

    # bind prefixes to namespaces
    g = Graph()
    g.bind('so', SO)
    g.bind('faldo', FALDO)

    # map feature types and strandedness to ontology terms
    feature2onto = {'gene' : SO.SO_0000704,
                    'mRNA' : SO.SO_0000234,
                    'CDS'  : SO.SO_0000316,
                    'exon' : SO.SO_0000147,
                    'intron' : SO.SO_0000188,
                    'five_prime_UTR' : SO.SO_0000204,
                    'three_prime_UTR' : SO.SO_0000205,
                    '+' : FALDO.ForwardStrandPosition,
                    '-' : FALDO.ReverseStrandPosition,
                    '?' : FALDO.StrandedPosition,
                    '.' : FALDO.Position}

    # Note: so far there is no optimal/universal solution to resolve features in SGN:
    # e.g., a gene feature identified by 'gene:Solyc00g005000.2' resolves to
    # https://solgenomics.net/locus/Solyc00g005000.2/view
    # However, its related features such as mRNA, exon or intron do not resolve the same way
    # e.g., mRNA:Solyc00g005000.2.1 resolves to https://solgenomics.net/feature/17660840/details
    base_uri ='https://solgenomics.net/locus/'

    for feature in db.all_features():
        if feature.strand not in feature2onto:
            raise KeyError("Incorrect strand information for feature ID '%s'." % feature.id)
        try:
            feature_type = URIRef(feature2onto[feature.featuretype])
            strand = feature2onto[feature.strand]
            subject = URIRef(base_uri + normalize_feature_id(feature.id) + '/view')
            start = BNode()
            end = BNode()

            # add triples to graph
            g.add( (subject, RDF.type, feature_type) )
            g.add( (subject, RDFS.label, Literal(feature.featuretype, datatype=XSD.string)) )
            g.add( (subject, RDF.type, FALDO.Region) )

            # add feature start/end positions and strand info
            g.add( (subject, FALDO.begin, start) )
            g.add( (start, RDF.type, FALDO.ExactPosition) )
            g.add( (start, RDF.type, strand) )
            g.add( (start, FALDO.position, Literal(feature.start, datatype=XSD.nonNegativeInteger)) )
            g.add( (start, FALDO.reference, Literal(feature.seqid, datatype=XSD.string)) )
            g.add( (subject, FALDO.end, end) )
            g.add( (end, RDF.type, FALDO.ExactPosition) )
            g.add( (end, RDF.type, strand) )
            g.add( (end, FALDO.position, Literal(feature.end, datatype=XSD.nonNegativeInteger)) )
            g.add( (end, FALDO.reference, Literal(feature.seqid, datatype=XSD.string)) )
        except KeyError:
            pass

    outfile = os.path.splitext(db.dbfn)[0] + format2fext[format]
    with open(outfile, 'w') as fout:
        fout.write(g.serialize(format=format))
 
if __name__ == '__main__':
    args = docopt(__doc__, version=__version__)
    #print(args)

    format = args['-o']
    debug = args['--verbose']
    fk_constraints = 'ON' if args['-c'] is True else 'OFF'
    pragmas = dict(foreign_keys=fk_constraints)

    # loop through GFF files, populate databases and write RDF graphs
    for gff_file in args['GFF_FILE']:
        if args['-d']:
            # populate a single db from all GFF files
            db_file = args['-d']
            if os.path.exists(db_file):
                db = gff.FeatureDB(db_file)
                db.update(gff_file)
            else:
                db = gff.create_db(gff_file, db_file, verbose=debug, pragmas=pragmas, force=False)
        else:
            # populate one db per GFF file
            base_name = os.path.splitext(gff_file)[0]
            db_file = base_name + normalize_filext(args['-e'])
            try:
                db = gff.create_db(gff_file, db_file, verbose=debug, pragmas=pragmas, force=False)
            except sql.OperationalError:
                raise IOError("Database file '%s' already exists." % db_file)
            except ValueError:
                raise IOError("GFF file '%s' not found." % gff_file)
            except sql.IntegrityError, e:
                raise IOError("%s in database '%s'." % (e, db_file))
            triplify(db, format)

    if args['-d']:
        db_file = args['-d']
        db = gff.FeatureDB(db_file)
        triplify(db, format)
