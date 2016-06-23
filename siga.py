"""
Generate semantically interoperable genome annotations from GFF[1] files according to the RDF[2] specification.

References:
[1] Generic Feature Format specification, http://www.sequenceontology.org/gff3.shtml
[2] Resource Description Framework, https://www.w3.org/TR/rdf11-concepts/

Usage:
  siga.py -h | --help
  siga.py -v | --version
  siga.py [-cV ] [ -d DB_FILE | -e DB_FILE_EXT ] [ -o FORMAT ] GFF_FILE...

Arguments:
  GFF_FILE...      Input file(s) in GFF (versions 2/3).

Options:
  -h, --help
  -v, --version
  -V, --verbose    Use for debugging.
  -d DB_FILE       Populate GFF database(s) in SQLite.
  -e DB_FILE_EXT   Database file extension. [default: .db].
  -o FORMAT        Select RDF graph serialization: turtle (.ttl), n3 (.n3) or xml (.rdf). [default: turtle]
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
__version__ = '0.1.2'
__status__  = 'Prototype'
__license__ = 'Apache License, Version 2.0'


def normalize_filext(s):
    dot = '.'
    if s.startswith(dot) is False:
        s = dot + s
    return s

#def normalize_feature_id(id):
#    return replace(id + ':', '')

def triplify(db, format):
    format2fext = dict(turtle = '.ttl', xml = '.rdf', n3 = '.n3')
    assert(format in format2fext), "Unsupported RDF serialization '%s'." % format

    # define additional namespaces
    SO = Namespace('http://purl.obolibrary.org/obo/') # Sequence Ontology (Feature Annotation)
    FALDO = Namespace('http://biohackathon.org/resource/faldo#')
    BASE = Namespace('https://solgenomics.net/search/quick?term=') # so far no better solution to resolve features

    # bind prefixes to namespaces
    g = Graph()
    g.bind('so', SO)
    g.bind('faldo', FALDO)
    g.bind(None, BASE)

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

    for feature in db.all_features():
        try:
            s = URIRef(BASE + feature.id)
            p = URIRef(RDF.type)
            o = URIRef(feature2onto[feature.featuretype])
            start = BNode()
            end = BNode()

            if feature.strand not in feature2onto:
                raise ValueError("Incorrect strand information for feature ID '%s'." % feature.id)
            strand = feature2onto[feature.strand]

            g.add( (s, p, o) )
            g.add( (s, RDFS.label, Literal(feature.featuretype, datatype=XSD.string)) )
            g.add( (s, RDF.type, FALDO.Region) )

            # add feature start/end positions and strand info
            g.add( (s, FALDO.begin, start) )
            g.add( (start, RDF.type, FALDO.ExactPosition) )
            g.add( (start, RDF.type, strand) )
            g.add( (start, FALDO.position, Literal(feature.start, datatype=XSD.nonNegativeInteger)) )
            g.add( (start, FALDO.reference, Literal(feature.seqid, datatype=XSD.string)) )
            g.add( (s, FALDO.begin, end) )
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

    # loop through GFF files, populate GFF databases and write RDF graphs
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
