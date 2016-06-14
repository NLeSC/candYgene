"""
Generate semantically interoperable genome annotations from GFF[1] files according to the RDF[2] specification.

References:
[1] Generic Feature Format specification, http://www.sequenceontology.org/gff3.shtml
[2] Resource Description Framework, https://www.w3.org/TR/rdf11-concepts/

Usage:
  siga.py -h | --help
  siga.py -v | --version
  siga.py [-iV ] [ -d DB_FILE | -e DB_FILE_EXT ] [ -o FORMAT ] GFF_FILE...

Arguments:
  GFF_FILE...      Input file(s) in GFF (versions 2/3).

Options:
  -h, --help
  -v, --version
  -V, --verbose    Use for debugging.
  -d DB_FILE       Populate GFF database(s) in SQLite.
  -e DB_FILE_EXT   Database file extension. [default: .db].
  -i               Check the referential integrity of database(s).
  -o FORMAT        Select RDF graph serialization: turtle (.ttl), n3 (.n3) or xml (.rdf). [default: turtle]

"""

from __future__ import print_function
from docopt import docopt

import os
import gffutils as gff
import sqlite3 as sql

from rdflib import Graph, URIRef, Literal
from rdflib.namespace import Namespace, RDF, RDFS, XSD

__author__  = 'Arnold Kuzniar'
__version__ = '0.1.1'
__status__  = 'Prototype'
__license__ = 'Apache License, Version 2.0'


def normalize_filext(s):
    dot = '.'
    if s.startswith(dot) is False:
        s = dot + s
    return s

def triplify(db, format):
    format2fext = dict(turtle = '.ttl', xml = '.rdf', n3 = '.n3')
    assert(format in format2fext), "Unsupported RDF serialization '%s'." % format

    # define additional namespaces
    SO = Namespace('http://purl.obolibrary.org/obo/') # Sequence Ontology (Feature Annotation)
    BASE = Namespace('https://solgenomics.net/search/quick?term=') # so far no better solution to resolve features

    # bind prefixes to namespaces
    g = Graph()
    g.bind('so', SO)
    g.bind(None, BASE)

    # map feature types to SO accessions
    ft2id = dict(gene = 'SO_0000704',
                 mRNA = 'SO_0000234',
                 CDS = 'SO_0000316',
                 exon = 'SO_0000147',
                 intron = 'SO_0000188',
                 five_prime_UTR = 'SO_0000204',
                 three_prime_UTR = 'SO_0000205')

    query = 'SELECT * FROM features'
    for row in db.execute(query):
        ft = row['featuretype']
        id = row['id'].replace(ft + ':', '')
        if ft in ft2id:
            s = URIRef(BASE + id)
            p = URIRef(RDF.type)
            o = URIRef(SO + ft2id[ft])
            g.add( (s, p, o) )
            g.add( (s, RDFS.label, Literal(ft, datatype=XSD.string)))

    outfile = os.path.splitext(db.dbfn)[0] + format2fext[format]
    with open(outfile, 'w') as fout:
        fout.write(g.serialize(format=format))
 
if __name__ == '__main__':
    args = docopt(__doc__, version=__version__)
    #print(args)

    format = args['-o']
    fk_constraints = 'ON' if args['-i'] is True else 'OFF'
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
                db = gff.create_db(gff_file, db_file, pragmas=pragmas, force=False)
        else:
            # populate one db per GFF file
            base_name = os.path.splitext(gff_file)[0]
            db_file = base_name + normalize_filext(args['-e'])
            try:
                db = gff.create_db(gff_file, db_file, verbose=args['--verbose'], pragmas=pragmas, force=False)
                triplify(db, format)
            except sql.OperationalError:
                raise IOError("Database file '%s' already exists." % db_file)
            except ValueError:
                raise IOError("GFF file '%s' not found." % gff_file)
            except sql.IntegrityError, e:
                raise IOError("%s in database '%s'." % (e, db_file))

    if args['-d']:
        db_file = args['-d']
        db = gff.FeatureDB(db_file)
        triplify(db, format)
