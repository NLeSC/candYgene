"""
SIGA.py is a command-line tool to generate Semantically Interoperable Genome Annotations from
GFF files [1] according to the RDF specification [2].

References:
[1] Generic Feature Format specification, http://www.sequenceontology.org/gff3.shtml
[2] Resource Description Framework, https://www.w3.org/TR/rdf11-concepts/

Usage:
  SIGA.py -h | --help
  SIGA.py -v | --version
  SIGA.py [-cV ] [ -d DB_FILE | -e DB_FILE_EXT ] [ -o FORMAT ] -b BASE_URI GFF_FILE...

Arguments:
  GFF_FILE...      Input file(s) in GFF version 2 or 3.

Options:
  -h, --help
  -v, --version
  -V, --verbose    Use for debugging.
  -b BASE_URI      Set base URI (e.g., https://solgenomics.net/).
  -d DB_FILE       Populate a single SQLite database from one or more GFF files.
  -e DB_FILE_EXT   Database file extension [default: .db].
  -o FORMAT        Select RDF serialization format: xml (.rdf), ntriples (.nt), n3 (.n3) or turtle (.ttl) [default: turtle].
  -c               Check the referential integrity of database(s).

"""

from __future__ import print_function
from docopt import docopt
from rdflib import Graph, URIRef, Literal, BNode
from rdflib.namespace import Namespace, RDF, RDFS, XSD
from urllib2 import urlparse, unquote

import os
import re
import gffutils as gff
import sqlite3 as sql

__author__  = 'Arnold Kuzniar'
__version__ = '0.1.7'
__status__  = 'Prototype'
__license__ = 'Apache License, Version 2.0'


def normalize_filext(s):
    dot = '.'
    if s.startswith(dot) is False:
        s = dot + s
    return s


def validate_uri(uri):
    u = urlparse.urlparse(uri)
    if u.scheme not in ('http', 'https', 'ftp'):
        raise ValueError("Invalid URI scheme used in '%s'." % uri)
    if u.netloc is '':
        raise ValueError('No host specified.')
    return u.geturl()


def get_resolvable_uri(uri):
    # Note: There is no optimal/universal solution to resolve all features in SGN (https://solgenomics.net/):
    # e.g., a feature ID in a GFF file 'gene:Solyc00g005000.2' corresponds to three URLs:
    #   https://solgenomics.net/locus/Solyc00g005000.2/view !!! note the use of 'locus' instead of 'gene' !!!
    #   https://solgenomics.net/locus/Solyc00g005000/view
    #   https://solgenomics.net/feature/17660839/details !!! note the use of internal IDs !!!
    #
    # The search term 'Solyc00g005000' returns a page with links to:
    # tomato locus   https://solgenomics.net/locus/8377/view !!! where the term is referred to as locus (name|symbol)
    # gene feature   https://solgenomics.net/feature/17660839/details
    #
    # However, related features such as mRNA, exon, intron do not resolve same way:
    #   mRNA:Solyc00g005000.2.1       https://solgenomics.net/feature/17660840/details
    #   exon:Solyc00g005000.2.1.1     https://solgenomics.net/feature/17660841/details
    #   exon:Solyc00g005000.2.1.2     https://solgenomics.net/feature/17660843/details
    #   intron:Solyc00g005000.2.1.1   https://solgenomics.net/feature/17660842/details
    #   five_prime_UTR and three_prime_UTR do not seem to have corresponding URLs.
    #
    # In principle, feature IDs in the 'attributes' column of a GFF file should be opaque.
    # Currently, the IDs are prefixed with feature type, e.g. 'gene:Solyc00g005000.2'.
    #
    # "Normalizing" feature ID by removing the prefixes seems a reasonable option for most feature types, except
    # for the UTRs, which would have ambiguous feature IDs, e.g., Solyc00g005000.2.1.0 for both
    # 'five_prime_UTR:Solyc00g005000.2.1.0' and 'three_prime_UTR:Solyc00g005000.2.1.0'
    #
    return re.sub('gene:|mRNA:|CDS:|exon:|intron:|\w+UTR:', '', validate_uri(uri)) # this is ugly


def triplify(db, format, base_uri):
    format2filext = dict(xml = '.rdf',
                         ntriples = '.nt',
                         turtle = '.ttl',
                         n3 = '.n3')

    if format not in format2filext:
        raise IOError("Unsupported RDF serialization '%s'." % format)

    # define additional namespace prefixes
    SO = Namespace('http://purl.obolibrary.org/obo/so.owl#')     # Sequence Ontology Feature Annotation
    FALDO = Namespace('http://biohackathon.org/resource/faldo#') # Feature Annotation Location Description Ontology
    #EDAM = Namespace('http://edamontology.org/')                 # EMBRACE Data And Methods ontology

    g = Graph()
    g.bind('so', SO)
    g.bind('faldo', FALDO)
    #g.bind('edam', EDAM)

    # map feature types and DNA strandedness to ontology classes
    feature_onto_class = {
        'gene' : SO.SO_0000704,
        'mRNA' : SO.SO_0000234,
        'CDS'  : SO.SO_0000316,
        'exon' : SO.SO_0000147,
        'intron' : SO.SO_0000188,
        'five_prime_UTR' : SO.SO_0000204,
        'three_prime_UTR' : SO.SO_0000205,
        'chromosome' : SO.SO_0000340,
        '+' : FALDO.ForwardStrandPosition,
        '-' : FALDO.ReverseStrandPosition,
        '?' : FALDO.StrandedPosition,
        '.' : FALDO.Position
    }

    for feature in db.all_features():
        if feature.strand not in feature_onto_class:
            raise KeyError("Incorrect strand information for feature ID '%s'." % feature.id)
        try:
            strand = feature_onto_class[feature.strand]
            feature_type = URIRef(feature_onto_class[feature.featuretype])
            feature_parent = URIRef(get_resolvable_uri(os.path.join(base_uri, feature.featuretype, feature.id)))
            start = BNode()
            end = BNode()

            # add feature types to graph
            g.add( (feature_parent, RDF.type, feature_type) )
            gff.constants.always_return_list = False # return GFF attributes (9th column) as string
            label = feature.attributes.get('Name')
            comment = feature.attributes.get('Note')
            if label is not None: # map Name attr to rdfs:label
                label = "{0} {1}".format(feature.featuretype, label) # concat feature type and name
                g.add( (feature_parent, RDFS.label, Literal(label, datatype=XSD.string)) )

            if comment is not None: # map Note attr to rdfs:comment
                g.add( (feature_parent, RDFS.comment, Literal(unquote(comment), datatype=XSD.string)) )

            # add feature source to graph
            #g.add( (feature_parent, RDF.type, EDAM.data_3034) ) # sequence feature identifier
            #g.add( (feature_parent, EDAM.is_output_of, EDAM.operation_2454) ) # Gene prediction

            # add chromosome, feature start/end coordinates and strand info to graph
            chrom = URIRef(os.path.join(base_uri, 'chromosome', feature.seqid))
            g.add( (chrom, RDF.type, feature_onto_class['chromosome']) )
            g.add( (chrom, RDFS.label, Literal('chromosome %s' % feature.seqid, datatype=XSD.string)) )
            g.add( (feature_parent, RDF.type, FALDO.Region) )
            g.add( (feature_parent, FALDO.begin, start) )
            g.add( (start, RDF.type, FALDO.ExactPosition) )
            g.add( (start, RDF.type, strand) )
            g.add( (start, FALDO.position, Literal(feature.start, datatype=XSD.nonNegativeInteger)) )
            g.add( (start, FALDO.reference, Literal(feature.seqid, datatype=XSD.string)) )
            g.add( (feature_parent, FALDO.end, end) )
            g.add( (end, RDF.type, FALDO.ExactPosition) )
            g.add( (end, RDF.type, strand) )
            g.add( (end, FALDO.position, Literal(feature.end, datatype=XSD.nonNegativeInteger)) )
            g.add( (end, FALDO.reference, chrom) )

            # add parent-child relationships between features to graph
            for child in db.children(feature, level=1):
                feature_child = URIRef(get_resolvable_uri(os.path.join(base_uri, child.featuretype, child.id)))
                g.add( (feature_child, SO.part_of, feature_parent) )
        except KeyError:
            pass

    outfile = os.path.splitext(db.dbfn)[0] + format2filext[format]
    with open(outfile, 'w') as fout:
        fout.write(g.serialize(format=format))
 
if __name__ == '__main__':
    args = docopt(__doc__, version=__version__)
    #print(args)
    base_uri = validate_uri(args['-b'])
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
            triplify(db, format, base_uri)

    if args['-d']:
        db_file = args['-d']
        db = gff.FeatureDB(db_file)
        triplify(db, format, base_uri)
