"""
SIGA.py is a command-line tool to generate Semantically Interoperable Genome Annotations from
GFF files [1] according to the RDF specification [2].

References:
[1] Generic Feature Format specification, http://www.sequenceontology.org/
[2] Resource Description Framework, https://www.w3.org/TR/rdf11-concepts/

Usage:
  SIGA.py -h|--help
  SIGA.py -v|--version
  SIGA.py db [-cuV] [-d DB_FILE|-e DB_FILEXT] GFF_FILE...
  SIGA.py rdf [-V] [-o FORMAT] -b BASE_URI -D URL -g NAME -t ID DB_FILE...

Arguments:
  GFF_FILE...      Input file(s) in GFF version 2 or 3.
  DB_FILE...       Input database file(s) in SQLite.

Options:
  -h, --help
  -v, --version
  -V, --verbose    Use for debugging.
  -b BASE_URI      Set the base URI (e.g. https://solgenomics.net/).
  -D URL           Set the download URL of GFF file(s).
  -g NAME          Set the genome or species name (e.g. Solanum lycopersicum).
  -t ID            Set the NCBI Taxonomy ID of the species (e.g. 4081).
  -d DB_FILE       Create a database from GFF file(s).
  -e DB_FILEXT     Set the database file extension [default: .db].
  -c               Check the referential integrity of the database(s).
  -u               Generate unique feature IDs if duplicates found.
  -o FORMAT        Select RDF output format:
                     turtle (.ttl) [default: turtle]
                     xml (.rdf),
                     nt (.nt),
                     n3 (.n3)
"""

#
# Supported feature types:
#   genome, chromosome, gene, prim_transcript, mRNA, intron, exon, CDS,
#   three_prime_UTR, five_prime_UTR, polyA_site, polyA_sequence
#
# Defined URI data space relative to base URI:
#   ../<feature type>/<feature ID> + [#<begin|end>|#<start>-<end>] only for chromosome
#

from __future__ import print_function
from docopt import docopt
from rdflib import Graph, URIRef, Literal, BNode
from rdflib.namespace import Namespace, RDF, RDFS, XSD, DCTERMS
from urllib2 import urlparse, unquote
from datetime import datetime

import os
import re
import gffutils as gff
import sqlite3 as sql

__author__  = 'Arnold Kuzniar'
__version__ = '0.3.5'
__status__  = 'alpha'
__license__ = 'Apache License, Version 2.0'


def normalize_filext(s):
    """Prefix file extension with a dot '.' if not done already."""
    dot = '.'
    if s.startswith(dot) is False:
        s = dot + s
    return s


def remove_file(fn):
    """Remove file if exists."""
    if os.path.exists(fn) is True:
        os.remove(fn)


def validate_uri(uri):
    """Validate input URI (scheme and host)."""
    u = urlparse.urlparse(uri)
    if u.scheme not in ('http', 'https', 'ftp'):
        raise ValueError("Invalid URI scheme used in '%s'." % uri)
    if u.netloc is '':
        raise ValueError('No host specified.')
    return u.geturl()


def normalize_feature_id(id):
    """Ad-hoc function to normalize feature IDs."""

    # Note: There is no optimal/universal solution to resolve all features in SGN (https://solgenomics.net/):
    # e.g., a feature ID in a GFF file 'gene:Solyc00g005000.2' corresponds to three URLs:
    #   https://solgenomics.net/locus/Solyc00g005000.2/view !!! note the use of 'locus' instead of 'gene' !!!
    #   https://solgenomics.net/locus/Solyc00g005000/view
    #   https://solgenomics.net/feature/17660839/details !!! note the use of internal IDs !!!
    #
    # The search term 'Solyc00g005000' returns a page with links to:
    #   tomato locus  https://solgenomics.net/locus/8377/view !!! where the term is referred to as locus (name|symbol)
    #   gene feature  https://solgenomics.net/feature/17660839/details
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
    # "Normalizing" feature IDs by removing the prefixes seems reasonable for most feature types, except
    # for the UTRs, which would have ambiguous feature IDs, e.g., Solyc00g005000.2.1.0 for both
    # 'five_prime_UTR:Solyc00g005000.2.1.0' and 'three_prime_UTR:Solyc00g005000.2.1.0'
    #
    return re.sub('gene:|mRNA:|CDS:|exon:|intron:|\w+UTR:', '', id)


def get_feature_attr(feature, attr):
    """Get feature attribute."""
    try:
        return feature[attr]
    except KeyError:
        try:
            return feature[attr.lower()]
        except KeyError:
            return None


def triplify(db, rdf_format, base_uri, download_url, species_name, taxon_id):
    """Generate RDF triples from RDB using Direct Mapping approach."""
    fmt2fext = dict(xml = '.rdf',
                    nt = '.nt',
                    turtle = '.ttl',
                    n3 = '.n3')

    if rdf_format not in fmt2fext:
        raise IOError("Unsupported RDF serialization '{0}'.".format(rdf_format))

    # define additional namespace prefixes
    SO = Namespace('http://purl.obolibrary.org/obo/so.owl#') # FIXME: URI resolution of classes/properties
    FALDO = Namespace('http://biohackathon.org/resource/faldo#')
    TAXON = Namespace('http://purl.obolibrary.org/obo/ncbitaxon.owl#')
    g = Graph()
    g.bind('so', SO)
    g.bind('faldo', FALDO)
    g.bind('taxon', TAXON)
    g.bind('dct', DCTERMS)

    # map GFF feature types and DNA strandedness to ontology classes
    # Note: The 'mRNA' feature key is often used (incorrectly) in place of 'prim_transcript'
    # in genome annotations. The former feature MUST NOT contain introns while the latter
    # MAY contain introns (see DDBJ/ENA/GenBank Feature Table Definition, http://www.insdc.org/)
    # FT to SO mappings:
    #   prim_transcript -> SO:0000120 refers to a protein-coding primary (unprocessed) transcript
    #   mRNA            -> SO_0000234 refers to a mature transcript
    #
    feature_onto_class = {
        'genome'          : SO.SO_0001026,
        'chromosome'      : SO.SO_0000340,
        'gene'            : SO.SO_0000704,
        'prim_transcript' : SO.SO_0000120,
        'mRNA'            : SO.SO_0000120, # N.B.
        'CDS'             : SO.SO_0000316,
        'exon'            : SO.SO_0000147,
        'intron'          : SO.SO_0000188,
        'five_prime_UTR'  : SO.SO_0000204,
        'three_prime_UTR' : SO.SO_0000205,
        'polyA_site'      : SO.SO_0000553,
        'polyA_sequence'  : SO.SO_0000610,
        '+' : FALDO.ForwardStrandPosition,
        '-' : FALDO.ReverseStrandPosition,
        '?' : FALDO.StrandedPosition,
        '.' : FALDO.Position
    }
    gff.constants.always_return_list = False # return GFF attributes as string

    # add genome info to graph
    genome = URIRef(os.path.join(base_uri, 'genome', species_name.replace(' ', '_')))
    taxon = TAXON.term(str(taxon_id))

    g.add( (genome, RDF.type, feature_onto_class['genome']) )
    g.add( (genome, DCTERMS.created, Literal(datetime.now().strftime("%Y-%m-%d"), datatype=XSD.date )) )
    for pred in (RDFS.label, DCTERMS.title):
        g.add( (genome, pred, Literal('genome of {0}'.format(species_name), datatype=XSD.string)) )
    g.add( (genome, DCTERMS.source, URIRef(download_url)) )
    g.add( (genome, SO.genome_of, taxon) )
    g.add( (taxon, RDFS.label, Literal('NCBI Taxonomy ID: {0}'.format(taxon_id), datatype=XSD.string)) )
    g.add( (taxon, RDF.value, Literal(taxon_id, datatype=XSD.positiveInteger)) )

    for feature in db.all_features():
        if feature.strand not in feature_onto_class:
            raise KeyError("Incorrect strand information for feature ID '{0}'.".format(feature.id))
        try: # skip GFF feature types not in feature_onto_class dict
            strand = feature_onto_class[feature.strand]
            feature_id = normalize_feature_id(feature.id)
            feature_type = URIRef(feature_onto_class[feature.featuretype])
            feature_parent = URIRef(os.path.join(base_uri, feature.featuretype, feature_id))
            seqid = URIRef(os.path.join(base_uri, 'chromosome', str(feature.seqid)))
            region = URIRef('{0}#{1}-{2}'.format(seqid, feature.start, feature.end))
            start = URIRef('{0}#{1}'.format(seqid, feature.start))
            end = URIRef('{0}#{1}'.format(seqid, feature.end))
            label = '{0} {1}'.format(feature.featuretype.replace('_', ' '), feature_id)

            # add genome and chromosome info to graph
            # Note: the assumption is that the `seqid` field refers to chromosome
            g.add( (seqid, RDF.type, feature_onto_class['chromosome']) )
            g.add( (seqid, RDFS.label, Literal('chromosome {0}'.format(feature.seqid), datatype=XSD.string)) )
            g.add( (seqid, SO.part_of, genome) )

            # add feature types to graph
            g.add( (feature_parent, RDF.type, feature_type) )
            g.add( (feature_parent, RDFS.label, Literal(label, datatype=XSD.string)) )

            # add feature description to graph
            arr = []
            for attr in ('Note', 'Name'):
                el = get_feature_attr(feature, attr)
                if el is not None:
                    arr.append(str(el))

            if len(arr) != 0:
                desc = unquote(' '.join(arr))
                if desc != feature_id:
                    g.add( (feature_parent, RDFS.comment, Literal(desc, datatype=XSD.string)) )

            # add feature start/end coordinates and strand info to graph
            g.add( (feature_parent, FALDO.location, region) )
            g.add( (region, RDF.type, FALDO.Region) )
            g.add( (region, RDFS.label, Literal('region {0}-{1} on chromosome {2}'.format(feature.start,
                                                                                          feature.end,
                                                                                          feature.seqid))) )
            g.add( (region, FALDO.begin, start) )
            g.add( (start, RDF.type, FALDO.ExactPosition) )
            g.add( (start, RDF.type, strand) )
            g.add( (start, RDFS.label, Literal('position at {0} on chromosome {1}'.format(feature.start, feature.seqid))) )
            g.add( (start, FALDO.position, Literal(feature.start, datatype=XSD.positiveInteger)) )
            g.add( (start, FALDO.reference, seqid) )
            g.add( (region, FALDO.end, end) )
            g.add( (end, RDF.type, FALDO.ExactPosition) )
            g.add( (end, RDF.type, strand) )
            g.add( (end, RDFS.label, Literal('position at {0} on chromosome {1}'.format(feature.end, feature.seqid))) )
            g.add( (end, FALDO.position, Literal(feature.end, datatype=XSD.positiveInteger)) )
            g.add( (end, FALDO.reference, seqid) ) 
            # TODO: phase info is mandatory for CDS feature types but can't find a corresponding ontology term

            # add parent-child relationships between features to graph
            for child in db.children(feature, level=1):
                feature_id = normalize_feature_id(child.id)
                feature_child = URIRef(os.path.join(base_uri, child.featuretype, feature_id))
                g.add( (feature_parent, SO.has_part, feature_child) ) # use the inverse of part_of

                if feature.featuretype == 'gene' and child.featuretype in ('mRNA', 'prim_transcript'):
                    g.add( (feature_parent, SO.transcribed_to, feature_child) )

        except KeyError:
            pass

    outfile = os.path.splitext(db.dbfn)[0] + fmt2fext[rdf_format]
    with open(outfile, 'w') as fout:
        fout.write(g.serialize(format=rdf_format))


if __name__ == '__main__':
    args = docopt(__doc__, version=__version__)
    #print(args)
    debug = args['--verbose']

    if args['db'] is True: # in db mode
        unique_keys = 'create_unique' if args['-u'] is True else 'error'
        fk_check = 'ON' if args['-c'] is True else 'OFF'
        pragmas = dict(foreign_keys=fk_check)

        # populate database(s) from GFF file(s)
        for gff_file in args['GFF_FILE']:
            if os.path.exists(gff_file) is False:
               remove_file(db_file)
               raise IOError("GFF file '{0}' not found.".format(gff_file))

            if args['-d']: # one db for all GFF files
                db_file = args['-d']
                if os.path.exists(db_file) is True:
                    db = gff.FeatureDB(db_file)
                    db.update(gff_file)
                else:
                    db = gff.create_db(gff_file, db_file, merge_strategy=unique_keys, verbose=debug, pragmas=pragmas, force=False)
            else: # one db per GFF file
                base_name = os.path.splitext(gff_file)[0]
                db_file = base_name + normalize_filext(args['-e'])
                try:
                    db = gff.create_db(gff_file, db_file, merge_strategy=unique_keys, verbose=debug, pragmas=pragmas, force=False)
                except sql.OperationalError:
                    raise IOError("Database file '{0}' already exists.".format(db_file))
                except sql.IntegrityError, err:
                    remove_file(db_file)
                    raise IOError("{0} in database '{1}'.".format(err, db_file))
    else: # in rdf mode
        base_uri = validate_uri(args['-b'])
        download_url = validate_uri(args['-D'])
        rdf_format = args['-o']
        species_name = args['-g']
        taxon_id = args['-t']
        try:
            taxon_id = int(taxon_id)
        except:
            raise ValueError('Enter valid NCBI Taxonomy ID.')

        # serialize RDF graphs from db files
        for db_file in args['DB_FILE']:
            db = gff.FeatureDB(db_file)
            triplify(db, rdf_format, base_uri, download_url, species_name, taxon_id)

