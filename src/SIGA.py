"""
SIGA.py is a command-line tool to generate Semantically Interoperable Genome Annotations from
GFF files [1,2] according to the RDF specification [3].

References:
[1] Generic Feature Format specification, http://www.sequenceontology.org/
[2] DDBJ/ENA/GenBank Feature Table Definition, http://www.insdc.org/
[3] Resource Description Framework, https://www.w3.org/TR/rdf11-concepts/

Usage:
  SIGA.py -h|--help
  SIGA.py -v|--version
  SIGA.py db [-ruV] [-d DB_FILE | -e DB_FILEXT] GFF_FILE...
  SIGA.py rdf [-V] [-o FORMAT] (-C CFG_FILE | -b URI -c URI -s URL -n NAME -t ID) DB_FILE...

Arguments:
  GFF_FILE...      Input file(s) in GFF version 2 or 3.
  DB_FILE...       Input database file(s) in SQLite.

Options:
  -h, --help
  -v, --version
  -V, --verbose    Show verbose output (debug mode).
  -b URI           Set the base URI of the generated RDF graph (e.g. https://solgenomics.net/).
  -c URI           Set the URI of a person, organization or service making the resource
                   available in RDF (e.g. http://orcid.org/0000-0003-1711-7961).
  -C FILE          Set the config file path.
  -s URL           Set the source or download URL of GFF file(s).
  -n NAME          Set the genome or species name (e.g. Solanum lycopersicum).
  -t ID            Set the NCBI Taxonomy ID of the species (e.g. 4081).
  -d DB_FILE       Create features database from GFF file(s).
  -e DB_FILEXT     Set the database file extension [default: .db].
  -r               Check the referential integrity of the database(s).
  -u               Generate unique feature IDs if duplicates found.
  -o FORMAT        Set the format of the generated RDF graph:
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
#   ../genome/<species name>/<feature type>/<feature ID> + [#<begin|end>|#<start>-<end>] only for chromosome
#
# FIXME: Sequence Ontology (SO) properties:
#   SO_genome_of
#   SO_part_of
#   SO_has_part
#   SO_transcribed_to
#   SO_translated_to
# do not resolve via http://purl.obolibrary.org/obo/.
#

from __future__ import print_function
from docopt import docopt
from rdflib import Graph, URIRef, Literal
from rdflib.namespace import Namespace, RDF, RDFS, XSD, DCTERMS
from urllib2 import urlparse, unquote
from datetime import datetime
from ConfigParser import SafeConfigParser

import os
import re
import gffutils as gff
import sqlite3 as sql

__author__  = 'Arnold Kuzniar'
__version__ = '0.3.9'
__status__  = 'alpha'
__license__ = 'Apache License, Version 2.0'


def init_config():
    """Initialize config dictionary with mandatory sections and attributes (keys)."""
    config = dict(URIs = dict(rdf_base = None,
                              rdf_creator = None,
                              gff_source = None),
                  Dataset = dict(species_name = None,
                                 ncbi_taxon_id = None))
    return config


def read_config_file(fname):
    """Read config file and return dataset-related setting."""
    parser = SafeConfigParser()
    parser.read(fname)

    if os.path.isfile(fname) is False:
        raise IOError("Config file '{0}' not found.".format(fname))

    for section in parser.sections():
        if section not in config.keys():
            raise ValueError("Unsupported section '{0}' in the config file '{1}.".format(section, fname))
        for attr,val in parser.items(section):
            if attr not in config[section].keys():
                raise ValueError("Unsupported attribute '{0}' in config file '{1}'.".format(attr, section, fname))
            config[section][attr] = val
    return config


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

    # Note: There is no optimal solution to resolve all features in SGN (https://solgenomics.net/):
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


def get_feature_attrs(ft):
    """Get feature attributes and concatenate the the values into a single string."""
    attrs = {} # selected attributes
    des = []

    for attr in ['Name', 'Note', 'Alias', 'Ontology_term', 'Interpro2go_term', 'Sifter_term']:
        attrs[attr] = None
        attrs[attr.lower()] = None

    for key in ft.attributes.keys():
        if key in attrs:
            val = ft[key]
            val = ', '.join(val) if type(val) == list else str(val)
            des.append('{0}: {1}'.format(key, val.encode('utf-8')))

    if len(des) == 0:
        return None
    else:
        return unquote('; '.join(sorted(des)))


def amend_feature_type(ft):
    """Change the original feature type to correct or more specific type
       according to the DDBJ/ENA/GenBank Feature Table Definition.
    """
    # TODO: add mappings to a config file
    # 'match' or 'match_part' features are sometimes (mis)used to indicate polymorphic sites
    # instead of using the 'variation' key according to [2].
    feature_types = dict(
        mRNA = 'prim_transcript',
        match = 'variation',
        match_part = 'variation')
    if ft in feature_types:
        return feature_types[ft]
    return ft


def triplify(db, rdf_format, config):
    """Generate RDF triples from RDB using Direct Mapping approach."""
    fmt2fext = dict(xml = '.rdf',
                    nt = '.nt',
                    turtle = '.ttl',
                    n3 = '.n3')

    if rdf_format not in fmt2fext:
        raise IOError("Unsupported RDF serialization '{0}'.".format(rdf_format))

    base_uri = config['URIs']['rdf_base']
    creator_uri = config['URIs']['rdf_creator']
    download_url = config['URIs']['gff_source']
    species_name = config['Dataset']['species_name']
    taxon_id = config['Dataset']['ncbi_taxon_id']

    # define additional namespace prefixes
    # TODO: add namespaces to a config file
    OBO = Namespace('http://purl.obolibrary.org/obo/')
    FALDO = Namespace('http://biohackathon.org/resource/faldo#')
    DCMITYPE = Namespace('http://purl.org/dc/dcmitype/')

    g = Graph()
    g.bind('obo', OBO)
    g.bind('faldo', FALDO)
    g.bind('dcterms', DCTERMS)
    g.bind('dcmitype', DCMITYPE)

    # map GFF feature types and DNA strandedness to ontology classes
    # Note: The 'mRNA' feature key is often used (incorrectly) in place of 'prim_transcript'
    # in genome annotations. The former feature MUST NOT contain introns while the latter
    # MAY contain introns [2].
    # Feature type to SO mappings:
    #   prim_transcript -> SO_0000120 refers to a protein-coding primary (unprocessed) transcript
    #   mRNA            -> SO_0000234 refers to a mature transcript
    #
    feature_onto_class = {
        'genome'          : OBO.SO_0001026,
        'chromosome'      : OBO.SO_0000340,
        'gene'            : OBO.SO_0000704,
        'prim_transcript' : OBO.SO_0000120,
        'mRNA'            : OBO.SO_0000234,
        'CDS'             : OBO.SO_0000316,
        'exon'            : OBO.SO_0000147,
        'intron'          : OBO.SO_0000188,
        'five_prime_UTR'  : OBO.SO_0000204,
        'three_prime_UTR' : OBO.SO_0000205,
        'polyA_site'      : OBO.SO_0000553,
        'polyA_sequence'  : OBO.SO_0000610,
        'variation'       : OBO.SO_0001645
    }

    strand_onto_class = {
        '+' : FALDO.ForwardStrandPosition,
        '-' : FALDO.ReverseStrandPosition,
        '?' : FALDO.StrandedPosition,
        '.' : FALDO.Position
    }

    # add genome info to graph
    genome_uri = URIRef(os.path.join(base_uri, 'genome', species_name.replace(' ', '_')))
    taxon_uri = OBO.term('NCBITaxon_%d' % taxon_id)
    g.add( (genome_uri, RDF.type, feature_onto_class['genome']) )
    g.add( (genome_uri, RDF.type, DCMITYPE.Dataset) )
    g.add( (genome_uri, RDFS.label, Literal('genome of {0}'.format(species_name), datatype=XSD.string)) )
    g.add( (genome_uri, DCTERMS.created, Literal(datetime.now().strftime("%Y-%m-%d"), datatype=XSD.date )) )
    g.add( (genome_uri, DCTERMS.creator, URIRef(creator_uri)) )
    g.add( (genome_uri, DCTERMS.title, Literal('genome of {0}'.format(species_name), datatype=XSD.string)) )
    g.add( (genome_uri, DCTERMS.source, URIRef(download_url)) )
    g.add( (genome_uri, OBO.SO_genome_of, taxon_uri) ) # N.B.: predicate has no domain/range defined
    g.add( (genome_uri, OBO.RO_0002162, taxon_uri) )   # use 'in taxon' alternatively
    g.add( (taxon_uri, RDFS.label, Literal('NCBI Taxonomy ID: {0}'.format(taxon_id), datatype=XSD.string)) )
    g.add( (taxon_uri, DCTERMS.identifier, Literal(taxon_id, datatype=XSD.positiveInteger)) )

    for feature in db.all_features():
        if feature.strand not in strand_onto_class:
            raise KeyError("Incorrect strand information for feature ID '{0}'.".format(feature.id))
        try: # skip GFF feature types not in feature_onto_class dict
            strand_uri = strand_onto_class[feature.strand]
            feature_id = normalize_feature_id(feature.id)
            feature_type = amend_feature_type(feature.featuretype)
            feature_type_uri = URIRef(feature_onto_class[feature_type])
            feature_uri = URIRef(os.path.join(genome_uri, feature_type, feature_id))
            seqid_uri = URIRef(os.path.join(genome_uri, 'chromosome', str(feature.seqid)))
            region_uri = URIRef('{0}#{1}-{2}'.format(seqid_uri, feature.start, feature.end))
            start_uri = URIRef('{0}#{1}'.format(seqid_uri, feature.start))
            end_uri = URIRef('{0}#{1}'.format(seqid_uri, feature.end))
            label = '{0} {1}'.format(feature_type, feature_id)

            # add genome and chromosome info to graph
            # Note: the assumption is that the seqid field refers to chromosome
            g.add( (seqid_uri, RDF.type, feature_onto_class['chromosome']) )
            g.add( (seqid_uri, RDFS.label, Literal('chromosome {0}'.format(feature.seqid), datatype=XSD.string)) )
            g.add( (seqid_uri, OBO.SO_part_of, genome_uri) )

            # add feature types and IDs to graph
            g.add( (feature_uri, RDF.type, feature_type_uri) )
            g.add( (feature_uri, RDFS.label, Literal(label, datatype=XSD.string)) )
            g.add( (feature_uri, DCTERMS.identifier, Literal(feature_id, datatype=XSD.string)) )

            # add feature descriptions (from the attributes field) to graph
            des = get_feature_attrs(feature)
            if des is not None:
                g.add( (feature_uri, RDFS.comment, Literal(des, datatype=XSD.string)) )

            # add feature start/end coordinates and strand info to graph
            g.add( (feature_uri, FALDO.location, region_uri) )
            g.add( (region_uri, RDF.type, FALDO.Region) )
            g.add( (region_uri, RDFS.label, Literal('region {0}-{1} on chromosome {2}'.format(feature.start,
                                                                                              feature.end,
                                                                                              feature.seqid))) )
            g.add( (region_uri, FALDO.begin, start_uri) )
            g.add( (start_uri, RDF.type, FALDO.ExactPosition) )
            g.add( (start_uri, RDF.type, strand_uri) )
            g.add( (start_uri, RDFS.label, Literal('position at {0} on chromosome {1}'.format(feature.start,
                                                                                              feature.seqid))) )
            g.add( (start_uri, FALDO.position, Literal(feature.start, datatype=XSD.positiveInteger)) )
            g.add( (start_uri, FALDO.reference, seqid_uri) )
            g.add( (region_uri, FALDO.end, end_uri) )
            g.add( (end_uri, RDF.type, FALDO.ExactPosition) )
            g.add( (end_uri, RDF.type, strand_uri) )
            g.add( (end_uri, RDFS.label, Literal('position at {0} on chromosome {1}'.format(feature.end,
                                                                                            feature.seqid))) )
            g.add( (end_uri, FALDO.position, Literal(feature.end, datatype=XSD.positiveInteger)) )
            g.add( (end_uri, FALDO.reference, seqid_uri) )
            # TODO: phase info is mandatory for CDS feature types but can't find a corresponding ontology term

            # add parent-child relationships between features to graph
            for child in db.children(feature, level=1):
                child_feature_id = normalize_feature_id(child.id)
                child_feature_type = amend_feature_type(child.featuretype)
                child_feature_uri = URIRef(os.path.join(genome_uri, child_feature_type, child_feature_id))
                g.add( (feature_uri, OBO.SO_has_part, child_feature_uri) ) # use the inverse of part_of

                if feature_type == 'gene' and child_feature_type == 'prim_transcript':
                    g.add( (feature_uri, OBO.SO_transcribed_to, child_feature_uri) )

        except KeyError:
            pass

    outfile = os.path.splitext(db.dbfn)[0] + fmt2fext[rdf_format]
    with open(outfile, 'w') as fout:
        fout.write(g.serialize(format=rdf_format))


if __name__ == '__main__':
    args = docopt(__doc__, version=__version__)
    debug = args['--verbose']

    if args['db'] is True: # in db mode
        unique_keys = 'create_unique' if args['-u'] is True else 'error'
        fk_check = 'ON' if args['-r'] is True else 'OFF'
        pragmas = dict(foreign_keys=fk_check)

        # populate database(s) from GFF file(s)
        for gff_file in args['GFF_FILE']:
            if os.path.exists(gff_file) is False:
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
        rdf_format = args['-o']
        cfg_file = args['-C']
        config = init_config()

        if cfg_file is None:
	    config['URIs']['rdf_base'] = validate_uri(args['-b'])
            config['URIs']['rdf_creator'] = validate_uri(args['-c'])
            config['URIs']['gff_source'] = validate_uri(args['-s'])
            config['Dataset']['species_name'] = args['-n']
            config['Dataset']['ncbi_taxon_id'] = args['-t']
        else:
            config = read_config_file(cfg_file)

        try:
            taxon_id = config['Dataset']['ncbi_taxon_id']
            config['Dataset']['ncbi_taxon_id'] = int(taxon_id)
        except:
            raise ValueError('Enter valid NCBI Taxonomy ID.')

        # serialize RDF graphs from db files
        for db_file in args['DB_FILE']:
            db = gff.FeatureDB(db_file)
            triplify(db, rdf_format, config)

