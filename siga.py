"""
Generate semantically interoperable genome annotations from GFF[1] files according to the RDF[2] specification.

References:
[1] Generic Feature Format specification, http://www.sequenceontology.org/gff3.shtml
[2] Resource Description Framework, https://www.w3.org/TR/rdf11-concepts/

Usage:
  siga.py -h | --help
  siga.py -v | --version
  siga.py [-V ] [ -d DB_FILE | -e DB_FILE_EXT ] [ -o FORMAT ] GFF_FILE...

Arguments:
  GFF_FILE...      Input file(s) in GFF (versions 2/3).

Options:
  -h, --help
  -v, --version
  -V, --verbose    Use for debugging.
  -d DB_FILE       Populate GFF database(s) in SQLite.
  -e DB_FILE_EXT   Database file extension. [default: .db].
  -o FORMAT        Select RDF serialization/format. [default: turtle]

"""

from __future__ import print_function
from docopt import docopt

import os
import gffutils as gff
import sqlite3 as sql

__author__  = 'Arnold Kuzniar'
__version__ = 0.1
__status__  = 'Prototype'
__license__ = 'Apache License, Version 2.0'


if __name__ == '__main__':
    args = docopt(__doc__, version=__version__)
    #print(args)

    def normalize_filext(s):
        dot = '.'
        if s.startswith(dot) is False:
            s = dot + s
        return s

    #def print_all_features(db):
    #    for ft in db.all_features():
    #        print(ft)

    for gff_file in args['GFF_FILE']:
        if args['-d']:
            # populate a single db from all GFF files
            db_file = args['-d']
            if os.path.exists(db_file):
                db = gff.FeatureDB(db_file)
                db.update(gff_file)
            else:
                db = gff.create_db(gff_file, db_file, force=False)
        else:
            # populate one db per GFF file
            base_name = os.path.splitext(gff_file)[0]
            db_file = base_name + normalize_filext(args['-e'])
            try:
                db = gff.create_db(gff_file, db_file, verbose=args['--verbose'], force=False)
                #print_all_features(db)
            except sql.OperationalError:
                raise IOError("Database file '%s' already exists." % db_file)
            except ValueError:
                raise IOError("GFF file '%s' not found." % gff_file)
