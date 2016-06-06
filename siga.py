"""
Generate semantically interoperable genome annotations from GFF[1] files according to the RDF[2] specification.

References:
[1] Generic Feature Format specification, http://www.sequenceontology.org/gff3.shtml
[2] Resource Description Framework, https://www.w3.org/TR/rdf11-concepts/

Usage:
  siga.py -h | --help
  siga.py -v | --version
  siga.py [ -d DB_FILE | -e DB_FILE_EXT ] [ -o FORMAT ] GFF_FILE...

Arguments:
  GFF_FILE...      Input file(s) in GFF (versions 2/3).

Options:
  -h, --help
  -v, --version
  -d DB_FILE       Populate GFF database(s) in SQLite.
  -e DB_FILE_EXT   Database file extension [default: .db].
  -o FORMAT        Select RDF serialization/format [default: turtle].

"""

from __future__ import print_function
from docopt import docopt

import os
import gffutils as gff

__author__  = 'Arnold Kuzniar'
__version__ = 0.1
__status__  = 'Prototype'
__license__ = 'Apache License, Version 2.0'


if __name__ == '__main__':
    args = docopt(__doc__, version=__version__)
    print(args)

    db_file = args['-d']
    for gff_file in args['GFF_FILE']:
        # generate one db per GFF file
        if args['-d'] is None:
            base_name = os.path.splitext(gff_file)[0]
            file_ext = args['-e']
            char = '.'
            if file_ext.startswith(char) is False:
                file_ext = char + file_ext
            db_file = base_name + file_ext
            db = gff.create_db(gff_file, db_file, force=True)

        # import all GFF files into a single db
        elif os.path.exists(db_file) is True:
            db = gff.FeatureDB(db_file)
            db.update(gff_file)
        else:
            db = gff.create_db(gff_file, db_file, force=False)

