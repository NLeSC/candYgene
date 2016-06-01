"""
Make genome annotations in GFF[1] files semantically interoperable according to the RDF[2] specification.

References:
[1] Generic Feature Format specification, http://www.sequenceontology.org/gff3.shtml
[2] Resource Description Framework, https://www.w3.org/TR/rdf11-concepts/

Usage:
  siga.py -h | --help
  siga.py -v | --version
  siga.py [-f FORMAT] [-o PATH] [-e EXT] [-r] FILE...

Arguments:
  FILE...       Input file(s) in GFF.

Options:
  -h, --help
  -v, --version
  -f FORMAT, --format=FORMAT  Select RDF serialization/format for output file(s) [default: turtle].
  -o PATH, --outpath=PATH     Select output directory [default: ./].
  -e EXT, --extension=EXT     GFF database file extension [default: .sqlite]
  -r, --replace               Replace existing GFF database file.
"""

__author__  = 'Arnold Kuzniar'
__version__ = 0.1
__status__  = 'Prototype'
__license__ = 'Apache License, Version 2.0'

import gffutils as gff
import os
from docopt import docopt


if __name__ == '__main__':
    args = docopt(__doc__, version=__version__)
    #print(args)

    for gff_file in args['FILE']:
        base_name = os.path.splitext(gff_file)[0]
        file_ext = args['--extension']
        if not file_ext.startswith('.'):
            file_ext = '.' + file_ext
        db_file = os.path.join(args['--outpath'], base_name + file_ext)
        db = gff.create_db(gff_file, db_file, force=args['--replace'])

