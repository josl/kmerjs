#!/usr/bin/env python
import pickle
import subprocess
from collections import defaultdict
import io
import json
import subprocess
import sys

# python kmerPyToMongo.py /Users/cisneror/code/genomic-git/kmerfinder/kmer-env/kmerdb-bitbucket/database/complete_genomes.ATGAC.p \
# /Users/cisneror/code/genomic-git/kmerfinder/kmer-env/kmerdb-bitbucket/database/complete_genomes.ATGAC.len.p \
# /Users/cisneror/code/genomic-git/kmerfinder/kmer-env/kmerdb-bitbucket/database/complete_genomes.ATGAC.ulen.p \
# /Users/cisneror/code/genomic-git/kmerfinder/kmer-env/kmerdb-bitbucket/database/complete_genomes.ATGAC.desc.p

templates = pickle.load(open(sys.argv[1], 'r'))
templates_lenghts = pickle.load(open(sys.argv[2], 'r'))
templates_ulengths = pickle.load(open(sys.argv[3], 'r'))
templates_descriptions = pickle.load(open(sys.argv[4], 'r'))


total_dna = defaultdict(list)
for dna in templates:
    for template in list(set(templates[dna].split(','))):
        total_dna[template].append(dna)

keys = set().union(
    templates_ulengths, templates_descriptions, templates_lenghts, total_dna)
dictionaries = [templates_ulengths,
                templates_descriptions,
                templates_lenghts,
                total_dna]

total = {k: [d.get(k, '') for d in dictionaries] for k in keys}
database = []
for template in total:
    entry = {}
    entry['ulenght'] = total[template][0]
    entry['species'] = total[template][1]
    entry['lengths'] = total[template][2]
    entry['reads'] = total[template][3]
    entry['sequence'] = template
    database.append(entry)
with io.open('../test_data/myDB.json', 'w', encoding='utf-8') as f:
    f.write(unicode(json.dumps(database, ensure_ascii=False)))

# System call
# mongoimport -d Kmers -c complete_genomes_2 --jsonArray < test_data/my_db_seq_2.json
# retcode = subprocess.call([
#     'mongoimport', '-d', 'Kmers', '-c', 'complete_genomes_2', '--jsonArray',
#     '<', 'test_data/my_db_seq_2.json'
# ])
