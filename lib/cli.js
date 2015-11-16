var cli = require('cli');
import {
    kmerjs
}
from '../lib/index';

cli.parse({
    fastq: ['f', 'FASTQ file to parse', 'file', 'test_data/test_short.fastq'],
    preffix: ['p', 'Kmer preffix', 'string', 'ATGAC'],
    lenght: ['l', 'Kmer lenght', 'number', 16],
    step: ['s', 'Kmer step', 'number', 1],
    output: ['o', 'Output file', 'file', 'test_data/out.json']
});

cli.main(function (args, options) {
    kmerjs(options.fastq, options.preffix, options.lenght, options.step, options.output);
});
