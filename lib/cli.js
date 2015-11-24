var cli = require('cli');
import {
    KmerJS
} from '../lib/index';
import Console from 'console';

cli.parse({
    fastq: ['f', 'FASTQ file to parse', 'file', 'test_data/test_long.fastq'],
    preffix: ['p', 'Kmer preffix', 'string', 'ATGAC'],
    lenght: ['l', 'Kmer lenght', 'number', 16],
    step: ['s', 'Kmer step', 'number', 1],
    coverage: ['c', 'Min coverage', 'number', 1],
    output: ['o', 'Output file', 'file', 'test_data/out.json']
});

cli.main(function (args, options) {
    let kmerjs = new KmerJS(
        options.fastq, options.preffix, options.lenght, options.step,
        options.coverage, options.output, 'mongo');
    kmerjs.findKmers().then(function (kmerMap) {
        kmerjs.findMatches(kmerMap).then(function (output) {
            Console.log('Template\tScore\tExpected\tz\tp_value\tquery\tcoverage [%]\ttemplate coverage [%]\tdepth\tKmers in Template\tDescription');
            output.forEach(function (data) {
                let seq = data.get('template');
                let score = data.get('score');
                let expec = data.get('expected');
                let z = data.get('z');
                let p = data.get('probability');
                let fracQ = data.get('frac-q');
                let fracD = data.get('frac-d');
                let cov = data.get('coverage');
                let ulen = data.get('ulength');
                let spec = data.get('species');
                Console.log(`${seq}\t${score}\t${expec}\t${z}\t${p}\t${fracQ}\t${fracD}\t${cov}\t${ulen}\t${spec}\t`);
            });
        });
    });
});
