var cli = require('cli');
import {
    KmerFinderServer
} from '../lib/kmerFinderServer';
import Console from 'console';

cli.parse({
    fastq: ['f', 'FASTQ file to parse', 'file', 'test_data/test_long.fastq'],
    preffix: ['p', 'Kmer preffix', 'string', 'ATGAC'],
    length: ['l', 'Kmer lenght', 'number', 16],
    step: ['s', 'Kmer step', 'number', 1],
    coverage: ['c', 'Min coverage', 'number', 1],
    output: ['o', 'Print info', 'number', 1],
    program: ['P', 'Program to execute: [findKmers, findMatches]', 'string', 'findKmers'],
    score: ['S', 'Score to execute: [standard, winner]', 'string', 'winner']
});

cli.main(function (args, options) {
    let kmerjs = new KmerFinderServer(
        options.fastq, options.preffix, options.length,
        options.step, options.coverage, options.output,
        'mongo', 'mongodb://localhost:27017/Kmers',
        'complete_genomes_4', options.score
    );
    if (options.program === 'findKmers') {
        kmerjs.findKmers()
            .then(function (kmerMap) {
                let keys = [...kmerMap.keys()];
                Console.log('\n' + 'Unique Kmers: ', keys.length);
                Console.log('Let\'s look at the top 10');
                let i = 0;
                let sortedArray = [...kmerMap]
                    .sort(function(a, b) {
                        return b[1] - a[1];
                    });
                for (let [k, v] of sortedArray) {
                    // We don’t escape the key '__proto__'
                    // which can cause problems on older engines
                    Console.log(k,v);
                    i += 1;
                    if (i === 10) {break;}
                }
                kmerjs.db.connection.close();
            });
    }else if (options.program === 'findMatches') {
        kmerjs.findKmers()
            .then(function (kmerMap) {
                let keys = [...kmerMap.keys()];
                Console.log('Kmers: ', keys.length);
                kmerjs.findMatches(kmerMap)
                    .then(function (output) {
                        Console.log('Template\tScore\tExpected\tz\tp_value\tquery\tcoverage [%]\ttemplate coverage [%]\tdepth\tKmers in Template\tDescription');
                        output.forEach(function (data) {
                            let seq = data.get('template');
                            let score = data.get('score');
                            let expec = data.get('expected');
                            let z = data.get('z');
                            let p = data.get('probability');
                            let fracQ = data.get('frac-q');
                            let fracD = data.get('frac-d');
                            let cov = data.get('depth');
                            let ulen = data.get('kmers-template');
                            let spec = data.get('species');
                            Console.log(`${seq}\t${score}\t${expec}\t${z}\t${p}\t${fracQ}\t${fracD}\t${cov}\t${ulen}\t${spec}\t`);
                        });
                        kmerjs.db.connection.close();
                        kmerjs.db.connection.close();
                    });
            });
    }else{
        Console.log(options.program + ' is not a valid option! [findKmers, findMatches]');
    }

});
