#!/usr/bin/env node

let cli = require('cli');
import {
    KmerFinderServer
} from './kmerFinderServer';
import Console from 'console';

cli.parse({
    fastq: ['f', 'FASTQ file to parse', 'file', 'test_data/test_long.fastq'],
    preffix: ['p', 'Kmer preffix', 'string', 'ATGAC'],
    length: ['l', 'Kmer lenght', 'number', 16],
    step: ['s', 'Kmer step', 'number', 1],
    coverage: ['c', 'Min coverage', 'number', 1],
    output: ['o', 'Print info', 'number', 1],
    program: ['P', 'Program to execute: [findKmers, findMatches]', 'string', 'findMatches'],
    score: ['S', 'Score to execute: [standard, winner]', 'string', 'winner'],
    database: ['d', 'Database to query: [Bacteria, Virus...]', 'string', 'KmerMap'],
    url: ['u', 'Url of the database', 'string', 'mongodb://localhost:27017/Kmers']
});

cli.main(function (args, options) {
    let kmerjs = new KmerFinderServer(
        options.fastq, options.preffix, options.length,
        options.step, options.coverage, options.output,
        'mongo', options.url,
        options.database, options.score
    );
    if (options.program === 'findKmers') {
        let kmers = kmerjs.findKmers();
        kmers.then(function (kmerMap) {
            // let keys = [...kmerMap.keys()];
            // Console.log('\n' + 'Unique Kmers: ', keys.length);
            // Console.log('Let\'s look at the top 10');
            // let i = 0;
            // let sortedArray = [...kmerMap]
            //     .sort(function(a, b) {
            //         return b[1] - a[1];
            //     });
            // for (let [k, v] of sortedArray) {
            //     // We don’t escape the key '__proto__'
            //     // which can cause problems on older engines
            //     Console.log(k,v);
            //     i += 1;
            //     if (i === 10) {break;}
            // }
            kmerjs.close();
            process.exit();
        });
    }else if (options.program === 'findMatches') {
        kmerjs.findKmers()
            .then(function (kmerMap) {
                let keys = [...kmerMap.keys()];
                Console.log('Kmers: ', keys.length);
                let matches = kmerjs.findMatches(kmerMap);
                matches.event
                    .on('winner', function (winner) {
                        let seq = winner.get('template');
                        let score = winner.get('score');
                        let expec = winner.get('expected');
                        let z = winner.get('z');
                        let p = winner.get('probability');
                        let fracQ = winner.get('frac-q');
                        let fracD = winner.get('frac-d');
                        let cov = winner.get('depth');
                        let ulen = winner.get('kmers-template');
                        let spec = winner.get('species');
                        let out = `${seq}\t${score}\t${expec}\t${z}\t${p}\t${fracQ}\t${fracD}\t${cov}\t${ulen}\t${spec}\n`;
                        process.stdout.write(out);
                    });
                matches.promise.then(function () {
                        kmerjs.close();
                        process.exit();
                    })
                    .catch(function (err) {
                        console.log(err);
                        kmerjs.close();
                        process.exit();
                    });
            });
    }else{
        Console.log(options.program + ' is not a valid option! [findKmers, findMatches]');
    }

});
