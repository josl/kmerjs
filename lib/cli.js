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
    console.log(options.url);
    let kmerjs = new KmerFinderServer(
        options.fastq, options.preffix, options.length,
        options.step, options.coverage, options.output,
        'mongo', options.url,
        options.database, options.score
    );
    if (options.program === 'findKmers') {
        let kmers = kmerjs.findKmers();
        kmers.then(function () {
            process.exit();
        });
    }else if (options.program === 'findMatches') {
        kmerjs.findKmers()
            .then(function (kmerMap) {
                let keys = [...kmerMap.keys()];
                Console.log('Kmers: ', keys.length);
                kmerjs.findMatches(kmerMap)
                    .then(function () {
                        process.exit();
                    })
                    .catch(function (err) {
                        console.log(err);
                        process.exit();
                    });
            });
    }else{
        Console.log(options.program + ' is not a valid option! [findKmers, findMatches]');
    }

});
