/* eslint no-underscore-dangle: [2, { "allow": ["_id", "_transform", "_lastLineData", "_flush"] }] */
// let fs = require('browserify-fs');
import Console from 'console';

import _ from 'lodash';
import {zScore, fastp, etta} from './stats';
import BN from 'bignumber.js';
BN.config({ ROUNDING_MODE: 2 }) ;
import mongoose from 'mongoose';
import {
    KmerJS,
} from './kmers.js';

let Schema = mongoose.Schema;

let kmerSchema = new Schema({
    lengths: Number,
    ulengths: Number,
    description: String,
    reads: Array,
    sequence: String
});

let kmerSummarySchema = new Schema({
    templates: Number,
    uniqueLens: Number,
    totalLen: Number
});

/**
 * [findMatchesMongo description]
 * @param  {[type]} kmerMap    [Dictionary of presence of Kmers]
 * @param  {[type]} url        [description]
 * @param  {[type]} collection [description]
 * @return {[type]}            [description]
 */
function findMatchesMongo(kmerMap, url, collection) {
    console.log('Let\'s find some matches!');
    kmerSchema.index({reads: 1});
    kmerSchema.index({sequence: 1});
    let kmerDB = mongoose.model(collection, kmerSchema, collection);
    let query = {
        reads: {
            $in: Array.from(kmerMap.keys())
        }
    };
    // Get unique and total matches
    let templates = new Map();
    let hits;
    return kmerDB
            // Query to return matched templates
            .aggregate([{$match: query}])
            .unwind('$reads')
            .match(query)
            .group({
                _id: { sequence: '$sequence', read: '$reads'},
                uScore: { $sum: 1}
            })
            .exec()
            .then(function(matches){
                console.log('Executed!');
                hits = _.reduce(matches, function (total, n) {
                    return total + n.uScore;
                }, 0);
                if (hits !== 0){
                    for (let temp of matches){
                        let read = temp._id.read;
                        let sequence = temp._id.sequence;
                        let match = templates.get(sequence);
                        let tCount = match ? match.tScore : 0;
                        let uCount = match ? match.uScore : 0;
                        tCount += kmerMap.get(read);
                        uCount += 1;
                        templates.set(sequence, {
                            tScore: tCount, uScore: uCount
                        });
                    }
                    // Get unique, total lenghts & species per match
                    query = {
                        sequence: {
                            $in: Array.from(templates.keys())
                        }
                    };
                    let projection = {
                        id: '$sequence', lengths: 1, ulenght: 1,
                        species: 1, _id: 0
                    };
                    console.log('More queries!')
                    return kmerDB
                            .aggregate({$match: query})
                            .project(projection)
                            .exec();
                }else {
                    throw new Error('No hits were found!');
                }
            })
            .then(function(sequences){
                for (let seq of sequences){
                    let vals = templates.get(seq.id);
                    vals.lengths = seq.lengths;
                    vals.ulength = seq.ulenght;
                    vals.species = seq.species;
                    templates.set(seq.id, vals);
                }
                return {
                    templates: templates,
                    hits: hits
                };
            })
            .catch(function(error){
                throw new Error(error);
            });
}

/**
 * [matchSummary description]
 * @param  {[type]} kmerObject [description]
 * @param  {[type]} sequence   [description]
 * @param  {[type]} match      [description]
 * @param  {[type]} results    [description]
 * @param  {[type]} summary    [description]
 * @return {[type]}            [description]
 */
function matchSummary(kmerObject, sequence, match, results, summary){
    let minScore = 0;
    if (match.uScore > minScore){
        let z = zScore(match.uScore, match.ulength, results.hits, summary.uniqueLens);
        let probability = fastp(z).times(summary.templates);
        let allow = kmerObject.evalue.cmp(probability); // p <= evalue
        if (allow >= 0) {
            let fracQ = new BN(100).times(2).times(match.uScore).dividedBy(
                new BN(kmerObject.uKmers).plus(etta)
            );
            let fracD = new BN(100).times(match.uScore).dividedBy(
                new BN(match.ulength).plus(etta)
            );
            // Console.log(fracQ);
            return new Map([
                ['template', sequence],
                ['score', match.uScore],
                ['expected', results.hits * match.ulength / summary.uniqueLens],
                ['z', z.toNumber()],
                ['probability', probability.toNumber()],
                ['frac-q', fracQ.toNumber()],
                ['frac-d', fracD.toNumber()],
                ['coverage', match.tScore / match.lengths],
                ['ulength', match.ulength],
                // ['total-frac-q', 100 * 2 * match.uScore / (kmerObject.uKmers + etta)],
                // ['total-frac-d', 100 * match.uScore / (match.ulength + etta)],
                // ['total-coverage', match.tScore / match.lengths],
                ['species', match.species]
            ]);
        }
    }
}

/**
 * [sortKmerResults description]
 * @param  {[type]} a [Match Summary entry]
 * @param  {[type]} b [description]
 * @return {[type]}   [description]
 */
function sortKmerResults(a, b) {
    if (a.get('score') > b.get('score')) {
        return -1;
    }
    if (a.get('score') < b.get('score')) {
        return 1;
    }
    // a must be equal to b
    return 0;
}

/**
 * [standardScoring description]
 * @param  {[Object]} kmerObject [Instance of KmerServer]
 * @param  {[Map]}    kmerMap    [Dictionary of presence of Kmers]
 * @return {[type]}              [description]
 */
function standardScoring(kmerObject, kmerMap) {
    return Promise.all([
        mongoose
            .model('Summary', kmerSummarySchema, 'Summary')
            .findOne({templates: {$gt: 1}}),
        findMatchesMongo(kmerMap, kmerObject.db.url, kmerObject.collection),
    ])
    .then(([summary, results]) => {
        let kmerResults = [];
        for (let [sequence, match] of results.templates){
            kmerResults.push(
                matchSummary(kmerObject, sequence, match, results, summary)
            );
        }
        return kmerResults.sort(sortKmerResults);
    })
    .catch(function(error) {
        // Receives first rejection among the Promises
        throw new Error(error);
    });
}

function myRecursive() {

}

/**
 * [winnerScoring Winner takes All scoring scheme]
 * @param  {[type]} kmerObject [KmerJS Object]
 * @param  {[type]} kmerMap    [DNA sequence in Kmer space]
 * @return {[type]}            [Promise]
 */
function winnerScoring(kmerObject, kmerMap) {
    let minScore = 0;
    let maxhits = 100;
    let hitcounter = 1
    let finishedProcessing = false;

    while (!finishedProcessing && hitcounter < maxhits){
        findMatchesMongo(kmerMap, kmerObject.db.url, kmerObject.collection)
            .then(function () {

            })
        hitcounter++;
    }

    return new Promise(function (resolve, reject) {

    });
}

export class KmerFinderServer extends KmerJS {
    /**
     * [constructor description]
     * @param  {[type]} fastq                                 [description]
     * @param  {[type]} preffix                               =             'ATGAC'   [description]
     * @param  {[type]} length                                =             16        [description]
     * @param  {[type]} step                                  =             1         [description]
     * @param  {[type]} coverage                              =             1         [description]
     * @param  {[type]} out                                   =             ''        [description]
     * @param  {[type]} db                                    =             'mongo'   [description]
     * @param  {[type]} method                                =             'standar' [description]
     * @param  {[type]} url='mongodb://localhost:27017/Kmers' [description]
     * @param  {[type]} collection='genomes'                  [description]
     * @return {[type]}                                       [description]
     */
    constructor(fastq, preffix = 'ATGAC', length = 16, step = 1, coverage = 1,
                progress = true, db = 'mongo',
                url='mongodb:\/\/localhost:27017/Kmers', collection='genomes', method = 'standard') {
        super(fastq, preffix, length, step, coverage, progress, 'node');
        this.db = {
            connection: mongoose.connect(url, {
                server: { socketOptions: {keepAlive: 120}},
                replset: { socketOptions: {keepAlive: 120}}
            }),
            type: db
        };
        this.method = method;
        this.collection = collection;
    }
    /**
     * [findKmers Wrapper around reading file function]
     * @return {[Map]} [Kmer Map]
     */
    findKmers() {
        return this.readFile();
    }
    /**
     * [findMatches description]
     * @param  {[Map]}     kmerMap [Kmer Map]
     * @return {[Promise]}         [description]
     */
    findMatches(kmerMap) {
        if (this.method === 'standard'){
            return standardScoring(this, kmerMap);
        }else if (this.method === 'winner') {
            return winnerScoring(this, kmerMap);
        }else{
            throw new Error('Scoring scheme unknown');
        }
    }
}
