/* eslint no-underscore-dangle: [2, { "allow": ["_id", "_transform", "_lastLineData", "_flush"] }] */
// let fs = require('browserify-fs');
import Console from 'console';

import _ from 'lodash';
import {zScore, fastp, etta} from './stats';
import BN from 'bignumber.js';
import mongoose from 'mongoose';

import {
    KmerJS,
} from './kmers.js';

let Schema = mongoose.Schema;

let kmerSchema = new Schema({
    lengths: Number,
    ulengths: Number,
    description: String,
    reads: Array
});

let kmerSummarySchema = new Schema({
    templates: Number,
    uniqueLens: Number,
    totalLen: Number,
    reads: Array
});


function findMatchesMongo(kmerMap, url, collection) {
    return new Promise(function (resolve, reject) {
        // let db = mongoose.connect(url, {
        //     server: { socketOptions: {keepAlive: 1}},
        //     replset: { socketOptions: {keepAlive: 1}}
        // });
        kmerSchema.index({reads: 1});
        kmerSchema.index({sequence: 1});
        let kmerDB = mongoose.model('Kmer', kmerSchema, collection);
        let query = {
            reads: {
                $in: Array.from(kmerMap.keys())
            }
        };
        // Get unique and total matches
        let cursor = kmerDB.aggregate([{$match: query}, {$limit: 30}]);

        cursor.unwind('$reads')
            .match(query)
            .group({
                _id: { sequence: '$sequence', read: '$reads'},
                uScore: { $sum: 1}
            })
            // .sort({uScore: -1})
            .exec(function (err, matches) {
                if (err === null){
                    let hits = _.reduce(matches, function (total, n) {
                        return total + n.uScore;
                    }, 0);
                    if (hits !== 0){
                        let templates = new Map();
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
                        kmerDB.aggregate({$match: query})
                            .project({id: '$sequence', lengths: 1, ulenght: 1, species: 1, _id: 0})
                            .exec(function (error, sequences) {
                                if (error === null){
                                    for (let seq of sequences){
                                        let vals = templates.get(seq.id);
                                        vals.lengths = seq.lengths;
                                        vals.ulength = seq.ulenght;
                                        vals.species = seq.species;
                                        templates.set(seq.id, vals);
                                    }
                                    // db.disconnect();
                                    resolve({
                                        templates: templates,
                                        hits: hits
                                    });
                                }else {
                                    // db.disconnect();
                                    reject(error.message);
                                }
                            });
                    }else {
                        // db.disconnect();
                        reject('No hits were found!');
                    }
                }else {
                    // db.disconnect();
                    reject(err.message);
                }
            });
    });
}



export class KmerFinderServer extends KmerJS {
    constructor(fastq, preffix = 'ATGAC', length = 16, step = 1, coverage = 1,
                out = '', db = 'mongo',
                url='mongodb://localhost:27017/Kmers', collection='genomes') {
        super(preffix, length, step, coverage, out, 'node');
        this.fastq = fastq;
        // this.db = {
        //     type: db,
        //     url: url
        // };
        this.db = {
            connection: mongoose.connect(url, {
                server: { socketOptions: {keepAlive: 120}},
                replset: { socketOptions: {keepAlive: 120}}
            }),
            type: db
        };
        this.collection = collection;
    }
    findKmers() {
        // this.db.connection.disconnect();
        return this.readLines();
    }
    findMatches(kmerMap) {
        let minScore = 0;
        let url = this.db.url;
        let that = this;
        return new Promise(function (resolve, reject) {
            findMatchesMongo(kmerMap, url, that.collection).then(function (results){
                // let db = mongoose.connect(url);
                let KmerSummary = mongoose.model('Summary', kmerSummarySchema, 'Summary');
                return KmerSummary.findOne({templates: {$gt: 1}}).exec(function (error, summary) {
                    if (error){
                        Console.log(error);
                        reject(error);
                    }
                    // that.db.connection.close();
                    let kmerResults = [];
                    for (let [sequence, match] of results.templates){
                        if (match.uScore > minScore){
                            let z = zScore(match.uScore, match.ulength, results.hits, summary.uniqueLens);
                            let probability = fastp(z).times(summary.templates);
                            let allow = that.evalue.cmp(probability); // p <= evalue
                            if (allow >= 0) {
                                let fracQ = new BN(100).times(2).times(match.uScore).dividedBy(new BN(that.uKmers).plus(etta));
                                let fracD = new BN(100).times(match.uScore).dividedBy(new BN(match.ulength).plus(etta));
                                // Console.log(fracQ);
                                kmerResults.push(new Map([
                                    ['template', sequence],
                                    ['score', match.uScore],
                                    ['expected', results.hits * match.ulength / summary.uniqueLens],
                                    ['z', z.toNumber()],
                                    ['probability', probability.toNumber()],
                                    ['frac-q', fracQ.toNumber()],
                                    ['frac-d', fracD.toNumber()],
                                    ['coverage', match.tScore / match.lengths],
                                    ['ulength', match.ulength],
                                    // ['total-frac-q', 100 * 2 * match.uScore / (that.uKmers + etta)],
                                    // ['total-frac-d', 100 * match.uScore / (match.ulength + etta)],
                                    // ['total-coverage', match.tScore / match.lengths],
                                    ['species', match.species]
                                ]));
                            }
                        }
                    }
                    let sortedResults = kmerResults.sort(function (a, b) {
                        if (a.get('score') > b.get('score')) {
                            return -1;
                        }
                        if (a.get('score') < b.get('score')) {
                            return 1;
                        }
                        // a must be equal to b
                        return 0;
                    });
                    resolve(sortedResults);
                });
            }, function(err) {
              console.log(err); // Error: "It broke"
            });
        });
    }
}
