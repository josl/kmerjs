/* eslint no-underscore-dangle: [2, { "allow": ["_id"] }] */
const fs = require('fs');
let fileReaderStream = require('filereader-stream');
// let fs = require('browserify-fs');
import Console from 'console';
import readline from 'readline';
import stream from 'stream';

import _ from 'lodash';
import {zScore, fastp, etta} from './stats';
import BN from 'bignumber.js';
import mongoose from 'mongoose';
let complementMap = new Map([
    ['A', 'T'],
    ['T', 'A'],
    ['G', 'C'],
    ['C', 'G']
]);

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

function objToStrMap(obj) {
    let strMap = new Map();
    for (let k of Object.keys(obj)) {
        strMap.set(k, obj[k]);
    }
    return strMap;
}

function jsonToStrMap(jsonStr) {
    return objToStrMap(jsonStr);
}

/**
 * [findKmers description]
 * @param  {[string]} line [Sequence read: ATGACCTGAGAGCCTT]
 * @param  {[Map]} kmerMap [Kmers as keys and number of times seen as value]
 * @param  {[integer]} length  [Lenght of kmer: 16-mer]
 * @param  {[string]} preffix [Sequence to downsample the map: ATGAC]
 * @param  {[integer]} step    [Overlapping kmers of step]
 * @return {[true]}         [description]
 */
export function kmers(line, kmerMap, length, preffix, step) {
    let ini = 0;
    let end = length;
    let stop = line.length - length + 1;
    for (let index = 0; index < stop; index += 1) {
        let key = line.substring(ini, end);
        if (key.startsWith(preffix)) {
            let count = kmerMap.get(key) || 0;
            kmerMap.set(key, count + 1);
        }
        ini += step;
        end = ini + length;
    }
}
/**
 * [standarScore description]
 * @param  {[type]} templateentries    [description]
 * @param  {[type]} templates_ulengths [description]
 * @param  {[type]} template_tot_ulen  [description]
 * @param  {[type]} hits               [description]
 * @return {[type]}                    [description]
 */
// export function standarScore(templateentries, templates_ulengths, template_tot_ulen, hits) {
//     // z = (score - expected)/sqrt(score + expected+etta)
//     // p  = fastp(z)
//     //
//     // If expected < 1 the above poisson approximation is a poor model
//     // Use instead: probabilyty of seing X hits is p**X if probability
//     // of seing one hit is p (like tossing a dice X times)
//     //
//     // if expected <1:
//     //   p = expected**score
//     //
//     //  Comparison of two fractions, Statistical methods in medical research, Armitage et al. p. 125:
//
// }

export function complement(string) {
    return string.replace(/[ATGC]/g, function (match) {
        return complementMap.get(match);
    }).split('').reverse().join('');
}

function findMatchesMongo(kmerMap, url) {
    return new Promise(function (resolve, reject) {
        Console.log('Finding matches...');
        let db = mongoose.connect(url, {
            server: { socketOptions: {keepAlive: 1}},
            replset: { socketOptions: {keepAlive: 1}}
        });
        kmerSchema.index({reads: 1});
        kmerSchema.index({sequence: 1});
        let Kmer = mongoose.model('Kmer', kmerSchema, 'complete_genomes_2');
        let query = {
            reads: {
                $in: Array.from(kmerMap.keys())
            }
        };
        // Get unique and total matches
        Kmer.aggregate({$match: query})
            .unwind('$reads')
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
                        Kmer.aggregate({$match: query})
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
                                    db.disconnect();
                                    resolve({
                                        templates: templates,
                                        hits: hits
                                    });
                                }else {
                                    Console.log(error);
                                    reject('Error found in second aggregate!');
                                }
                            });
                    }else {
                        Console.log('new erroe//');
                        reject('No hits were found!');
                    }
                }else {
                    reject('Error found in first aggregate!');
                }
            });
    });
}

function findMatchesJSON(kmerMap, coverage, url) {
    let db = jsonToStrMap(require(url));
    let hits = 0;
    let matches = new Map();
    let totalMatches = new Map();
    for (let [k, v] of kmerMap) {
        if (v >= coverage) {
            let dbMatch = (db.get(k) || '').split(',');
            let templates = new Set(dbMatch.length > 1 ? dbMatch : []);
            hits += templates.size;
            for (let match of templates) {
                let count = matches.get(match) || 0;
                matches.set(match, count + 1);
                count = totalMatches.get(match) || v - 1;
                totalMatches.set(match, v + 1);
            }
        }
    }
    return {
        templateentries: matches,
        templateentriestot: matches,
        hits: hits
    };
}


function readLineNode(kmerObj) {
    return new Promise(function (resolve, reject) {
        let kmerMap = new Map();
        let rl = readline.createInterface({
            input: fs.createReadStream(kmerObj.fastq),
            // output: fs.createWriteStream(that.out),
            terminal: false
        });
        let i = 0;
        rl.on('line', function (line) {
            if (i === 1 && line.length > 1) {
                kmers(line, kmerMap, kmerObj.length, kmerObj.preffix, kmerObj.step);
                kmers(complement(line), kmerMap, kmerObj.length, kmerObj.preffix, kmerObj.step);
            } else if (i === 3) {
                i = -1;
            }
            i += 1;
        })
        .on('close', function () {
            if (kmerObj.out) {
                let out = fs.createWriteStream(kmerObj.out);
                out.write('{\n');
                for (let [k, v] of kmerMap) {
                    out.write(`${k}: ${v},`);
                }
                out.write('}\n');
            }
            kmerObj.uKmers = kmerMap.size;
            resolve(kmerMap);
        })
        .on('SIGINT', function () {
            reject('SIGINT');
        });
    });
}

function readLineBrowser(kmerObj) {
    return new Promise(function (resolve, reject) {
        // Source: https://strongloop.com/strongblog/practical-examples-of-the-new-node-js-streams-api/
        let liner = new stream.Transform({ objectMode: true });
        liner._transform = function (chunk, encoding, done) {
            var data = chunk.toString();
            if (this._lastLineData) {
                data = this._lastLineData + data;
            }

            var lines = data.split('\n');
            this._lastLineData = lines.splice(lines.length-1,1)[0];

            lines.forEach(this.push.bind(this));
            done();
        };
        liner._flush = function (done) {
            if (this._lastLineData) {
                this.push(this._lastLineData);
            }
            this._lastLineData = null;
            done();
        };
        fileReaderStream(kmerObj.fastq).pipe(liner);
        let i = 0;
        let kmerMap = new Map();
        liner.on('readable', function () {
             let line;
             while (null !== (line = liner.read())) {
                 if (i === 1 && line.length > 1) {
                     kmers(line, kmerMap, kmerObj.length, kmerObj.preffix, kmerObj.step);
                     kmers(complement(line), kmerMap, kmerObj.length, kmerObj.preffix, kmerObj.step);
                 } else if (i === 3) {
                     i = -1;
                 }
                 i += 1;
             }
        });

        liner.on('end', function () {
            resolve(kmerMap);
        });
    });
}

let dbs = new Map([
    ['mongo', 'mongodb://localhost:27017/Kmers'],
    ['json', '../test_data/db.json']
]);


export class KmerJS {
    constructor(fastq, preffix = 'ATGAC', length = 16, step = 1, coverage = 1, out = '', db = 'mongo') {
        this.fastq = fastq;
        this.preffix = preffix;
        this.length = length;
        this.step = step;
        this.out = out === '' ? undefined : out;
        this.coverage = coverage;
        this.uKmers = 0;
        this.db = {
            type: db,
            url: dbs.get(db)
        };
        this.evalue = evalue || new BN(0.05);
    }

    findMatches(kmerMap) {
        let that = this;
        let type = that.db.type;
        let url = that.db.url;
        let minScore = 0;
        if (type === 'mongo'){
            return new Promise(function (resolve, reject) {
                findMatchesMongo(kmerMap, url).then(function (results){
                    Console.log('Computing stats...');
                    let db = mongoose.connect(url);
                    let KmerSummary = mongoose.model('Summary', kmerSummarySchema, 'Summary');
                    return KmerSummary.findOne({templates: {$gt: 1}}).exec(function (error, summary) {
                        if (error){
                            Console.log(error);
                            reject(error);
                        }
                        db.disconnect();
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
                });
            });
            // return [];
            // return findMatchesMongo(kmerMap, url);
        }else if (type === 'json'){
            return findMatchesJSON(kmerMap, this.coverage, url);
        }else {
            throw 'Wrong database source!';
        }

    }
    convertToMap(object){
        return objToStrMap(object);
    }
    findKmers(env = 'node') {
        if (env === 'node'){
            return readLineNode(this);
        }else {
            return readLineBrowser(this);
        }

    }
}
