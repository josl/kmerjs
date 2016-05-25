/* eslint no-underscore-dangle: [2, { "allow": ["_id", "_transform", "_lastLineData", "_flush"] }] */
// let fs = require('browserify-fs');
import Console from 'console';
import readline from 'readline';
import stream from 'stream';
import request from 'request';

import from2 from 'from2';
const fs = require('fs');
import fileReaderStream from 'filereader-stream';

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


let dbs = new Map([
    ['mongo', 'mongodb://localhost:27017/Kmers'],
    ['json', '../test_data/db.json']
]);


class KmerJS {
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
        // this.evalue = evalue || new BN(0.05);
        this.evalue = new BN(0.05);
    }
    stringToMap(string) {
        return objToStrMap(JSON.parse(string));
    }
    objectToMap(object) {
        return objToStrMap(object);
    }
    mapToJSON(strMap) {
        let obj = Object.create(null);
        for (let [k, v] of strMap) {
            // We don’t escape the key '__proto__'
            // which can cause problems on older engines
            obj[k] = v;
        }
        // return obj;
        return JSON.stringify(obj);
    }
}

function readLinesBrowser(kmerObj, readerStream) {
    return new Promise(function (resolve) {
        // Source: https://strongloop.com/strongblog/practical-examples-of-the-new-node-js-streams-api/
        let liner = new stream.Transform({ objectMode: true });
        liner._transform = function (chunk, encoding, done) {
            let data = chunk.toString();
            if (this._lastLineData) {
                data = this._lastLineData + data;
            }

            let lines = data.split('\n');
            this._lastLineData = lines.splice(lines.length - 1, 1)[0];

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
        readerStream().pipe(liner);
        let i = 0;
        let kmerMap = new Map();
        let lines = 0;
        liner.on('readable', function () {
            let line;
            while (null !== (line = liner.read())) {
                if (i === 1 && line.length > 1) {
                    kmers(line, kmerMap, kmerObj.length, kmerObj.preffix, kmerObj.step);
                    kmers(complement(line), kmerMap, kmerObj.length, kmerObj.preffix, kmerObj.step);
                    process.stdout.write(`Lines: ${lines} / Kmers: ${kmerMap.size}\r`);
                } else if (i === 3) {
                    i = -1;
                } else if (i === 2) {
                    if (line !== '+') {
                        Console.log('ERROR!');
                    }
                } else if (i === 0) {
                    if (line[0] !== '@') {
                        Console.log('ERROR2!');
                    }
                }
                i += 1;
                lines += 1;

            }
        });

        liner.on('end', function () {
            process.stdout.write(`\n`);
            resolve(kmerMap);
        });
    });
}

function fromString(string) {
    return from2(function (size, next) {
        // if there's no more content
        // left in the string, close the stream.
        if (string.length <= 0) {
            return next(null, null);
        }
        // Pull in a new chunk of text,
        // removing it from the string.
        var chunk = string.slice(0, size);
        string = string.slice(size);
        // Emit "chunk" from the stream.
        next(null, chunk);
    });
}

export class KmerJSClient extends KmerJS {
    constructor(fastq, preffix = 'ATGAC', length = 16, step = 1, coverage = 1, out = '', db = 'mongo', url = 'http://localhost:3000/kmers') {
        super(fastq, preffix, length, step, coverage, out, db);
        this.url = url;
    }
    findKmers() {
        let that = this;
        return readLinesBrowser(this, function () {
            return fileReaderStream(that.fastq);
        });
    }
    findMatches(kmerQuery) {
        let that = this;
        if (this.db.type === 'json') {
            return findMatchesJSON(kmerQuery, that.coverage);
        }else if (this.db.type === 'mongo') {
            return new Promise(function (resolve, reject) {
                fromString(kmerQuery)
                .pipe(
                    request({url: that.url, method: 'POST'})
                    .on('response', function (response) {
                        resolve(response);
                    })
                    .on('error', function (err) {
                        reject(err);
                    })
                );
            });
        }
    }
}

function readLineNode(that) {
    return new Promise(function (resolve, reject) {
        Console.log('Empezamos...');
        let kmerMap = new Map();
        let rl = readline.createInterface({
            input: fs.createReadStream(that.fastq),
            // output: fs.createWriteStream(that.out),
            terminal: true
        });
        Console.log('let\'s see...');
        let i = 0;
        let lines = 0;
        rl.on('line', function (line) {
            if (i === 1 && line.length > 1) {
                // let kmersBefore = [...kmerMap.keys()];
                kmers(line, kmerMap, that.length, that.preffix, that.step);
                kmers(complement(line), kmerMap, that.length, that.preffix, that.step);
            } else if (i === 3) {
                i = -1;
            }
            i += 1;
            lines += 1;
            process.stdout.write(`Lines: ${lines} / Kmers: ${kmerMap.size}\r`);
        })
        .on('close', function () {
            if (that.out) {
                let out = fs.createWriteStream(that.out);
                out.write('{\n');
                for (let [k, v] of kmerMap) {
                    out.write(`${k}: ${v},`);
                }
                out.write('}\n');
            }
            that.uKmers = kmerMap.size;
            resolve(kmerMap);
        })
        .on('SIGINT', function () {
            reject('SIGINT');
        });
    });
}

export class KmerJSServer extends KmerJS {
    constructor(fastq, preffix = 'ATGAC', length = 16, step = 1, coverage = 1, out = '', db = 'mongo') {
        super(fastq, preffix, length, step, coverage, out, db);
    }
    findKmers() {
        let that = this;
        return readLinesBrowser(that, function () {
            return fs.createReadStream(that.fastq);
        });
        // return readLineNode(that);
    }
    findMatches(kmerMap) {
        let that = this;
        let minScore = 0;
        let url = that.db.url;
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
    }
}
