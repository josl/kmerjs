/* eslint no-underscore-dangle: [2, { "allow": ["_id", "_transform", "_lastLineData", "_flush"] }] */
// let fs = require('browserify-fs');
import Console from 'console';
// import stream from 'stream';
import request from 'request';
import BN from 'bignumber.js';
import from2 from 'from2';
// import fileReaderStream from 'filereader-stream';

import {
    KmerJS,
    jsonToStrMap,
    mapToJSON
} from './kmers.js';

import {zScore, fastp, etta} from './stats';


// let dbs = new Map([
//     ['mongo', 'mongodb://localhost:27017/Kmers'],
//     ['json', '../test_data/db.json']
// ]);

// function readLinesBrowser(kmerObj, readerStream) {
//     Console.log(kmerObj.fastq);
//     return new Promise(function (resolve) {
//         // Source: https://strongloop.com/strongblog/practical-examples-of-the-new-node-js-streams-api/
//         let liner = new stream.Transform({ objectMode: true });
//         Console.log(liner);
//         liner._transform = function (chunk, encoding, done) {
//             let data = chunk.toString();
//             if (this._lastLineData) {
//                 data = this._lastLineData + data;
//             }
//
//             let lines = data.split('\n');
//             this._lastLineData = lines.splice(lines.length - 1, 1)[0];
//
//             lines.forEach(this.push.bind(this));
//             done();
//         };
//         liner._flush = function (done) {
//             if (this._lastLineData) {
//                 this.push(this._lastLineData);
//             }
//             this._lastLineData = null;
//             done();
//         };
//         Console.log('LALALALALA');
//         readerStream(kmerObj.fastq).pipe(liner);
//         Console.log('LOLOLOLOLO');
//         let i = 0;
//         let lines = 0;
//         liner.on('readable', function () {
//             let line;
//             Console.log(line);
//             while (null !== (line = liner.read())) {
//                 if (i === 1 && line.length > 1) {
//                     kmerObj.kmersInLine(line, kmerObj.kmerMap, kmerObj.length,
//                                  kmerObj.preffix, kmerObj.step);
//                     kmerObj.kmersInLine(complement(line), kmerObj.kmerMap,
//                                  kmerObj.length, kmerObj.preffix, kmerObj.step);
//                     process.stdout.write(`Lines: ${lines} /
//                                           Kmers: ${kmerObj.kmerMap.size}\r`);
//                 } else if (i === 3) {
//                     i = -1;
//                 } else if (i === 2) {
//                     if (line !== '+') {
//                         Console.log('ERROR!');
//                     }
//                 } else if (i === 0) {
//                     if (line[0] !== '@') {
//                         Console.log('ERROR2!');
//                     }
//                 }
//                 i += 1;
//                 lines += 1;
//
//             }
//         });
//         Console.log('LELELELELELLE');
//         liner.on('end', function () {
//             process.stdout.write(`\n`);
//             resolve(kmerObj.kmerMap);
//         });
//     });
// }

// FIXME: BROWSER DATABASE NEEDS TO BE FIXED!!!!
function computeStats(results, query) {
    let summary = jsonToStrMap(require(query.db.summary));
    let kmerResults = [];
    let minScore = 0;
    Console.log('Computing stats...', results);
    for (let [sequence, match] of results.templates){
        Console.log(match.uScore);
        if (match.uScore > minScore){
            let z = zScore(match.uScore, match.ulength, results.hits, summary.uniqueLens);
            let probability = fastp(z).times(summary.templates);
            let allow = query.evalue.cmp(probability); // p <= evalue
            if (allow >= 0) {
                let fracQ = new BN(100).times(2).times(match.uScore).dividedBy(new BN(query.uKmers).plus(etta));
                let fracD = new BN(100).times(match.uScore).dividedBy(new BN(match.ulength).plus(etta));
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
    return sortedResults;
}

function findMatchesJSON(kmerMap, query) {
    let kmerDB = jsonToStrMap(require(query.db.url));
    let hits = 0;
    let matches = new Map();
    let totalMatches = new Map();
    for (let [k, v] of kmerMap) {
        if (v >= query.coverage) {
            let dbMatch = (kmerDB.get(k) || '').split(',');
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
    return computeStats({
        templates: matches,
        templateentriestot: totalMatches,
        hits: hits
    }, query);
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

export class KmerFinderClient extends KmerJS {
    constructor(fastq, env, preffix = 'ATGAC', length = 16, step = 1,
                coverage = 1, out = '', db = 'mongo',
                url = 'http://localhost:3000/kmers',
                summary = '/Users/cisneror/code/genomic-git/kmerjs/test_data/summary.json') {
        super(preffix, length, step, coverage, out, env);
        this.fastq = fastq;
        this.db = {
            type: db,
            url: url,
            summary: summary
        };
    }
    findKmers() {
        return this.readLines();
    }
    findMatches(kmerQuery) {
        let that = this;
        return new Promise(function (resolve, reject) {
            if (that.db.type === 'json') {
                resolve(findMatchesJSON(kmerQuery, that));
            }else if (that.db.type === 'mongo') {
                    fromString(JSON.stringify(mapToJSON(kmerQuery)))
                    .pipe(
                        request({url: that.db.url, method: 'POST'})
                        .on('response', function (response) {
                            resolve(response);
                        })
                        .on('error', function (err) {
                            reject(err);
                        })
                    );
            }
        });
    }
}
