/* eslint no-underscore-dangle: [2, { "allow": ["_id"] }] */
import fs from 'fs';
import Console from 'console';
import readline from 'readline';
import _ from 'lodash';

import mongoose from 'mongoose';
let complementMap = new Map([
    ['A', 'T'],
    ['T', 'A'],
    ['G', 'G'],
    ['C', 'C']
]);
let url = 'mongodb://localhost:27017/Kmers';

let Schema = mongoose.Schema;

let kmerSchema = new Schema({
    lengths: Number,
    ulengths: Number,
    description: String,
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
export function findKmers(line, kmerMap, length, preffix, step) {
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
    });
}

export class KmerJS {
    constructor(fastq, preffix, length, step, out) {
        this.fastq = fastq;
        this.preffix = preffix;
        this.length = length;
        this.step = step;
        this.out = out;
    }
}

export function findMatches(kmerMap) {
    // Use connect method to connect to the Server
    return new Promise(function (resolve, reject) {
        mongoose.connect(url);
        kmerSchema.index({reads: 1});
        let Kmer = mongoose.model('Kmer', kmerSchema, 'complete_genomes_2');
        let query = {
            reads: {
                $in: Array.from(kmerMap.keys())
            }
        };
        Kmer.aggregate({$match: query})
            .unwind('$reads')
            .match(query)
            .group({
                _id: '$sequence',
                uScore: { $sum: 1}
            })
            .sort({uScore: -1})
            .exec(function (err, result) {
                if (err === null){
                    let hits = _.reduce(result, function (total, n) {
                        return total + n.uScore;
                    }, 0);
                    if (hits !== 0){
                        let templateentries = new Map();
                        for (let temp of result){
                            templateentries.set(temp._id, temp.uScore);
                        }
                        resolve({
                            templateentries: templateentries,
                            hits: hits
                        });
                        mongoose.disconnect();
                    }else {
                        Console.log('new erroe//');
                        reject('No hits were found!');
                    }
                }else {
                    reject('Error found!');
                }
                mongoose.disconnect();
            });
    });
}

export function findMatchesJSON(kmerMap, coverage) {
    let db = jsonToStrMap(require('../test_data/db.json'));
    let hits = 0;
    let matches = new Map();
    let totalMatches = new Map();
    for (let [k, v] of kmerMap) {
        if (v >= coverage) {
            let dbMatch = (db.get(k) || '')
                .split(',');
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


export function kmerjs(fastq, preffix, length, step, out) {
    let kmerMap = new Map();
    return new Promise(function (resolve, reject) {
        let rl = readline.createInterface({
            input: fs.createReadStream(fastq),
            output: fs.createWriteStream(out)
        });
        let i = 0;
        rl.on('line', function (line) {
            if (i === 1) {
                findKmers(line, kmerMap, length, preffix, step);
                findKmers(complement(line), kmerMap, length, preffix, step);
            } else if (i === 4) {
                i = 0;
            }
            i += 1;
        })
        .on('close', function () {
            if (out !== null) {
                rl.write('{\n');
                for (let [k, v] of kmerMap) {
                    rl.write(`${k}: ${v},`);
                }
                rl.write('}\n');
            }
            resolve(kmerMap);
        })
        .on('SIGINT', function () {
            reject(false);
        });
    });
}

// export default function (fastq, preffix, length) {
//     return new Promise(function (resolve, reject) {
//         let init = moment();
//         let stream = fs.createReadStream(fastq);
//         let lines = 0;
//         let rl = readline.createInterface({
//           input: stream
//         });
//         let i = 0;
//         rl.on('line', function(line) {
//             if (i === 1) {
//                 let blocks = Math.ceil(line.length / length);
//                 let ini = 0;
//                 let end = length;
//                 for (let b =0; b < blocks; b++) {
//                     let key = line.substring(ini, end);
//                     let count = kmer_map.get(key) || 0;
//                     kmer_map.set(key, count + 1);
//                     ini = b + length;
//                     end = ini + length;
//                 }
//             } else if (i === 4){
//                 i = 0;
//             }
//             i += 1;
//             lines += 1;
//         }).on('close', function() {
//             let end = moment();
//             console.log('Time passed: ', end.diff(init), ' ms.');
//             resolve(kmer_map.size);
//         }).on('SIGINT', function() {
//             reject(false);
//         });
//     });
// }
