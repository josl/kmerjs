/* eslint no-underscore-dangle: [2, { "allow": ["_id"] }] */
import fs from 'fs';
import Console from 'console';
import readline from 'readline';
import _ from 'lodash';

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

function findMatchesJSON(kmerMap, coverage, url) {
    let db = jsonToStrMap(require(url));
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

export class KmerJS {
    constructor(fastq, preffix, length, step, coverage, out, db) {
        this.fastq = fastq;
        this.preffix = preffix;
        this.length = length;
        this.step = step;
        this.out = out;
        this.coverage = coverage;
        if (db === 'mongo'){
            this.db = {
                type: db,
                url: 'mongodb://localhost:27017/Kmers'
            };
        }else if (db === 'json'){
            this.db = {
                type: db,
                url: '../test_data/db.json'
            };
        }
    }

    findMatches(kmerMap) {
        let type = this.db.type;
        let url = this.db.url;
        if (type === 'mongo'){
            return findMatchesMongo(kmerMap, url);
        }else if (type === 'json'){
            return findMatchesJSON(kmerMap, this.coverage, url);
        }else {
            throw 'Wrong database source!';
        }

    }

    findKmers() {
        let kmerMap = new Map();
        let that = this;
        return new Promise(function (resolve, reject) {
            let rl = readline.createInterface({
                input: fs.createReadStream(that.fastq),
                // output: fs.createWriteStream(that.out),
                terminal: false
            });
            let i = 0;
            rl.on('line', function (line) {
                if (i === 1 && line.length > 1) {
                    kmers(line, kmerMap, that.length, that.preffix, that.step);
                    kmers(complement(line), kmerMap, that.length, that.preffix, that.step);
                } else if (i === 3) {
                    i = -1;
                }
                i += 1;
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
                resolve(kmerMap);
            })
            .on('SIGINT', function () {
                reject(false);
            });
        });
    }

}
