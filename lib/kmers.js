/* eslint no-underscore-dangle: [2, { "allow": ["_id", "_transform", "_lastLineData", "_flush"] }] */

import BN from 'bignumber.js';
import Console from 'console';
import stream from 'stream';
const fs = require('fs');
import fileReaderStream from 'filereader-stream';

export let complementMap = new Map([
    ['A', 'T'],
    ['T', 'A'],
    ['G', 'C'],
    ['C', 'G']
]);

function objToStrMap(obj) {
    let strMap = new Map();
    for (let k of Object.keys(obj)) {
        strMap.set(k, obj[k]);
    }
    return strMap;
}

export function jsonToStrMap(jsonStr) {
    return objToStrMap(jsonStr);
}

export function complement(string) {
    return string.replace(/[ATGC]/g, function (match) {
        return complementMap.get(match);
    })
    .split('')
    .reverse()
    .join('');
}

export function stringToMap(string) {
    return objToStrMap(JSON.parse(string));
}
export function objectToMap(object) {
    return objToStrMap(object);
}
export function mapToJSON(strMap) {
    let obj = Object.create(null);
    for (let [k, v] of strMap) {
        // We don’t escape the key '__proto__'
        // which can cause problems on older engines
        obj[k] = v;
    }
    return obj;
}

export class KmerJS {
    /**
     * [constructor description]
     * @param  {[type]} preffix  =             'ATGAC' [Filter kmers starting with]
     * @param  {[type]} length   =             16      [Lenght of kmer: 16-mer]
     * @param  {[type]} step     =             1       [Overlapping kmers of step STEP]
     * @param  {[type]} coverage =             1       [minimun coverage]
     * @param  {[type]} out      =             ''      [description]
     * @param  {[type]} env      =             'node'  [description]
     * @return {[type]}          [description]
     */
    constructor(fastq ='', preffix = 'ATGAC', length = 16, step = 1,
                coverage = 1, progress = false, env = 'node') {
        this.fastq = fastq;
        this.preffix = preffix;
        this.length = length;
        this.step = step;
        this.progress = progress;
        this.coverage = coverage;
        this.uKmers = 0;
        this.evalue = new BN(0.05);
        this.kmerMap = new Map(); // [Map object: {16-mer: times found in line}]
        this.env = env;
    }
    /**
     * [kmersInLine description]
     * @param  {[string]} line [Sequence read: ATGACCTGAGAGCCTT]
     * @return {[type]}      [description]
     */
    kmersInLine(line) {
        let ini = 0;
        let end = this.length;
        let stop = line.length - this.length + 1;
        for (let index = 0; index < stop; index += 1) {
            let key = line.substring(ini, end);
            if (key.startsWith(this.preffix)) {
                this.kmerMap.set(key, (this.kmerMap.get(key) || 0) + 1);
            }
            ini += this.step;
            end = ini + this.length;
        }
    }
    /**
     * [readFile extract Kmers from file.]
     * @return {[type]} [description]
     */
    readFile() {
        var kmerObj = this;
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
            if (kmerObj.env === 'node'){
                fs.createReadStream(kmerObj.fastq).pipe(liner);
            }else if (kmerObj.env === 'browser') {
                fileReaderStream(kmerObj.fastq).pipe(liner);
            }
            let i = 0;
            let lines = 0;
            liner.on('readable', function () {
                let line;
                while (null !== (line = liner.read())) {
                    if (i === 1 && line.length > 1) {
                        [line, complement(line)].forEach(function (line) {
                            kmerObj.kmersInLine(line, kmerObj.kmerMap,
                                kmerObj.length,kmerObj.preffix, kmerObj.step);
                        });
                        let progress = `Lines: ${lines} / Kmers: ${kmerObj.kmerMap.size}\r`;
                        if (kmerObj.env === 'node' && kmerObj.progress){
                            process.stdout.write(progress);
                        }else if (kmerObj.env === 'browser') {
                            // Console.log(progress);
                        }
                    } else if (i === 3) {
                        i = -1;
                    }
                    i += 1;
                    lines += 1;

                }
            });

            liner.on('end', function () {
                kmerObj.uKmers = kmerObj.kmerMap.size;
                resolve(kmerObj.kmerMap);
            });
        });
    }

}
