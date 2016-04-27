/* eslint no-underscore-dangle: [2, { "allow": ["_id", "_transform", "_lastLineData", "_flush"] }] */
// let fs = require('browserify-fs');

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
    // return JSON.stringify(obj);
}

export class KmerJS {
    constructor(preffix = 'ATGAC', length = 16, step = 1,
                coverage = 1, out = '', env = 'node') {
        this.preffix = preffix;
        this.length = length;
        this.step = step;
        this.out = out === '' ? undefined : out;
        this.coverage = coverage;
        this.uKmers = 0;
        // this.evalue = evalue || new BN(0.05);
        this.evalue = new BN(0.05);
        this.kmerMap = new Map();
        this.env = env;
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
    kmersInLine(line) {
        let ini = 0;
        let end = this.length;
        let stop = line.length - this.length + 1;
        for (let index = 0; index < stop; index += 1) {
            let key = line.substring(ini, end);
            if (key.startsWith(this.preffix)) {
                let count = this.kmerMap.get(key) || 0;
                this.kmerMap.set(key, count + 1);
            }
            ini += this.step;
            end = ini + this.length;
        }
    }
    readLines() {
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
                        kmerObj.kmersInLine(line, kmerObj.kmerMap, kmerObj.length,
                                     kmerObj.preffix, kmerObj.step);
                        kmerObj.kmersInLine(complement(line), kmerObj.kmerMap,
                                     kmerObj.length, kmerObj.preffix, kmerObj.step);
                        let progress = `Lines: ${lines} / Kmers: ${kmerObj.kmerMap.size}\r`;
                        if (kmerObj.env === 'node' && kmerObj.out !== undefined){
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
                resolve(kmerObj.kmerMap);
            });
        });
    }

}
