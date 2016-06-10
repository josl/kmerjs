/* eslint no-underscore-dangle: [2, { "allow": ["_id", "_transform", "_lastLineData", "_flush"] }] */

import BN from 'bignumber.js';
import Console from 'console';
import stream from 'stream';
const fs = require('fs');
import fileReaderStream from 'filereader-stream';
import events from 'events';

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
                coverage = 1, progress = true, env = 'node') {
        this.fastq = fastq;
        this.preffix = preffix;
        this.kmerLength = length;
        this.step = step;
        this.progress = progress;
        this.coverage = coverage;
        this.evalue = new BN(0.05);
        this.kmerMap = new Map(); // [Map object: {16-mer: times found in line}]
        this.env = env;
        if (env === 'browser') {
            this.fileDataRead = 0;
        }
    }
    /**
     * [kmersInLine description]
     * @param  {[string]} line [Sequence read: ATGACCTGAGAGCCTT]
     * @return {[type]}      [description]
     */
    kmersInLine(line) {
        let ini = 0;
        let end = this.kmerLength;
        let stop = line.length - this.kmerLength;
        for (let index = 0; index <= stop; index += 1) {
            let kmer = line.substring(ini, end);
            if (kmer.startsWith(this.preffix)) {
                this.kmerMap.set(kmer, (this.kmerMap.get(kmer) || 0) + 1);
            }
            ini += this.step;
            end = ini + this.kmerLength;
        }
    }
    /**
     * [readFile extract Kmers from file.]
     * @return {[type]} [description]
     */
    readFile() {
        let kmerObj = this;
        let eventEmmitter = new events.EventEmitter();
        let promise = new Promise(function (resolve) {

            // Source: https://strongloop.com/strongblog/practical-examples-of-the-new-node-js-streams-api/
            let liner = new stream.Transform({ objectMode: true });
            liner._transform = function (chunk, encoding, done) {
                let data = chunk.toString();
                if (this._lastLineData) {
                    data = this._lastLineData + data;
                }
                let lines = data.split('\n');
                kmerObj.bytesRead += chunk.length;
                kmerObj.linesPerChunk += lines.length;
                if (kmerObj.env === 'browser') {
                    kmerObj.fileDataRead += chunk.length;
                    // kmerObj.fileDataRead += (kmerObj.lines * kmerObj.bytesRead) / kmerObj.linesPerChunk;
                    eventEmmitter.emit('progress');
                }
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
            kmerObj.lines = 0;
            kmerObj.bytesRead = 0;
            kmerObj.linesPerChunk = 0;
            liner.on('readable', function () {
                let line;
                while (null !== (line = liner.read())) {
                    if (i === 1 && line.length > 1) {
                        [line, complement(line)].forEach(function (kmerLine) {
                            kmerObj.kmersInLine(kmerLine, kmerObj.kmerMap,
                                kmerObj.length,kmerObj.preffix, kmerObj.step);
                        });
                    } else if (i === 3) {
                        i = -1;
                    }
                    i += 1;
                    lines += 1;
                    kmerObj.lines = lines;
                    if (kmerObj.env === 'node' && kmerObj.progress){
                        let progress = `Lines: ${lines} / Kmers: ${kmerObj.kmerMap.size}\r`;
                        process.stdout.write(progress);
                    }
                }
            });
            liner.on('end', function () {
                // Clean up progress output
                if (kmerObj.env === 'node' && kmerObj.progress){
                    process.stdout.write('\n                               \n');
                }
                resolve(kmerObj.kmerMap);
            });
        });
        return {
            promise: promise,
            event: eventEmmitter
        };
    }
}
