/* eslint no-underscore-dangle: [2, { "allow": ["_id", "_transform", "_lastLineData", "_flush"] }] */
import Console from 'console';
import request from 'request';
import BN from 'bignumber.js';
import from2 from 'from2';

import {
    KmerJS,
    jsonToStrMap,
    mapToJSON
} from './kmers.js';

import {zScore, fastp, etta} from './stats';


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
    let kmerQuerySize = kmerObject.kmerMapSize;
    let sequenceHit = kmerObject.firstMatches.get(sequence);
    let originalUScore = sequenceHit.uScore;
    let originalTScore = sequenceHit.tScore;
    let matchUScore = match.uScore;
    if (match.uScore > minScore){
        let z = zScore(match.uScore, match.ulength, results.hits, summary.uniqueLens);
        let probability = fastp(z).times(summary.templates);
        let allow = kmerObject.evalue.cmp(probability); // p <= evalue
        if (allow >= 0) {
            let fracQ = new BN(100).times(2).times(matchUScore).dividedBy(
                new BN(kmerQuerySize).plus(etta)
            );
            let fracD = new BN(100).times(matchUScore).dividedBy(
                new BN(match.ulength).plus(etta)
            );
            let totFracQ = new BN(100)
                .times(2)
                .times(originalUScore)
                .dividedBy(
                    new BN(kmerQuerySize).plus(etta)
                );
            let totFracD = new BN(100)
                .times(originalUScore)
                .dividedBy(
                    new BN(match.ulength).plus(etta)
                );
            let totFracCov = new BN(originalTScore)
                .dividedBy(match.lengths).round(2, 6).toNumber();
            let expected = new BN(results.hits)
                .times(match.ulength)
                .dividedBy(summary.uniqueLens);
            return new Map([
                ['template', sequence],
                ['score', matchUScore],
                ['expected', expected.round(0, 6).toNumber()],
                ['z', z.round(2).toNumber()],
                ['probability', probability.toNumber()],
                ['frac-q', fracQ.round(2, 6).toNumber()],
                ['frac-d', fracD.round(2,6).toNumber()],
                ['depth', new BN(match.tScore).dividedBy(match.lengths).round(2,6).toNumber()],
                ['kmers-template', match.ulength],
                ['total-frac-q', totFracQ.round(2, 6).toNumber()],
                ['total-frac-d', totFracD.round(2,6).toNumber()],
                ['total-temp-cover', totFracCov],
                ['species', match.species]
            ]);
        }
    }
}

/**
 * [sortKmerMatches Sort Matches by Hits (= Score)]
 * @param  {[type]} a [Match]
 * @param  {[type]} b [Match]
 * @return {[type]}   [description]
 */
function sortKmerMatches(a, b) {
    if (a[1].uScore> b[1].uScore) {
        return -1;
    }
    if (a[1].uScore < b[1].uScore) {
        return 1;
    }
    // a must be equal to b
    return 0;
}

export class KmerFinderClient extends KmerJS {
    constructor(fastq, env, preffix = 'ATGAC', length = 16, step = 1,
                coverage = 1, out = true, db = 'server',
                url = 'http://localhost:3000/kmers',
                summary = '/Users/cisneror/code/genomic-git/kmerjs/test_data/summary.json',
                collection = 'genomes',
                dbName = 'Kmers') {
        super(fastq, preffix, length, step, coverage, out, env);
        this.dbLocation = db;
        this.dbURL = url;
        this.collection = collection;
        this.dbName = dbName;
        this.maxHits = 100;
    }
    findKmers() {
        return this.readFile();
    }
    findFirstMatch(kmerQuery) {
        console.log(this);
        let that = this;
        let promise =  new Promise(function (resolve, reject) {
            kmerQuery.set('db', that.dbName);
            kmerQuery.set('collection', that.collection);
            console.log(kmerQuery.get('db'), kmerQuery.get('collection'));
            fromString(JSON.stringify(mapToJSON(kmerQuery)))
                .pipe(
                    request
                        .post(that.dbURL)
                        .on('response', function (response) {
                             if (response.statusCode === 201 || response.statusCode === 200 || response.statusCode === 202) {
                                 console.log('resolving!');
                                 let winnerData = '';
                                 let totalData = 0;
                                 response
                                     .on('data', function(chunk) {
                                         winnerData += chunk.toString();
                                         totalData += chunk.length;
                                     })
                                     .on('end', function() {
                                         let winner = JSON.parse(winnerData);
                                         console.log('Total data of reduced DB (bytes) ', totalData);
                                         winner.templates = jsonToStrMap(winner.templates);
                                         winner.templates.forEach(function (hit, sequence) {
                                             hit.kmers = new Set(hit.kmers);
                                         });
                                         console.log('resolving!');
                                         resolve(winner);
                                     });
                            }else if (response.statusCode === 204){
                                console.log('we get an empty dataset!');
                                reject('No hits were found!');
                            }else{
                                reject('error');
                            }
                        })
                        .on('error', function (err) {
                            console.log('we get an error!', err);
                            reject(err);
                        })
                );
        });
        return promise;
    }
    findMatches(winner, kmerMap) {
        let hitCounter = 0;
        let notFound = true;
        // let kmerResults = [];
        let kmerObject = this;
        function findWinner(results) {
            console.log('Let\'s find the winner!');
            let templates = [...results.templates].sort(sortKmerMatches);
            if (hitCounter === 0) {
                kmerObject.firstMatches = results.templates;
                if (kmerObject.progress) {
                    let out = 'Template\tScore\tExpected\tz\tp_value\tquery\tcoverage [%]\ttemplate coverage [%]\tdepth\tKmers in Template\tDescription\n';
                    process.stdout.write(out);
                }
            }
            let sequence = templates[0][0];
            let match = templates[0][1];
            let winner = matchSummary(kmerObject, sequence, match, results, kmerObject.summary);

            // Winner is undefined if match's score < minScore
            if (winner && kmerObject.evalue.cmp(winner.get('probability')) >= 0){
                hitCounter += 1;
                if (kmerObject.progress) {
                    let seq = winner.get('template');
                    let score = winner.get('score');
                    let expec = winner.get('expected');
                    let z = winner.get('z');
                    let p = winner.get('probability');
                    let fracQ = winner.get('frac-q');
                    let fracD = winner.get('frac-d');
                    let cov = winner.get('depth');
                    let ulen = winner.get('kmers-template');
                    let spec = winner.get('species');
                    let out = `${seq}\t${score}\t${expec}\t${z}\t${p}\t${fracQ}\t${fracD}\t${cov}\t${ulen}\t${spec}\n`;
                    process.stdout.write(out);
                }
                return {
                    kmers: match.kmers,
                    match: winner
                };
            }else {
                notFound = false;
                return;
            }
        }

        function removeWinnerKmers(matchKmers){
            if (matchKmers){
                // remove all kmers in best hit from kmerMap
                matchKmers.forEach(function(kmer){
                    kmerMap.delete(kmer);
                });
                return;
            }else{
                return;
            }
        }

        function getMatches (kmerObject, kmerQueryMap) {
            let templates = new Map();
            let nHits = 0;

            kmerObject.firstMatches.forEach(function (hit, sequence) {
                let template = templates.get(sequence);
                for (const kmer of hit.kmers) {
                    if (kmerQueryMap.has(kmer)){
                        let kmerCoverage = kmerQueryMap.get(kmer);
                        if (template !== undefined){
                            template.tScore += kmerCoverage;
                            template.uScore += 1;
                            template.kmers.add(kmer);
                        }else {
                            templates.set(sequence, {
                                tScore: kmerCoverage,
                                uScore: 1,
                                lengths: hit.lengths,
                                ulength: hit.ulength,
                                species: hit.species,
                                kmers: new Set([kmer])
                            });
                            template = templates.get(sequence);
                        }
                    }
                }
                if (template !== undefined){
                    nHits += template.kmers.size;
                }else{
                    kmerObject.firstMatches.delete(sequence);
                }
            });
            if (nHits === 0) {
                throw new Error('No hits were found! (nHits === 0)');
            }
            return {
                templates: templates,
                hits: nHits
            };
        }

        function * loop() {
            while (notFound && hitCounter < kmerObject.maxHits){
                // Find new matches from first matches.
                let results = getMatches(kmerObject, kmerMap);
                let winner = findWinner(results);
                if (winner){
                    removeWinnerKmers(winner.kmers);
                    yield mapToJSON(winner.match);
                }
            }
            if (hitCounter === 0){
                throw new Error('No hits were found! (kmerResults.length === 0)');
            }
        }
        kmerObject.summary = winner.summary;
        kmerObject.firstMatches = winner.templates;
        return loop();
    }
}
