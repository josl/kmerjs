/* eslint no-underscore-dangle: [2, { "allow": ["_id", "_transform", "_lastLineData", "_flush", "emit"] }] */
// let fs = require('browserify-fs');
import Console from 'console';

import _ from 'lodash';
import {zScore, fastp, etta} from './stats';
import BN from 'bignumber.js';
BN.config({ ROUNDING_MODE: 2 }) ;
import mongoose from 'mongoose';
import {
    KmerJS,
} from './kmers.js';
import events from 'events';

let Schema = mongoose.Schema;

let kmerSchema = new Schema({
    lengths: Number,
    ulengths: Number,
    description: String,
    reads: Array,
    sequence: String
});

let kmerSummarySchema = new Schema({
    templates: Number,
    uniqueLens: Number,
    totalLen: Number
});

let kmerMapSchema = new Schema({
    kmer: String,
    templates: Array
});
/**
 * [extractKmers
     db.genomes.aggregate([
      {$match: {}},
      {$unwind: '$reads'},
      {$group: {
          _id: {read: '$reads'},
          templates: {
            $push: {
              sequence: '$sequence',
              lengths: '$lengths',
              ulengths: '$ulenght',
              species: '$species'
            }
          }
        }
      },
      {$project: {
          _id: 0, kmer: '$_id.read', templates: '$templates'
        }
      },
      {$out: 'KmerBacteria'}
    ], {
      allowDiskUse: true
    })
 * ]
 * @param  {[type]} kmerObject [description]
 * @return {[type]}            [description]
 */
function extractKmers(kmerObject) {
    return kmerObject.conn.model('Bacteria', kmerMapSchema, 'Bacteria')
        .aggregate([{
            $match: {
            }
        }])
        .unwind('reads')
        .group({
            _id: {read: '$reads'},
            templates: {
                $push: {
                    sequence: '$sequence',
                    lengths: '$lengths',
                    ulengths: '$ulenght',
                    species: '$species',
                }
            }
        })
        .project({
            _id: 0, kmer: '$_id.read', templates: '$templates'
        })
        .out('KmerBacteria')
        .allowDiskUse(true)
        .exec();
}

function createSummary(kmerObject) {
    kmerObject.conn.model.aggregate([
        { $unwind: "$templates" },
        { $group: {
            _id: {sequence: {"templates": "$templates.sequence"}},
            count: {$sum:1}
        }},
        {$group: {
            _id: null,
            count: {$sum:1}
        }
    }]);
    kmerObject.conn.model.aggregate([
        { $unwind: "$templates"},
        { $group: {
            _id: null,
            count: {$sum:"$templates.lengths"}
        }
    }]);
    kmerObject.conn.model.aggregate([
        { $unwind: "$templates"},
        { $group: {
            _id: null,
            count: {$sum:"$templates.ulengths"}
        }
    }]);
    return;
}

function findMatchesMapReduce(kmerMap, url, collection) {
    kmerMapSchema.index({kmer: 1});
    kmerMapSchema.index({sequence: 1});
    let kmerDB = mongoose.model(collection, kmerSchema, collection);
    let mapReduce = {};
    mapReduce.map = `function () {
        var matches = 0;
        for (var i = 0; i < kmerQuery.length; i++) {
            if (this.reads.indexOf(kmerQuery[i]) !== -1) {
                matches += 1;
            }
        }
        emit('test', matches);
    }`;
    mapReduce.reduce = `function (key, values) {
            return Array.sum(values);
    }`;
    mapReduce.scope = {
        kmerQuery: [...kmerMap.keys()]
    };
    mapReduce.query = {
        reads: {
            // The $in operator selects the documents where the value of a field
            // equals any value in the specified array
            $in: [...kmerMap.keys()]
        }
    };
    return kmerDB.mapReduce(mapReduce, function(err, data) {
       console.log(err, data);
    });
}

function updateTemplates(templates, hit, template, kmerCoverage){
    if (template !== undefined){
        template.tScore += kmerCoverage;
        template.uScore += 1;
    }else {
        templates.set(hit.template, {
            tScore: kmerCoverage,
            uScore: 1,
            lengths: hit.lengths,
            ulength: hit.ulengths,
            species: hit.species,
            kmers: new Set(hit.filteredKmers)
        });
    }
}

function findKmersMatches(kmerMap, conn, collection) {
    kmerMapSchema.index({reads: 1});
    kmerMapSchema.index({'templates.sequence': 1});
    let kmerDB = conn.model(collection, kmerMapSchema, collection);
    let kmerQuery = [...kmerMap.keys()];

    return new Promise(
        function (resolve, reject) {
        let cursor = kmerDB
            .aggregate([{
                $match: {
                    kmer: {
                        $in: kmerQuery,
                        // $type: 2 // BSON type (2) = String
                    }
                }
            }])
            .unwind('templates') // FIXME: Look into not unwind & postprocess in client
            .group({             // Apparently is faster than the alternative...
                _id: {template: '$templates'},
                filteredKmers:{
                    $push: '$kmer'
                }
            })
            .project({
                _id: 0,
                template: '$_id.template.sequence',
                lengths: '$_id.template.lengths',
                ulengths: '$_id.template.ulengths',
                species: '$_id.template.species',
                filteredKmers: '$filteredKmers'
            })
            .allowDiskUse(true)
            .cursor({ batchSize: 3000 })
            .exec();

            cursor.toArray(function(err, hits) {
                if (err) {
                    reject(err);
                }
                let templates = new Map();
                let nHits = 0;
                 for (let hit of hits){
                    nHits += hit.filteredKmers.length;
                    for (var i = 0; i < hit.filteredKmers.length; i++) {
                        updateTemplates(
                            templates, hit, templates.get(hit.template),
                            kmerMap.get(hit.filteredKmers[i])
                        );
                    }
                }
                console.log(nHits);
                if (nHits === 0) {
                    reject('No hits were found!');
                }
                resolve( {
                    templates: templates,
                    hits: nHits
                });
            });
        });
}


function findKmersMatches2(kmerMap, conn, collection) {
    kmerMapSchema.index({reads: 1});
    kmerMapSchema.index({'templates.sequence': 1});
    let kmerDB = conn.model(collection, kmerMapSchema, collection);
    let kmerQuery = [...kmerMap.keys()];

    return new Promise(
        function (resolve, reject) {
        let cursor = kmerDB
            .aggregate([{
                $match: {
                    kmer: {
                        $in: kmerQuery,
                        // $type: 2 // BSON type (2) = String
                    }
                }
            }])
            // .unwind('templates') // FIXME: Look into not unwind & postprocess in client
            // .group({             // Apparently is faster than the alternative...
            //     _id: {template: '$templates'},
            //     filteredKmers:{
            //         $push: '$kmer'
            //     }
            // })
            // .project({
            //     _id: 0,
            //     template: '$templates.sequence',
            //     lengths: '$templates.lengths',
            //     ulengths: '$templates.ulengths',
            //     species: '$templates.species',
            //     kmer: '$kmer'
            // })
            .allowDiskUse(true)
            .cursor({ batchSize: 3000 })
            .exec();

            cursor.toArray(function(err, kmerHits) {
                if (err) {
                    reject(err);
                }
                let templates = new Map();
                let nHits = 0;
                // let templateHits = new Map();
                kmerHits.forEach(function (hit) {
                    let kmerCoverage = kmerMap.get(hit.kmer);
                    hit.templates.forEach(function (kmerTemplate) {
                        let template = templates.get(kmerTemplate.sequence);
                        if (template !== undefined){
                            template.tScore += kmerCoverage;
                            template.uScore += 1;
                            template.kmers.add(hit.kmer);
                        }else {
                            templates.set(kmerTemplate.sequence, {
                                tScore: kmerCoverage,
                                uScore: 1,
                                lengths: kmerTemplate.lengths,
                                ulength: kmerTemplate.ulengths,
                                species: kmerTemplate.species,
                                kmers: new Set([hit.kmer])
                            });
                            template = templates.get(kmerTemplate.sequence);
                        }
                        nHits += 1;
                    });
                });

                //  for (let hit of hits){
                //     nHits += hit.filteredKmers.length;
                //     for (var i = 0; i < hit.filteredKmers.length; i++) {
                //         updateTemplates(
                //             templates, hit, templates.get(hit.template),
                //             kmerMap.get(hit.filteredKmers[i])
                //         );
                //     }
                // }
                console.log(nHits);
                if (nHits === 0) {
                    reject('No hits were found!');
                }
                resolve( {
                    templates: templates,
                    hits: nHits
                });
            });
        });
}


function findMatchesMongoAggregation(kmerMap, conn, collection) {
    kmerSchema.index({reads: 1});
    kmerSchema.index({sequence: 1});
    let kmerDB = conn.model(collection, kmerMapSchema, collection);
    let kmerQuery = [...kmerMap.keys()];

    return kmerDB
        .aggregate([{
            $match: {
                reads: {
                        $in: kmerQuery
                }
            }
        }])
        .project({
            _id: 0, sequence: 1, lengths: 1, ulength: "$ulenght", species: 1,
            filteredKmers: {
                $filter: {
                    input: "$reads",
                    as: "kmer",
                    cond:{
                        $setIsSubset: [
                            {
                                $map: {
                                    input: {
                                        "$literal": ["single_element"]
                                    },
                                    as: "el",
                                    in: "$$kmer"
                                }
                            },
                            kmerQuery
                        ]
                    }
                }
            }
        })
        .exec()
        .then(function (hits) {
            let templates = new Map();
            let nHits = 0;
            for (let hit of hits){
                nHits += hit.filteredKmers.length;
                hit.filteredKmers.forEach(function (read) {
                    let template = templates.get(hit.sequence);
                    let kmerCoverage = kmerMap.get(read);
                    if (template !== undefined){
                        template.tScore += kmerCoverage;
                        template.uScore += 1;
                    }else {
                        templates.set(hit.sequence, {
                            tScore: kmerCoverage,
                            uScore: 1,
                            lengths: hit.lengths,
                            ulength: hit.ulength,
                            species: hit.species,
                            kmers: new Set(hit.filteredKmers)
                        });
                    }
                });
            }
            if (nHits === 0) {
                throw new Error('No hits were found!');
            }
            let out = `Hits: ${nHits} / Templates: ${hits.length}\r`;
            return {
                templates: templates,
                hits: nHits
            };
        });
}

/**
 * [findMatchesMongo description]
 * @param  {[type]} kmerMap    [Dictionary of presence of Kmers]
 * @param  {[type]} url        [description]
 * @param  {[type]} collection [description]
 * @return {[type]}            [description]
 */
function findMatchesMongo(kmerMap, url, collection, limit=0) {
    kmerSchema.index({reads: 1});
    kmerSchema.index({sequence: 1});
    let kmerDB = mongoose.model(collection, kmerSchema, collection);
    let query = {
        reads: {
            // The $in operator selects the documents where the value of a field
            // equals any value in the specified array
            $in: [...kmerMap.keys()]
        }
    };
    // Get unique and total matches
    let templates = new Map();
    let hits;
    // Query to return matched templates
    let cursor = kmerDB
        .aggregate([{$match: query}]);
    return cursor
        .unwind('reads')
        .match(query)
        .group({
            _id: { sequence: '$sequence', kmer: '$reads'},
            uScore: { $sum: 1}
        })
        .exec()
        .then(function(matches){
            hits = _.reduce(matches, function (total, n) {
                return total + n.uScore;
            }, 0);
            console.log('Hits found: ', hits);
            if (hits !== 0){
                for (let temp of matches){
                    let kmer = temp._id.kmer;
                    let sequence = temp._id.sequence;
                    let match = templates.get(sequence);
                    let tCount = match ? match.tScore : 0;
                    let uCount = match ? match.uScore : 0;
                    tCount += kmerMap.get(kmer);
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
                let projection = {
                    id: '$sequence', lengths: 1, ulenght: 1,
                    species: 1, _id: 0
                };
                return kmerDB
                        .aggregate({$match: query})
                        .project(projection)
                        .exec();
            }else {
                throw new Error('No hits were found!');
            }
        })
        .then(function(sequences){
            for (let seq of sequences){
                let vals = templates.get(seq.id);
                vals.lengths = seq.lengths;
                vals.ulength = seq.ulenght;
                vals.species = seq.species;
                templates.set(seq.id, vals);
            }
            return {
                templates: templates,
                hits: hits
            };
        });
}
/**
 * [findKmersTemplate Queries DB for matche's kmers]
 * @param  {[String]} template [Sequence template]
 * @return {[Array]}           [Array of kmers]
 */
function findKmersTemplate(template, collection){
    // let kmerDB = mongoose.model(collection, kmerMapSchema, collection);
    // return kmerDB.find({'templates.sequence': template}).exec();
    let kmerDB = mongoose.model(collection, kmerSchema, collection);
    return kmerDB.findOne({sequence: template}, {reads: 1}).exec();
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
 * [sortKmerResults description]
 * @param  {[type]} a [Match Summary entry]
 * @param  {[type]} b [description]
 * @return {[type]}   [description]
 */
function sortKmerResults(a, b) {
    if (a.get('score') > b.get('score')) {
        return -1;
    }
    if (a.get('score') < b.get('score')) {
        return 1;
    }
    // a must be equal to b
    return 0;
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

/**
 * [winnerScoring]
 * @param  {[Map]}    kmerMap    [Dictionary of presence of Kmers]
 * @param  {[Object]} kmerObject [Instance of KmerServer]
 * @return {[type]}              [description]
 */
function winnerScoring(kmerObject, kmerMap) {
    let hitCounter = 0;
    let notFound = true;
    let kmerResults = [];
    let firstMatch = kmerObject.conn
            .model('Summary', kmerSummarySchema, 'Summary')
            .findOne({templates: {$gt: 1}}, { _id: 0})
            .then(function (summary) {
                kmerObject.summary = summary;
                return findKmersMatches(kmerMap, kmerObject.conn,
                        kmerObject.collection);
            });

    function findWinner(results) {
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
            kmerResults.push(winner);
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
            return match.kmers;
        }else {
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
            notFound = false;
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

    var loop = function() {
        if (!(notFound && hitCounter < kmerObject.maxHits)){
            if (kmerResults.length === 0){
                throw new Error('No hits were found! (kmerResults.length === 0)');
            }
            console.log('exiting...', kmerObject.kmerMapSize);
            console.log('kmerMapSize...', kmerObject.kmerMapSize);
            return kmerResults;
        }else{
            // Find new matches from first matches.
            let results = getMatches(kmerObject, kmerMap);
            let winnerKmers = findWinner(results);
            removeWinnerKmers(winnerKmers);
            return loop();
        }
    };

    return firstMatch
            .then(findWinner)
            .then(removeWinnerKmers)
            .then(loop);
}

/**
 * [standardScoring description]
 * @param  {[Object]} kmerObject [Instance of KmerServer]
 * @param  {[Map]}    kmerMap    [Dictionary of presence of Kmers]
 * @return {[type]}              [description]
 */
function standardScoring(kmerObject, kmerMap) {
    let kmerResults = [];
    return Promise.all([
        mongoose
            .model('Summary', kmerSummarySchema, 'Summary')
            .findOne({templates: {$gt: 1}}, { _id: 0}).exec(),
        findMatchesMongoAggregation(kmerMap, kmerObject.db.url, kmerObject.collection, kmerObject.progress),
    ])
    .then(([summary, results]) => {
        kmerObject.firstMatches = results.templates;
        for (let [sequence, match] of results.templates){
            kmerResults.push(
                matchSummary(kmerObject, sequence, match, results, summary)
            );
        }
        return kmerResults.sort(sortKmerResults);
    });
}

export class KmerFinderServer extends KmerJS {
    /**
     * [constructor description]
     * @param  {[type]} fastq                                 [description]
     * @param  {[type]} preffix                               =             'ATGAC'   [description]
     * @param  {[type]} length                                =             16        [description]
     * @param  {[type]} step                                  =             1         [description]
     * @param  {[type]} coverage                              =             1         [description]
     * @param  {[type]} out                                   =             ''        [description]
     * @param  {[type]} db                                    =             'mongo'   [description]
     * @param  {[type]} method                                =             'standar' [description]
     * @param  {[type]} url='mongodb://localhost:27017/Kmers' [description]
     * @param  {[type]} collection='genomes'                  [description]
     * @return {[type]}                                       [description]
     */
    constructor(fastq, preffix = 'ATGAC', length = 16, step = 1, coverage = 1,
                progress = true, db = 'mongo',
                url='mongodb:\/\/localhost:27017/Kmers', collection='genomes',
                method = 'standard', maxHits = 100) {
        super(fastq, preffix, length, step, coverage, progress, 'node');
        this.conn = mongoose.createConnection(url, {
            server: { socketOptions: {keepAlive: 120}},
            replset: { socketOptions: {keepAlive: 120}}
        });
        this.method = method;
        this.collection = collection;
        this.firstMatches = new Map();
        this.maxHits = maxHits;
    }
    /**
     * [findKmers Wrapper around reading file function]
     * @return {[Map]} [Kmer Map]
     */
    findKmers() {
        return this.readFile().promise;
    }
    /**
     * [findMatches description]
     * @param  {[Map]}     kmerMap [Kmer Map]
     * @return {[Promise]}         [description]
     */
    findMatches(kmerMap) {
        if (this.method === 'standard'){
            return standardScoring(this, kmerMap);
        }else if (this.method === 'winner') {
            return winnerScoring(this, kmerMap);
        }else{
            throw new Error('Scoring scheme unknown');
        }
    }
    findMatchesTest(kmerMap) {
        return winnerScoring(this, kmerMap);
    }
    close() {
        this.conn.close();
    }
}
