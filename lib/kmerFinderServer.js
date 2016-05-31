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

function extractKmers() {
    return mongoose.model('Bacteria', kmerMapSchema, 'Bacteria')
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

function findKmersMatches(kmerMap, url, collection) {
    kmerMapSchema.index({reads: 1});
    kmerMapSchema.index({'templates.sequence': 1});
    let kmerDB = mongoose.model(collection, kmerMapSchema, collection);
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
            .unwind('templates')
            .group({
                _id: {template: '$templates'},
                filteredReads:{
                    $push: '$kmer'
                }
            })
            .project({
                _id: 0,
                template: '$_id.template.sequence',
                lengths: '$_id.template.lengths',
                ulengths: '$_id.template.ulengths',
                species: '$_id.template.species',
                filteredReads: '$filteredReads'
            })
            .allowDiskUse(true)
            .cursor({ batchSize: 1000 })
            .exec();

            cursor.toArray(function(err, hits) {
                let templates = new Map();
                let nHits = 0;
                for (let hit of hits){
                    nHits += hit.filteredReads.length;
                    hit.filteredReads.forEach(function (read) {
                        let template = templates.get(hit.template);
                        let kmerCoverage = kmerMap.get(read);
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
                                reads: hit.filteredReads
                            });
                        }
                    });
                }
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
function findMatchesMongoAggregation(kmerMap, url, collection, progress) {
    kmerSchema.index({reads: 1});
    kmerSchema.index({sequence: 1});
    let kmerDB = mongoose.model(collection, kmerSchema, collection);
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
            filteredReads: {
                $filter: {
                    input: "$reads",
                    as: "read",
                    cond:{
                        $setIsSubset: [
                            {
                                $map: {
                                    input: {
                                        "$literal": ["single_element"]
                                    },
                                    as: "el",
                                    in: "$$read"
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
                nHits += hit.filteredReads.length;
                hit.filteredReads.forEach(function (read) {
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
                            reads: hit.filteredReads
                        });
                    }
                });
            }
            if (nHits === 0) {
                throw new Error('No hits were found!');
            }
            let out = `Hits: ${nHits} / Templates: ${hits.length}\r`;
            if (progress) {
                process.stdout.write(out);
            }
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
            _id: { sequence: '$sequence', read: '$reads'},
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
    if (match.uScore > minScore){
        let z = zScore(match.uScore, match.ulength, results.hits, summary.uniqueLens);
        let probability = fastp(z).times(summary.templates);
        let allow = kmerObject.evalue.cmp(probability); // p <= evalue
        if (allow >= 0) {
            let fracQ = new BN(100).times(2).times(match.uScore).dividedBy(
                new BN(kmerObject.kmerMap.size).plus(etta)
            );
            let fracD = new BN(100).times(match.uScore).dividedBy(
                new BN(match.ulength).plus(etta)
            );
            let totFracQ = new BN(100)
                .times(2)
                .times(kmerObject.firstMatches.get(sequence).uScore)
                .dividedBy(
                    new BN(kmerObject.kmerMap.size).plus(etta)
                );
            let totFracD = new BN(100)
                .times(kmerObject.firstMatches.get(sequence).uScore)
                .dividedBy(
                    new BN(match.ulength).plus(etta)
                );
            let totFracCov = new BN(kmerObject.firstMatches.get(sequence).tScore)
                .dividedBy(match.lengths).round(2).toNumber();
            let expected = new BN(results.hits)
                .times(match.ulength)
                .dividedBy(summary.uniqueLens);
            return new Map([
                ['template', sequence],
                ['score', match.uScore],
                ['expected', expected.round(0, 1).toNumber()],
                ['z', z.round(2).toNumber()],
                ['probability', probability.toNumber()],
                ['frac-q', fracQ.round(2, 2).toNumber()],
                ['frac-d', fracD.round(2).toNumber()],
                ['depth', new BN(match.tScore).dividedBy(match.lengths).round(2).toNumber()],
                ['kmers-template', match.ulength],
                ['total-frac-q', totFracQ.round(2, 2).toNumber()],
                ['total-frac-d', totFracD.round(2).toNumber()],
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
 * [winnerScoringRecursive]
 * @param  {[Object]} summary    [Summary query from DB]
 * @param  {[Map]}    kmerMap    [Dictionary of presence of Kmers]
 * @param  {[Object]} kmerObject [Instance of KmerServer]
 * @return {[type]}              [description]
 */
function winnerScoringRecursive(summary, kmerMap, kmerObject) {
    let hitCounter = 0;
    let notFound = true;
    let maxHits = 100;
    let kmerResults = [];


    function findWinner(results) {
        let templates = [...results.templates].sort(sortKmerMatches);
        if (hitCounter === 0) {
            kmerObject.firstMatches = results.templates;
            if (kmerObject.progress) {
                let out = 'Template\tScore\tExpected\tz\tp_value\tquery\tcoverage [%]\ttemplate coverage [%]\tdepth\tKmers in Template\tDescription\n'
                process.stdout.write(out);
            }
        }
        let winner = matchSummary(
            kmerObject, templates[0][0],
            templates[0][1], results, summary
        );
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
            return results.templates.get(winner.get('template')).reads;
        }else {
            return;
        }
    }

    function removeWinnerKmers(matchReads){
        if (matchReads){
            // remove all kmers in best hit from kmerMap
            matchReads.forEach(function(kmer){
                kmerMap.delete(kmer);
            });
            return;
        }else{
            notFound = false;
            return;
        }
    }

    function getMatches (kmerObject, kmerMap) {
        let templates = new Map();
        let nHits = 0;
        let queryKmers = [...kmerMap.keys()];

        kmerObject.firstMatches.forEach((hit, sequence) => {
            let filteredReads = _.intersection(hit.reads, queryKmers);
            nHits += filteredReads.length;
            filteredReads.forEach(function (read) {
                let template = templates.get(sequence);
                let kmerCoverage = kmerMap.get(read);
                if (template !== undefined){
                    template.tScore += kmerCoverage;
                    template.uScore += 1;
                }else {
                    templates.set(sequence, {
                        tScore: kmerCoverage,
                        uScore: 1,
                        lengths: hit.lengths,
                        ulength: hit.ulength,
                        species: hit.species,
                        reads: filteredReads
                    });
                }
            });
        });
        if (nHits === 0) {
            throw new Error('No hits were found!');
        }
        // let out = `Hits: ${nHits} / Templates: ${templates.size}\r`;
        // if (kmerObject.progress) {
        //     process.stdout.write(out);
        // }
        return {
            templates: templates,
            hits: nHits
        };
    }

    var loop = function() {
        if (!(notFound && hitCounter < maxHits)){
            if (kmerResults.length === 0){
                throw new Error('No hits were found!');
            }
            return kmerResults.sort(sortKmerResults);
        }else{
            // Find new matches from first matches.
            let results = getMatches(kmerObject, kmerMap);
            let winnerKmers = findWinner(results);
            removeWinnerKmers(winnerKmers);
            return loop();
        }
    };
    return findKmersMatches(kmerMap, kmerObject.db.url, kmerObject.collection,
                            kmerObject.progress)
        .then(findWinner)
        .then(removeWinnerKmers)
        .then(loop);
}
/**
 * [winnerScoring Winner takes All scoring scheme]
 * @param  {[type]} kmerObject [KmerJS Object]
 * @param  {[type]} kmerMap    [DNA sequence in Kmer space]
 * @return {[type]}            [Promise]
 */
function winnerScoring(kmerObject, kmerMap) {
    return mongoose
            .model('Summary', kmerSummarySchema, 'Summary')
            .findOne({templates: {$gt: 1}}, { _id: 0})
            .then(function (summary) {
                return winnerScoringRecursive(summary, kmerMap, kmerObject);
            });
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
                progress = true, db = 'mongo', url='mongodb:\/\/localhost:27017/Kmers', collection='genomes', method = 'standard') {
        super(fastq, preffix, length, step, coverage, progress, 'node');
        this.db = {
            connection: mongoose.connect(url, {
                server: { socketOptions: {keepAlive: 120}},
                replset: { socketOptions: {keepAlive: 120}}
            }),
            type: db
        };
        this.method = method;
        this.collection = collection;
        this.firstMatches = new Map();
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
        this.db.connection.disconnect();
    }
}
