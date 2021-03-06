import {
    KmerFinderServer,
} from '../lib/kmerFinderServer';
import {
    jsonToStrMap,
} from '../lib/kmers';
import chai from 'chai';
import chaiAsPromised from 'chai-as-promised';
import chaiJsonEqual from 'chai-json-equal';
import console from 'console';
// let kmer_json = require('test_data/test_long.json');
import kmer_json from '../test_data/test_long.json';
chai.use(chaiJsonEqual);
chai.use(chaiAsPromised);
chai.should();

// describe('kmerFinderServer findMatches test_short', function () {
//     this.timeout(500000000);
//
//     it('Winner takes All match NC_017625 (Escherichia coli DH1)', function () {
//         let kmerObj = new KmerFinderServer(
//             './test_data/test_short.fastq',
//             'ATGAC', 16, 1, 1, true, 'mongo', 'mongodb://localhost:27017/Kmers',
//             'complete_genomes_4', 'winner'
//         );
//         // let json = require('.././test_data/db_short_results.json');
//         return kmerObj.findKmers().then(function (kmers) {
//             return kmerObj.findMatches(kmers)
//                 .then(function (matches) {
//                     kmerObj.close();
//                     let bestMatch = matches[0];
//                     return bestMatch.get('template').should.equal('NC_017625') &&
//                            bestMatch.get('score').should.equal(2295) &&
//                            bestMatch.get('expected').should.equal(108) &&
//                            bestMatch.get('z').should.equal(211.00) &&
//                            bestMatch.get('probability').should.equal(5.03e-23) &&
//                            bestMatch.get('frac-q').should.equal(74.14) &&
//                            bestMatch.get('frac-d').should.equal(47.02) &&
//                            bestMatch.get('depth').should.equal(0.36) &&
//                            bestMatch.get('total-frac-q').should.equal(74.14) &&
//                            bestMatch.get('total-frac-d').should.equal(47.02) &&
//                            bestMatch.get('total-temp-cover').should.equal(0.36) &&
//                            bestMatch.get('kmers-template').should.equal(4881) &&
//                            bestMatch.get('species').should.equal('Escherichia coli DH1');
//                 })
//                 .catch(function(error) {
//                     console.log(error);
//                     return false;
//                 });
//         });
//     });
// });

describe('kmerFinderServer findMatches test_long', function () {
    this.timeout(500000000);

    it('Winner takes All match NC_017625 (Escherichia coli DH1)', function () {
        let kmerObj = new KmerFinderServer(
            './test_data/test_long.fastq',
            'ATGAC', 16, 1, 1, true, 'mongo', 'mongodb://localhost:27017/Kmers',
            'KmerMap', 'winner'
        );
        // let kmers = new Map([['A', 2]]);
        return kmerObj.findKmers().then(function (kmers) {
            return kmerObj.findMatches(kmers)
                .promise
                .then(function (matches) {
                    kmerObj.close();
                    let bestMatch = matches[0];
                    return bestMatch.get('template').should.equal('NC_017625') &&
                           bestMatch.get('score').should.equal(2295) &&
                           bestMatch.get('expected').should.equal(108) &&
                           bestMatch.get('z').should.equal(211.00) &&
                           bestMatch.get('probability').should.equal(5.03e-23) &&
                           bestMatch.get('frac-q').should.equal(74.14) &&
                           bestMatch.get('frac-d').should.equal(47.02) &&
                           bestMatch.get('depth').should.equal(0.36) &&
                           bestMatch.get('total-frac-q').should.equal(74.14) &&
                           bestMatch.get('total-frac-d').should.equal(47.02) &&
                           bestMatch.get('total-temp-cover').should.equal(0.36) &&
                           bestMatch.get('kmers-template').should.equal(4881) &&
                           bestMatch.get('species').should.equal('Escherichia coli DH1');
                })
                .catch(function(error) {
                    console.log(error);
                    kmerObj.close();
                    return false;
                });
        });
    });

    // it('Winner takes All match NC_017625 (Escherichia coli DH1)', function () {
    //     let kmerObj = new KmerFinderServer(
    //         './test_data/test_long.fastq',
    //         'ATGAC', 16, 1, 1, true, 'mongo', 'mongodb://localhost:27017/Kmers',
    //         'complete_genomes_4', 'winner'
    //     );
    //     let kmers = new Map([['A', 1]]);
    //     // let json = require('.././test_data/db_short_results.json');
    //     // return kmerObj.findKmers().then(function (kmers) {
    //         return kmerObj.findMatches(kmers)
    //             .then(function (matches) {
    //                 kmerObj.close();
    //                 let bestMatch = matches[0];
    //                 return bestMatch.get('template').should.equal('NC_017625') &&
    //                        bestMatch.get('score').should.equal(2295) &&
    //                        bestMatch.get('expected').should.equal(108) &&
    //                        bestMatch.get('z').should.equal(211.00) &&
    //                        bestMatch.get('probability').should.equal(5.03e-23) &&
    //                        bestMatch.get('frac-q').should.equal(74.14) &&
    //                        bestMatch.get('frac-d').should.equal(47.02) &&
    //                        bestMatch.get('depth').should.equal(0.36) &&
    //                        bestMatch.get('total-frac-q').should.equal(74.14) &&
    //                        bestMatch.get('total-frac-d').should.equal(47.02) &&
    //                        bestMatch.get('total-temp-cover').should.equal(0.36) &&
    //                        bestMatch.get('kmers-template').should.equal(4881) &&
    //                        bestMatch.get('species').should.equal('Escherichia coli DH1');
    //             })
    //             .catch(function(error) {
    //                 console.log(error);
    //                 return false;
    //             });
    //     // });
    // });
    // it('Standard should match NC_017625 (Escherichia coli DH1)', function () {
    //     let kmerObj = new KmerFinderServer(
    //         './test_data/test_long.fastq',
    //         'ATGAC', 16, 1, 1, true, 'mongo', 'mongodb://localhost:27017/Kmers',
    //         'complete_genomes_4', 'standard'
    //     );
    //     return kmerObj.findKmers().then(function (kmers) {
    //         return kmerObj.findMatches(kmers)
    //             .then(function (matches) {
    //                 kmerObj.close();
    //                 let bestMatch = matches[0];
    //                 return bestMatch.get('template').should.equal('NC_017625') &&
    //                        bestMatch.get('score').should.equal(2295) &&
    //                        bestMatch.get('expected').should.equal(108) &&
    //                        bestMatch.get('z').should.equal(211.00) &&
    //                        bestMatch.get('probability').should.equal(5.03e-23) &&
    //                        bestMatch.get('frac-q').should.equal(74.14) &&
    //                        bestMatch.get('frac-d').should.equal(47.02) &&
    //                        bestMatch.get('depth').should.equal(0.36) &&
    //                        bestMatch.get('total-frac-q').should.equal(74.14) &&
    //                        bestMatch.get('total-frac-d').should.equal(47.02) &&
    //                        bestMatch.get('total-temp-cover').should.equal(0.36) &&
    //                        bestMatch.get('kmers-template').should.equal(4881) &&
    //                        bestMatch.get('species').should.equal('Escherichia coli DH1');
    //             })
    //             .catch(function(error) {
    //                 console.log(error);
    //                 return false;
    //             });
    //     });
    // });
});
//
// TODO: FIX FASTA parser
// describe('kmerFinderServer findMatches (4_20_6_2012_1430_975_193318_contigs)', function () {
//     this.timeout(500000000);
//     it('Standard should match NC_008463 (Pseudomonas aeruginosa UCBPP-PA14)', function () {
//         let kmerObj = new KmerFinderServer(
//             './test_data/4_20_6_2012_1430_975_193318_contigs.fsa',
//             'ATGAC', 16, 1, 1, true, 'mongo', 'mongodb://localhost:27017/Kmers',
//             'complete_genomes_4', 'standard'
//         );
//         return kmerObj.findKmers().then(function (kmers) {
//             return kmerObj.findMatches(kmers)
//                 .then(function (matches) {
//                     kmerObj.close();
//                     let match = matches[0];
//                     return match.get('template').should.equal('NC_008463') &&
//                            match.get('score').should.equal(3502) &&
//                            match.get('expected').should.equal(40) &&
//                            match.get('z').should.equal(537.20) &&
//                            match.get('probability').should.equal(5.03e-23) &&
//                            match.get('frac-q').should.equal(47.35) &&
//                            match.get('frac-d').should.equal(97.28) &&
//                            match.get('coverage').should.equal(0.49) &&
//                            match.get('ulength').should.equal(3600) &&
//                            match.get('species').should.equal('Pseudomonas aeruginosa UCBPP-PA14');
//                 });
//         });
//     });
//     it('Winner should match NC_008463 (Pseudomonas aeruginosa UCBPP-PA14)', function () {
//         let kmerObj = new KmerFinderServer(
//             './test_data/4_20_6_2012_1430_975_193318_contigs.fsa',
//             'ATGAC', 16, 1, 1, true, 'mongo', 'mongodb://localhost:27017/Kmers',
//             'complete_genomes_4', 'winner'
//         );
//         return kmerObj.findKmers().then(function (kmers) {
//             return kmerObj.findMatches(kmers)
//                 .then(function (matches) {
//                     kmerObj.close();
//                     let match = matches[0];
//                     return match.get('template').should.equal('NC_008463') &&
//                            match.get('score').should.equal(3502) &&
//                            match.get('expected').should.equal(40) &&
//                            match.get('z').should.equal(537.20) &&
//                            match.get('probability').should.equal(5.03e-23) &&
//                            match.get('frac-q').should.equal(47.35) &&
//                            match.get('frac-d').should.equal(97.28) &&
//                            match.get('coverage').should.equal(0.49) &&
//                            match.get('ulength').should.equal(3600) &&
//                            match.get('species').should.equal('Pseudomonas aeruginosa UCBPP-PA14');
//                 });
//         });
//     });
// });

// describe('kmerFinderServer Long', function () {
//     this.timeout(500000000);
//     let kmerObj = new KmerFinderServer('/Users/cisneror/code/genomic-git/kmerjs/./test_data/GMI15-006-DNA_S6_L001_R1_001.fastq');
//     it('long fastq should contain 6191 kmers!', function () {
//         // let json = require('.././test_data/kmers_long.json');
//         return kmerObj.findKmers().then(function (kmers) {
//             Console.log('answer size', [...data].length, 'U Kmers: ', kmers.uKmers);
//             return true;
//             // return  json.should.jsonEqual(mapToJSON(data));
//         });
//     });
// });
