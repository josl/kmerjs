import {
    KmerFinderServer,
} from '../lib/kmerFinderServer';
import chai from 'chai';
import chaiAsPromised from 'chai-as-promised';
import chaiJsonEqual from 'chai-json-equal';
import Console from 'console';

chai.use(chaiJsonEqual);
chai.use(chaiAsPromised);
chai.should();

// describe('kmerFinderServer findKmers', function () {
//     this.timeout(500000000);
//     let kmers = new KmerFinderServer('/Users/cisneror/code/genomic-git/kmerjs/test_data/test_short.fastq');
//     let answer = kmers.findKmers();
//     it('short fastq should contain 2 kmers!', function () {
//         return answer.then(function (data) {
//             kmers.db.connection.disconnect();
//             let expected = [
//                 ['ATGACGCAATACTCCT', 1],
//                 ['ATGACCTGAGAGCCTT', 1]
//             ];
//             return [...data].should.deep.equal(expected);
//         });
//     });
// });

// describe('kmerFinderServer findMatches', function () {
//     this.timeout(500000000);
//     let kmers = new KmerFinderServer(
//         '/Users/cisneror/code/genomic-git/kmerjs/test_data/test_long.fastq',
//         'ATGAC', 16, 1, 1, '', 'mongo', 'mongodb://localhost:27017/Kmers',
//         'complete_genomes_2'
//     );
//     let answer = kmers.findKmers();
//     it('Long fastq should matche NC_010554 (Proteus mirabilis HI4320)', function () {
//         // let json = require('../test_data/db_short_results.json');
//         return answer.then(function (data) {
//             return kmers.findMatches(data)
//                 .then(function (matches) {
//                     kmers.db.connection.disconnect();
//                     return matches[0].get('template').should.equal('NC_010554');
//                 });
//         });
//     });
// });

// describe('kmerFinderServer findMatches', function () {
//     this.timeout(500000000);
//     let kmers = new KmerFinderServer(
//         '/Users/cisneror/Dropbox/PhD/CGE/ringtrials/4_20_6_2012_1430_975_193318_contigs.fsa',
//         'ATGAC', 16, 1, 1, 'test', 'mongo', 'mongodb://localhost:27017/Kmers',
//         'complete_genomes_2'
//     );
//     let answer = kmers.findKmers();
//     it('4_20_6_2012_1430_975_193318_contigs should match NC_008463 (Pseudomonas aeruginosa UCBPP-PA14)', function () {
//         // let json = require('../test_data/db_short_results.json');
//         return answer.then(function (data) {
//             return kmers.findMatches(data)
//                 .then(function (matches) {
//                     kmers.db.connection.disconnect();
//                     let match = matches[0];
//                     // Map {110481 / Kmers: 1375
//                     //     'template' => 'NC_008463',
//                     //     'score' => 675,
//                     //     'expected' => 2.112850798361426,
//                     //     'z' => 422.0645825021683,
//                     //     'probability' => 5.03e-23,
//                     //     'frac-q' => 13500000000000,
//                     //     'frac-d' => 18.749999999947917,
//                     //     'coverage' => 0.09251161519540858,
//                     //     'ulength' => 3600,
//                     //     'species' => 'Pseudomonas aeruginosa UCBPP-PA14' }
//                     return match.get('template').should.equal('NC_008463') &&
//                            match.get('score').should.equal(675) &&
//                            match.get('expected').should.equal(2.112850798361426) &&
//                            match.get('z').should.equal(422.0645825021683) &&
//                            match.get('probability').should.equal(5.03e-23) &&
//                            match.get('frac-q').should.equal(13500000000000) &&
//                            match.get('frac-d').should.equal(18.749999999947917) &&
//                            match.get('coverage').should.equal(0.09251161519540858) &&
//                            match.get('ulength').should.equal(3600) &&
//                            match.get('species').should.equal('Pseudomonas aeruginosa UCBPP-PA14');
//                 });
//         });
//     });
// });
//
describe('kmerFinderServer findMatches', function () {
    this.timeout(500000000);
    // fastq,
    // preffix = 'ATGAC',
    // length = 16,
    // step = 1,
    // coverage = 1,
    // progress = true,
    // db = 'mongo',
    // url='mongodb:\/\/localhost:27017/Kmers',
    // collection='genomes',
    // method = 'standard'
    let kmers = new KmerFinderServer(
        '/Users/cisneror/Dropbox/PhD/CGE/ringtrials/4_20_6_2012_1430_975_193318_contigs.fsa',
        'ATGAC', 16, 1, 1, true, 'mongo', 'mongodb://localhost:27017/Kmers',
        'complete_genomes_4', 'standard'
    );
    let answer = kmers.findKmers();
    it('4_20_6_2012_1430_975_193318_contigs should match NC_008463 (Pseudomonas aeruginosa UCBPP-PA14)', function () {
        // let json = require('../test_data/db_short_results.json');
        return answer.then(function (data) {
            return kmers.findMatches(data)
                .then(function (matches) {
                    kmers.db.connection.disconnect();
                    let match = matches[0];
                    // Map {'template' => 'NC_008463',
                    //     'score' => 3502,
                    //     'expected' => 40,
                    //     'z' => 537.20,
                    //     'probability' => 5.03e-23,
                    //     'frac-q' => 47.35,
                    //     'frac-d' => 97.28,
                    //     'coverage' => 0.49,
                    //     'ulength' => 3600,
                    //     'species' => 'Pseudomonas aeruginosa UCBPP-PA14' }
                    return match.get('template').should.equal('NC_008463') &&
                           match.get('score').should.equal(3502) &&
                           match.get('expected').should.equal(40) &&
                           match.get('z').should.equal(537.20) &&
                           match.get('probability').should.equal(5.03e-23) &&
                           match.get('frac-q').should.equal(47.35) &&
                           match.get('frac-d').should.equal(97.28) &&
                           match.get('coverage').should.equal(0.49) &&
                           match.get('ulength').should.equal(3600) &&
                           match.get('species').should.equal('Pseudomonas aeruginosa UCBPP-PA14');
                });
        });
    });
});

// describe('kmerFinderServer Long', function () {
//     this.timeout(500000000);
//     let kmers = new KmerFinderServer('/Users/cisneror/code/genomic-git/kmerjs/test_data/GMI15-006-DNA_S6_L001_R1_001.fastq');
//     it('long fastq should contain 6191 kmers!', function () {
//         // let json = require('../test_data/kmers_long.json');
//         return kmers.findKmers().then(function (data) {
//             Console.log('answer size', [...data].length, 'U Kmers: ', kmers.uKmers);
//             return true;
//             // return  json.should.jsonEqual(mapToJSON(data));
//         });
//     });
// });
