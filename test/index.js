import {
    KmerJS, complement, kmers
} from '../lib/index';
import chai from 'chai';
import chaiAsPromised from 'chai-as-promised';
import Console from 'console';
chai.use(chaiAsPromised);
chai.should();

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

// describe('find kmers', function () {
//     this.timeout(50000);
//     it('Fastq linw should contain ATGACGCAATACTCCT', function () {
//         let kmerMap = new Map();
//         let seq = 'NTTTATGACGCAATACTCCTCTCTCCTTCGTGGTCTTGCAGCGGGTTCTGCATTTTTATTCCTTTTTGCCCCAACGGCATTCGCGGCGGAACAAACCGTTG';
//         let kmer = 'ATGACGCAATACTCCT';
//         kmers(seq, kmerMap, 16, 'ATGAC', 1);
//         let expected = new Map([[kmer, 1]]);
//         Console.log(kmerMap, expected);
//         kmerMap.should.deep.equal(expected);
//     });
// });
//
describe('Reverse complement', function () {
    it('ATGACCTGAGAGCCTT should be AAGGCTCTCAGGTCAT', function () {
        let seq = 'ATGACCTGAGAGCCTT';
        let rev = 'AAGGCTCTCAGGTCAT';
        let comp = complement(seq);
        comp.should.equal(rev);
    });
});


describe('kmerjs', function () {
    this.timeout(50000000);
    describe('short kmers', function () {
        let kmerjs = new KmerJS(
            './test_data/test_short.fastq',
            'ATGAC', 16, 1, 1, 'mongo');
        // let answer = kmerjs.findKmers();
        it('short fastq should contain 2 kmers!', function () {
            let kmerMap = new Map(
                [
                    ['ATGACGCAATACTCCT', 1],
                    ['ATGACCTGAGAGCCTT', 1]
                ]);
            return kmerjs.findKmers().should.eventually.deep.equal(kmerMap);
        });
        //
        // it('158 hits found in short fastq (MongoDB)', function () {
        //     let json = require('../test_data/db_short_results.json');
        //     let kmerMap = jsonToStrMap(json.templateentries);
        //     let expected = {
        //         templateentries: kmerMap,
        //         hits: json.hits
        //     };
        //     answer.then(function (kmers){
        //         return kmerjs.findMatches(kmers).should.eventually.be.deep.equal(expected);
        //     });
        // });
        //
        // kmerjs = new KmerJS(
        //     './test_data/test_short.fastq',
        //     'ATGAC', 16, 1, 1, '../out.log', 'json');
        // answer = kmerjs.findKmers();
        //
        // it('158 hits found in short fastq (JSON)', function () {
        //     let json = require('../test_data/db_short_results.json');
        //     answer.then(function (kmers){
        //         kmerjs.findMatches(kmers).should.have.property('hits', 158);
        //         chai.assert.deepEqual(
        //             answer.templateentries, jsonToStrMap(json.templateentries));
        //         chai.assert.deepEqual(
        //             answer.templateentriestot, jsonToStrMap(json.templateentriestot)
        //         );
        //     });
        // });
    });
    // describe('long kmers', function () {
    //     let kmerjs = new KmerJS(
    //         './test_data/test_long.fastq',
    //         'ATGAC', 16, 1, 1, '../out.log', 'mongo');
    //     let kmers = kmerjs.findKmers();
    //
    //     it('long fastq should contain 3766 kmers!', function () {
    //         let json = require('../test_data/kmers_long.json');
    //         let kmerMap = jsonToStrMap(json);
    //         return kmers.should.eventually.deep.equal(kmerMap);
    //     });
    //
    //     it('179108 hits found in long fastq (MongoDB)', function () {
    //         let json = require('../test_data/db_long_results.json');
    //         let kmerMap = jsonToStrMap(json.templateentries);
    //         let expected = {
    //             templateentries: kmerMap,
    //             hits: json.hits
    //         };
    //         kmers.then(function (answer) {
    //             return answer.findMatches(kmers)
    //                     .should.eventually.be.deep.equal(expected);
    //         });
    //     });
    //
    //     kmerjs = new KmerJS(
    //         './test_data/test_long.fastq',
    //         'ATGAC', 16, 1, 1, '../out.log', 'json');
    //     kmers = kmerjs.findKmers();
    //
    //     it('179108 hits found in long fastq (JSON)', function () {
    //         let json = require('../test_data/db_long_results.json');
    //         kmers.then(function (answer){
    //             let out = kmerjs.findMatches(answer);
    //             out.should.have.property('hits', 179108);
    //             chai.assert.deepEqual(
    //                 out.templateentries, jsonToStrMap(json.templateentries));
    //             chai.assert.deepEqual(
    //                 out.templateentriestot, jsonToStrMap(json.templateentriestot)
    //             );
    //         });
    //     });
    // });
});
