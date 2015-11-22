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
//     it('Fastq line should contain ATGACGCAATACTCCT', function () {
//         let kmerMap = new Map();
//         let seq = `NTTTATGACGCAATACTCCTCTCTCCTTCGTGGTCTTGCAGCGGGTTCTGC
//                    ATTTTTATTCCTTTTTGCCCCAACGGCATTCGCGGCGGAACAAACCGTTG`;
//         let kmer = 'ATGACGCAATACTCCT';
//         kmers(seq, kmerMap, 16, 'ATGAC', 1);
//         let expected = new Map([[kmer, 1]]);
//         kmerMap.should.deep.equal(expected);
//     });
// });
//
// describe('Reverse complement', function () {
//     it('ATGACCTGAGAGCCTT should be AAGGCTCTCAGGTCAT', function () {
//         let seq = 'ATGACCTGAGAGCCTT';
//         let rev = 'AAGGCTCTCAGGTCAT';
//         let comp = complement(seq);
//         comp.should.equal(rev);
//     });
// });

describe('kmerjs', function () {
    this.timeout(50000000);
    describe('small fastq', function () {
        let kmerjs = new KmerJS(
            './test_data/test_short.fastq',
            'ATGAC', 16, 1, 1, '', 'mongo');
        let answer = kmerjs.findKmers();





        // it('4 entries after stat (MongoDB)', function () {
        //     let json = require('../test_data/db_short_results.json');
        //     let expected = {
        //         templateentries: jsonToStrMap(json.templateentries),
        //         templateentriestot: jsonToStrMap(json.templateentriestot),
        //         hits: json.hits
        //     };
        //     let kmerMap = new Map(
        //         [
        //             ['ATGACGCAATACTCCT', 1],
        //             ['ATGACCTGAGAGCCTT', 1]
        //         ]);
        //     kmerjs.uKmers = 0;
        //     // let answer = kmerjs.findMatches(kmerMap);
        //     // Console.log(answer);
        //     // answer.then(function (data) {
        //     //     Console.log(data);
        //     //     return true;
        //     // });
        //     return kmerjs.findMatches(kmerMap)
        //                 .should.eventually.be.empty;
        // });





        // it('short fastq should contain 2 kmers!', function () {
        //     let kmerMap = new Map(
        //         [
        //             ['ATGACGCAATACTCCT', 1],
        //             ['ATGACCTGAGAGCCTT', 1]
        //         ]);
        //     return answer.should.eventually.deep.equal(kmerMap);
        // });
        // it('158 hits found in short fastq (MongoDB)', function () {
        //     let json = require('../test_data/db_short_results.json');
        //     let expected = {
        //         templateentries: jsonToStrMap(json.templateentries),
        //         templateentriestot: jsonToStrMap(json.templateentriestot),
        //         hits: json.hits
        //     };
        //     let kmerMap = new Map(
        //         [
        //             ['ATGACGCAATACTCCT', 1],
        //             ['ATGACCTGAGAGCCTT', 1]
        //         ]);
        //     // let answer = kmerjs.findMatches(kmerMap);
        //     // Console.log(answer);
        //     // answer.then(function (data) {
        //     //     Console.log(data);
        //     //     return true;
        //     // });
        //     return kmerjs.findMatches(kmerMap)
        //                 .should.eventually.have
        //                 .property('hits', expected.hits);
        // });
        // kmerjs = new KmerJS(
        //     './test_data/test_short.fastq',
        //     'ATGAC', 16, 1, 1, '', 'json');
        // answer = kmerjs.findKmers();
        //
        // it('158 hits found in short fastq (JSON)', function () {
        //     let json = require('../test_data/db_short_results.json');
        //     answer.then(function (kmersMap){
        //         kmerjs.findMatches(kmersMap).should.have.property('hits', 158);
        //         chai.assert.deepEqual(
        //             kmersMap.templateentries, jsonToStrMap(json.templateentries));
        //     });
        // });
    });
    describe('big fastq', function () {
        let kmerjs = new KmerJS(
            './test_data/test_long.fastq',
            'ATGAC', 16, 1, 1, '', 'mongo');
        let answer = kmerjs.findKmers();
        it('kmer hit results', function () {
            return answer.then(function (kmerMap) {
                return kmerjs.findMatches(kmerMap);
            }).should.eventually.have.length(179);
            // return answer.then()
            // return answer.should.eventually.deep.equal(kmerMap);
        });

    //
    //     it('long fastq should contain 3766 kmers!', function () {
    //         let json = require('../test_data/kmers_long.json');
    //         let kmerMap = jsonToStrMap(json);
    //         return answer.should.eventually.deep.equal(kmerMap);
    //     });
    //
    //     it('179108 hits found in long fastq (MongoDB)', function () {
    //         let json = require('../test_data/db_long_results.json');
    //         let kmerMap = jsonToStrMap(json.templateentries);
    //         let expected = {
    //             templateentries: kmerMap,
    //             hits: json.hits
    //         };
    //         answer.then(function (kmersMap) {
    //             return kmerjs.findMatches(kmersMap)
    //                     .should.eventually.be.deep.equal(expected);
    //         });
    //     });
    //
    //     kmerjs = new KmerJS(
    //         './test_data/test_long.fastq',
    //         'ATGAC', 16, 1, 1, '', 'json');
    //     answer = kmerjs.findKmers();
    //
    //     it('179108 hits found in long fastq (JSON)', function () {
    //         let json = require('../test_data/db_long_results.json');
    //         answer.then(function (kmersMap){
    //             let out = kmerjs.findMatches(kmersMap);
    //             out.should.have.property('hits', 179108);
    //             chai.assert.deepEqual(
    //                 out.templateentries, jsonToStrMap(json.templateentries));
    //             chai.assert.deepEqual(
    //                 out.templateentriestot, jsonToStrMap(json.templateentriestot)
    //             );
    //         });
    //     });
    });
});
