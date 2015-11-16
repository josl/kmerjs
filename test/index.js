import {
    kmerjs, complement, findMatches, findMatchesJSON
} from '../lib/index';
import chai from 'chai';
import chaiAsPromised from 'chai-as-promised';

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

describe('complement', function () {
    this.timeout(50000);
    it('ATGACCTGAGAGCCTT should be TAGTCCAGTGTGCCAA', function () {
        let seq = 'ATGACCTGAGAGCCTT';
        let rev = 'TAGTCCAGTGTGCCAA';
        return complement(seq)
            .should.equal(rev);
    });
});

describe('kmerjs', function () {
    this.timeout(50000000);
    it('short fastq should contain 2 kmers!', function () {
        let kmers = new Map(
            [
                ['ATGACGCAATACTCCT', 1],
                ['ATGACCTGAGAGCCTT', 1]
            ]);
        return kmerjs('./test_data/test_short.fastq', 'ATGAC', 16, 1, '../out.log')
            .should.eventually.deep.equal(kmers);
    });

    it('long fastq should contain 3766 kmers!', function () {
        let json = require('../test_data/kmers_long.json');
        let kmerMap = jsonToStrMap(json);
        return kmerjs('./test_data/test_long.fastq', 'ATGAC', 16, 1, '../out.log')
            .should.eventually.deep.equal(kmerMap);
    });

    it('158 hits found in short fastq (MongoDB)', function () {
        let kmers = new Map(
            [
                ['ATGACGCAATACTCCT', 1],
                ['ATGACCTGAGAGCCTT', 1]
            ]);
        let json = require('../test_data/db_short_results.json');
        let kmerMap = jsonToStrMap(json.templateentries);
        let expected = {
            templateentries: kmerMap,
            hits: json.hits
        };
        return findMatches(kmers, 1).should.eventually.be.deep.equal(expected);
    });
    it('179108 hits found in long fastq (MongoDB)', function () {
        let json = require('../test_data/kmers_long.json');
        let kmers = jsonToStrMap(json);
        json = require('../test_data/db_long_results.json');
        let kmerMap = jsonToStrMap(json.templateentries);
        let expected = {
            templateentries: kmerMap,
            hits: json.hits
        };
        return findMatches(kmers, 1).should.eventually.be.deep.equal(expected);
    });

    it('158 hits found in short fastq (JSON)', function () {
        let kmers = new Map(
            [
                ['ATGACGCAATACTCCT', 1],
                ['ATGACCTGAGAGCCTT', 1]
            ]);
        let json = require('../test_data/db_short_results.json');
        let answer = findMatchesJSON(kmers, 1);
        chai.assert.deepEqual(answer.templateentries, jsonToStrMap(json.templateentries));
        chai.assert.deepEqual(
            answer.templateentriestot, jsonToStrMap(json.templateentriestot)
        );
        answer.should.have.property('hits', 158);
    });
    it('179108 hits found in long fastq (JSON)', function () {
        let json = require('../test_data/kmers_long.json');
        let kmerMap = jsonToStrMap(json);
        let answer = findMatchesJSON(kmerMap, 1);
        json = require('../test_data/db_long_results.json');
        chai.assert.deepEqual(answer.templateentries, jsonToStrMap(json.templateentries));
        chai.assert.deepEqual(
            answer.templateentriestot, jsonToStrMap(json.templateentriestot)
        );
        answer.should.have.property('hits', 179108);
    });

});
