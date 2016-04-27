import {
    KmerFinderClient,
} from '../lib/kmerFinderClient';
import chai from 'chai';
import chaiAsPromised from 'chai-as-promised';

chai.use(chaiAsPromised);
chai.should();

describe('kmerFinderClient findKmers', function () {
    this.timeout(50000);
    let kmers = new KmerFinderClient('/Users/cisneror/code/genomic-git/kmerjs/test_data/test_short.fastq', 'node');
    let answer = kmers.findKmers();
    it('short fastq should contain 2 kmers!', function () {
        return answer.then(function (data) {
            let expected = [['ATGACGCAATACTCCT', 1],
                            ['ATGACCTGAGAGCCTT', 1]];
            return [...data].should.deep.equal(expected);
        });
    });
});

// describe('kmerFinderClient findMatches', function () {
//     this.timeout(500000000);
//     let kmers = new KmerFinderClient('/Users/cisneror/code/genomic-git/kmerjs/test_data/test_long.fastq',
//                          'ATGAC', 16, 1, 1, '', 'json', '/Users/cisneror/code/genomic-git/kmerjs/test_data/db.json');
//     let answer = kmers.findKmers();
//     it('179 hits found in long fastq (JSON)', function () {
//         return answer.then(function (data){
//             return kmers.findMatches(data).then(function(matches){
//                 return matches.should.have.length(179);
//             });
//         });
//     });
// });
