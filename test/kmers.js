import {
    KmerJS,
    complement
} from '../lib/kmers';
import chai from 'chai';
import chaiAsPromised from 'chai-as-promised';

chai.use(chaiAsPromised);
chai.should();

describe('kmerjs', function () {
    it('Fastq line should contain ATGACGCAATACTCCT', function () {
        let kmers = new KmerJS();
        let seq = `NTTTATGACGCAATACTCCTCTCTCCTTCGTGGTCTTGCAGCGGGTTCTGC
                   ATTTTTATTCCTTTTTGCCCCAACGGCATTCGCGGCGGAACAAACCGTTG`;
        let kmer = 'ATGACGCAATACTCCT';
        kmers.kmersInLine(seq);
        return [...kmers.kmerMap][0][0].should.equal(kmer);
    });

    it('ATGACCTGAGAGCCTT should be AAGGCTCTCAGGTCAT', function () {
        let seq = 'ATGACCTGAGAGCCTT';
        let rev = 'AAGGCTCTCAGGTCAT';
        let comp = complement(seq);
        return comp.should.equal(rev);
    });
});
