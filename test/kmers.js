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

    it('test_short.fastq should have 2 kmer', function () {
        this.timeout(500000000);
        let file = '/Users/cisneror/code/genomic-git/kmerjs/test_data/test_short.fastq';
        let kmers = new KmerJS(file);
        return kmers.readFile().then(function (kmers) {

            return kmers.size.should.equal(2);
        });
    });

    it('test_long.fastq should have 6045 kmers', function () {
        this.timeout(500000000);
        let file = '/Users/cisneror/code/genomic-git/kmerjs/test_data/test_long.fastq';
        let kmers = new KmerJS(file);
        return kmers.readFile().then(function (kmers) {
            console.log(kmers);
            return kmers.size.should.equal(6045);
        });
    });

    it('4_20_6_2012_1430_975_193318_contigs.fsa should have 7196 kmers', function () {
        this.timeout(500000000);
        let file = '/Users/cisneror/Dropbox/PhD/CGE/ringtrials/4_20_6_2012_1430_975_193318_contigs.fsa';
        let kmers = new KmerJS(file);
        return kmers.readFile().then(function (kmers) {
            return kmers.size.should.equal(7196);
        });
    });
});
