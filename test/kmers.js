import {
    KmerJS,
    complement
} from '../lib/kmers';
import chai from 'chai';
import chaiAsPromised from 'chai-as-promised';

chai.use(chaiAsPromised);
chai.should();
//
describe('kmerjs', function () {
    it('ATGACGCAATACTCCT in kmersInLine()', function () {
        let kmers = new KmerJS();
        let seq = `NTTTATGACGCAATACTCCTCTCTCCTTCGTGGTCTTGCAGCGGGTTCTGC
                   ATTTTTATTCCTTTTTGCCCCAACGGCATTCGCGGCGGAACAAACCGTTG`;
        let kmer = 'ATGACGCAATACTCCT';
        kmers.kmersInLine(seq);
        return [...kmers.kmerMap][0][0].should.equal(kmer);
    });

    it('complement(ATGACCTGAGAGCCTT) = AAGGCTCTCAGGTCAT', function () {
        let seq = 'ATGACCTGAGAGCCTT';
        let rev = 'AAGGCTCTCAGGTCAT';
        let comp = complement(seq);
        return comp.should.equal(rev);
    });

    it('readFile(test_short.fastq) should have 2 kmer', function () {
        this.timeout(500000000);
        let file = './test_data/test_short.fastq';
        let kmers = new KmerJS(file);
        return kmers.readFile().then(function (kmers) {
            return kmers.size.should.equal(2);
        });
    });

    it('readFile(test_long.fastq) should have 6191 kmers', function () {
        this.timeout(500000000);
        let file = './test_data/test_long.fastq';
        let kmers = new KmerJS(file);
        return kmers.readFile().then(function (kmers) {
            return kmers.size.should.equal(6191);
        });
    });
    it('readFile(test_long.kmer.fastq) should have 401 kmers', function () {
        this.timeout(500000000);
        let file = './test_data/test_long.kmer.fastq';
        let kmers = new KmerJS(file);
        return kmers.readFile().then(function (kmers) {
            return kmers.size.should.equal(401);
        });
    });
    // TODO: FASTA tests missing!
    // it('readFile(4_20_6_2012_1430_975_193318_contigs.fsa) should have 7196 kmers', function () {
    //     this.timeout(500000000);
    //     let file = './test_data/4_20_6_2012_1430_975_193318_contigs.fsa';
    //     let kmers = new KmerJS(file);
    //     return kmers.readFile().then(function (kmers) {
    //         return kmers.size.should.equal(7196);
    //     });
    // });
});
