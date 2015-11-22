import BN from 'bignumber.js';
import Console from 'console';
export let etta = new BN(1.0e-8);

/**
* Conservative two sided p-value from z-score
* Comparison of two fractions,
* Statistical methods in medical research, Armitage et al. p. 125
 * @param  {[int]} r1 [positives in sample 1]
 * @param  {[int]} n1 [size of sample 1]
 * @param  {[int]} r2 [positives in sample 2]
 * @param  {[int]} n2 [size of sample 2]
 * @return {[float]}    [z-score]
 */

export function zScore(r1, n1, r2, n2){
    // Console.log(r1, n1, r2, n2);
    let p1 = new BN(r1).dividedBy(n1).plus(etta);
    // Console.log(p1);
    // let p2 = r2 / n2 + etta;
    let p2 = new BN(r2).dividedBy(n2).plus(etta);
    // Console.log(p2);
    // let q1 = 1 - p1;
    // let q2 = 1 - p2;
    // let p = (r1 + r2) / (n1 + n2 + etta);

    let p = new BN(r1).plus(r2).dividedBy(new BN(n1).plus(n2).plus(etta));
    let q = new BN(1).minus(p);
    // Console.log(p);
    // Console.log(q);
    // let square = Math.sqrt(p * q * (1 / (n1 + etta) + 1 / (n2 + etta)) + etta);
    // sqrt(p*q*(1/(n1+etta)+1/(n2+etta))+etta)
    let square = new BN(new BN(p).times(q).times(
        new BN(1).dividedBy(new BN(n1).plus(etta))
        .plus(new BN(1).dividedBy(new BN(n2).plus(etta))))
    ).plus(etta).sqrt();
    // Console.log(square);
    // let z = (p1 - p2) / square;
    let z = new BN(p1).minus(p2).dividedBy(square);
    // Console.log(r1, n1, r2, n2, p.toNumber(), q.toNumber(), square.toNumber(), z.toNumber());
    return z;
}

/**
 * [Conservative two sided p-value from z-score]
 * @param  {[float]} z [z-score]
 * @return {[float]}   [p-value]
 */
export function fastp(z){
    // Console.log(z, z.toNumber(), new BN(10.7016));
    // Console.log(z.comparedTo(new BN(10.7016)));
    let p = 0.0;
    if (z.comparedTo(new BN(10.7016)) > 0){
        p = 1e-26;
    }else if (z.comparedTo(new BN(10.4862)) > 0){
        p = 1e-25;
    }else if (z.comparedTo(new BN(10.2663)) > 0){
        p = 1e-24;
    }else if (z.comparedTo(new BN(10.0416)) > 0){
        p = 1e-23;
    }else if (z.comparedTo(new BN(9.81197)) > 0){
        p = 1e-22;
    }else if (z.comparedTo(new BN(9.5769)) > 0){
        p = 1e-21;
    }else if (z.comparedTo(new BN(9.33604)) > 0){
        p = 1e-20;
    }else if (z.comparedTo(new BN(9.08895)) > 0){
        p = 1e-19;
    }else if (z.comparedTo(new BN(8.83511)) > 0){
        p = 1e-18;
    }else if (z.comparedTo(new BN(8.57394)) > 0){
        p = 1e-17;
    }else if (z.comparedTo(new BN(8.30479)) > 0){
        p = 1e-16;
    }else if (z.comparedTo(new BN(8.02686)) > 0){
        p = 1e-15;
    }else if (z.comparedTo(new BN(7.73926)) > 0){
        p = 1e-14;
    }else if (z.comparedTo(new BN(7.4409)) > 0){
        p = 1e-13;
    }else if (z.comparedTo(new BN(7.13051)) > 0){
        p = 1e-12;
    }else if (z.comparedTo(new BN(6.8065)) > 0){
        p = 1e-11;
    }else if (z.comparedTo(new BN(6.46695)) > 0){
        p = 1e-10;
    }else if (z.comparedTo(new BN(6.10941)) > 0){
        p = 1e-9;
    }else if (z.comparedTo(new BN(5.73073)) > 0){
        p = 1e-8;
    }else if (z.comparedTo(new BN(5.32672)) > 0){
        p = 1e-7;
    }else if (z.comparedTo(new BN(4.89164)) > 0){
        p = 1e-6;
    }else if (z.comparedTo(new BN(4.41717)) > 0){
        p = 1e-5;
    }else if (z.comparedTo(new BN(3.89059)) > 0){
        p = 1e-4;
    }else if (z.comparedTo(new BN(3.29053)) > 0){
        p = 1e-3;
    }else if (z.comparedTo(new BN(2.57583)) > 0){
        p = 0.01;
    }else if (z.comparedTo(new BN(1.95996)) > 0){
        p = 0.05;
    }else if (z.comparedTo(new BN(1.64485)) > 0){
        p = 0.1;
    }else {
        p = 1.0;
    }
    // Console.log(p, new BN(p));
    return new BN(p);
}
