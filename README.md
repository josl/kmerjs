# kmerjs [![NPM version][npm-image]][npm-url] [![Build Status][travis-image]][travis-url] [![Dependency Status][daviddm-image]][daviddm-url] [![Coverage Status](https://coveralls.io/repos/josl/kmerjs/badge.svg?branch=master&service=github)](https://coveralls.io/github/josl/kmerjs?branch=master)
> Finds kmers in a fastq file

## Install

```sh
$ npm install --save kmerjs
```

## Usage

```js
import kmerjs from 'kmerjs';

let kmerMap = kmerjs('test_data/test_short.fastq', 'ATGAC', 16, 1, 'output');
```

## License
Apache-2.0 © [Jose Luis Bellod Cisneros](http://josl.github.io)

[npm-image]: https://badge.fury.io/js/kmerjs.svg
[npm-url]: https://npmjs.org/package/kmerjs
[travis-image]: https://travis-ci.org/josl/kmerjs.svg?branch=master
[travis-url]: https://travis-ci.org/josl/kmerjs
[daviddm-image]: https://david-dm.org/josl/kmerjs.svg?theme=shields.io
[daviddm-url]: https://david-dm.org/josl/kmerjs
[coveralls-image]: https://coveralls.io/repos/josl/kmerjs/badge.svg
[coveralls-url]: https://coveralls.io/r/josl/kmerjs
