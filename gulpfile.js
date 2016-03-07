var path = require('path');
var gulp = require('gulp');
var eslint = require('gulp-eslint');
var excludeGitignore = require('gulp-exclude-gitignore');
var mocha = require('gulp-mocha');
var istanbul = require('gulp-istanbul');
var nsp = require('gulp-nsp');
var plumber = require('gulp-plumber');
var coveralls = require('gulp-coveralls');
var babel = require('gulp-babel');
var del = require('del');
var isparta = require('isparta');
var browserify = require('gulp-browserify');
var babelify = require("babelify");
var fs = require('fs');

// Initialize the babel transpiler so ES2015 files gets compiled
// when they're loaded
require('babel-core/register');

gulp.task('static', function () {
    return gulp.src('**/*.js')
        .pipe(excludeGitignore())
        .pipe(eslint())
        .pipe(eslint.format())
        .pipe(eslint.failAfterError());
});

gulp.task('nsp', function (cb) {
    nsp({
        package: path.resolve('package.json')
    }, cb);
});

gulp.task('pre-test', function () {
    return gulp.src('lib/**/*.js')
        .pipe(istanbul({
            includeUntested: true,
            instrumenter: isparta.Instrumenter
        }))
        .pipe(istanbul.hookRequire());
});

gulp.task('test', ['pre-test'], function (cb) {
    var mochaErr;

    gulp.src('test/**/*.js')
        .pipe(plumber())
        .pipe(mocha({
            reporter: 'spec'
        }))
        .on('error', function (err) {
            mochaErr = err;
        })
        .pipe(istanbul.writeReports())
        .on('end', function () {
            cb(mochaErr);
        });
});

gulp.task('coveralls', ['test'], function () {
    if (!process.env.CI) {
        return;
    }

    return gulp.src(path.join(__dirname, 'coverage/lcov.info'))
        .pipe(coveralls());
});

gulp.task('babel', ['clean'], function () {
    return gulp.src('lib/**/*.js')
        .pipe(babel())
        .pipe(gulp.dest('dist'));
});

gulp.task('clean', function () {
    return del('dist');
});

// Browserify
gulp.task('browser', function () {
    // Single entry point to browserify
    //  browserify dist/index.js  -o uploader/app/index.js -t [ babelify --presets [ es2015] ] -s kmerModule
    //  browserify dist/stats.js -o uploader/app/stats.js -t [ babelify --presets [ es2015] ] -s statModule
    return gulp.src('dist/**/*.js')
        .pipe(browserify({ debug: true })
              .transform(babelify)
              .bundle()
            //   .on('error', function (err) { console.log('Error: ' + err.message); })
              .pipe(fs.createWriteStream('bundle.js')))
        .pipe(gulp.dest('./browser'));
});

gulp.task('prepublish', ['nsp', 'babel', 'browserify']);
gulp.task('default', ['static', 'test', 'coveralls']);
