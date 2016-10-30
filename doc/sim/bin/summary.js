#!/usr/bin/env node
var fs = require("fs");
var methods = process.argv[2].split(" ");
var dirs = process.argv.slice(3);
function mean_rmserr(n,m1,m2){m = m1/n;return [m.toPrecision(6),Math.sqrt(m2/n-2*m+1).toPrecision(6)]}
process.stdout.write(["method","insmean","insrmse","delmean","delrmse","insdelmean","insdelrmse"].join(" ")+"\n")
methods.forEach (function (method) {
    var n = 0, m1i = 0, m2i = 0, m1d = 0, m2d = 0;
    dirs.forEach (function (dir) {
        var match = /.*rate([\d\.]+).*/.exec(dir);
        var rate = match[1];
        var filenameRegex = new RegExp('\\d+\\.' + method + '\\.json');
        fs.readdirSync(dir)
            .filter (function (filename) { return filenameRegex.test(filename) })
            .map (function (filename) {
                try {
                    var json = JSON.parse(fs.readFileSync(dir+'/'+filename));
                    var ir = json.insrate/rate;
                    var dr = json.delrate/rate;
                    m1i += ir;
                    m2i += ir*ir;
                    m1d += dr;
                    m2d += dr*dr;
                    ++n;
                } catch (err) {
                    process.stderr.write(err.toString()+" ("+dir+"/"+filename+")\n")
                }
            })
    })
    if (n > 0)
        process.stdout.write([method]
                             .concat(mean_rmserr(n,m1i,m2i))
                             .concat(mean_rmserr(n,m1d,m2d))
                             .concat(mean_rmserr(n*2,m1i+m1d,m2i+m2d))
                             .join(" ")+"\n")
})
