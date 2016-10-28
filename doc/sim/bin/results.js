#!/usr/bin/env node

var fs = require("fs");

var binsize = process.argv[2],
    binmax = process.argv[3],
    nbins = (binmax/binsize)+1,
    togetherFlag = eval(process.argv[4]),
    methods = process.argv[5].split(" "),
    dirs = process.argv.slice(6)

function incbin(r,bins) {
    var bin = Math.round(r/binsize)
    bins[bin] = (bins[bin] || 0) + 1
}

var binStart = Array.apply(null,Array(nbins)).map (function(x,n) { return (n*binsize).toPrecision(2) })
function printbins(label,headers,bins) {
    binStart.forEach (function (start, n) {
        process.stdout.write ([label]
                              .concat(headers)
                              .concat([start,bins[n] || 0,"\n"])
                              .join(" "))
    })
}

process.stdout.write (["event","method"]
                      .concat (togetherFlag ? [] : ["rate","tree","multiplier"])
                      .concat(["bin","count","\n"])
                      .join (" "))

methods.forEach (function(method) {
    var ibins = {}, dbins = {}
    dirs.forEach (function(dir) {
        var match = /.*rate([\d\.]+)\/([a-z]+)(\d+)/.exec(dir),
            rate = match[1],
            tree = match[2],
            mul = match[3],
            filenameRegex = new RegExp('\\d+\\.' + method + '\\.json'),
            headers
        fs.readdirSync(dir)
            .filter (function (filename) { return filenameRegex.test(filename) })
            .map (function (filename) {
                try {
                    var json = JSON.parse(fs.readFileSync(dir+"/"+filename));
                    incbin(json.insrate/rate,ibins);
                    incbin(json.delrate/rate,dbins)
                } catch(err) {
                    process.stderr.write(err.toString()+" ("+dir+"/"+filename+")\n")
                }
            })
        if (!togetherFlag) {
            printbins("ins",[method,rate,tree,mul],ibins)
            printbins("del",[method,rate,tree,mul],dbins)
            ibins = {}
            dbins = {}
        }
    })
    if (togetherFlag) {
        printbins("ins",[method],ibins)
        printbins("del",[method],dbins)
    }
})
