#!/usr/bin/env node

var fs = require("fs");

var binsize = process.argv[2],
    binmax = process.argv[3],
    nbins = (binmax/binsize)+1,
    quartileFlag = eval(process.argv[4]),
    methods = process.argv[5].split(" "),
    dirs = process.argv.slice(6)

function incbin(r,bins,vals) {
    var bin = Math.round(r/binsize)
    bins[bin] = (bins[bin] || 0) + 1
    vals.push(r)
}

process.stdout.write ((quartileFlag
                       ? ["method","tree","lower","median","upper"]
                       : ["event","method","rate","tree","multiplier","bin","count"])
                      .join (" ") + "\n")

var binStart = Array.apply(null,Array(nbins)).map (function(x,n) { return (n*binsize).toPrecision(2) })
function printrange(method,tree,vals) {
    var sortedVals = vals.sort (function (a,b) { return a - b })
    var len = sortedVals.length
    var lq = sortedVals[Math.floor(len/4)], median = sortedVals[Math.floor(len/2)], uq = sortedVals[Math.floor(3*len/4)]
    process.stdout.write ([method,tree,lq,median,uq]
                          .join(" ") + "\n")
}

function printbins(label,headers,bins) {
    var total = 0
    Object.keys(bins).forEach (function (n) {
        total += bins[n]
    })
    binStart.forEach (function (start, n) {
        process.stdout.write ([label]
                              .concat(headers)
                              .concat([start,(bins[n] || 0) / total])
                              .join(" ") + "\n")
    })
}

methods.forEach (function(method) {
    var ibins = {}, dbins = {}, vals = {}
    dirs.forEach (function(dir) {
        var match = /.*rate([\d\.]+)\/([a-z]+)(\d+)/.exec(dir)
        if (match) {
            var rate = match[1],
                tree = match[2],
                mul = match[3],
                filenameRegex = new RegExp('\\d+\\.' + method + '\\.json'),
                headers
            if (!vals[tree])
                vals[tree] = []
            fs.readdirSync(dir)
                .filter (function (filename) { return filenameRegex.test(filename) })
                .map (function (filename) {
                    try {
                        var json = JSON.parse(fs.readFileSync(dir+"/"+filename))
                        incbin(json.insrate/rate,ibins,vals[tree])
                        incbin(json.delrate/rate,dbins,vals[tree])
                    } catch(err) {
                        process.stderr.write(err.toString()+" ("+dir+"/"+filename+")\n")
                    }
                })
            if (!quartileFlag) {
                printbins("ins",[method,rate,tree,mul],ibins)
                printbins("del",[method,rate,tree,mul],dbins)
                ibins = {}
                dbins = {}
            }
        }
    })
    if (quartileFlag) {
        Object.keys(vals).forEach (function (tree) {
            printrange(method,tree,vals[tree])
        })
    }
})

