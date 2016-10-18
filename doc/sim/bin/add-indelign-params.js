#!/usr/bin/env node

var fs = require("fs")
var len = process.argv[2], maxIndelLen = process.argv[3], indelProb = process.argv[4], treeFile = process.argv[5], branchMultiplier = process.argv[6] || 1
var tree = fs.readFileSync (treeFile, 'utf8')

tree = tree.replace (/:([\d\.]+)/g, function (str, p1, offset, s) {
    return ':' + (p1 * branchMultiplier)
})

process.stdout.write ("[" + len + "]{" + maxIndelLen + "," + indelProb + "}" + tree)
