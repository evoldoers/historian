#!/usr/bin/env node

var fs = require("fs")
var len = process.argv[2], maxIndelLen = process.argv[3], indelProb = process.argv[4], treeFile = process.argv[5]
var tree = fs.readFileSync (treeFile, 'utf8')
process.stdout.write ("[" + len + "]{" + maxIndelLen + "," + indelProb + "}" + tree)
