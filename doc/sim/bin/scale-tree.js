#!/usr/bin/env node

var fs = require("fs")
var treeFile = process.argv[2], branchMultiplier = process.argv[3] || 1
var tree = fs.readFileSync (treeFile, 'utf8')

tree = tree.replace (/:([\d\.]+)/g, function (str, p1, offset, s) {
    return ':' + (p1 * branchMultiplier)
})

process.stdout.write (tree)
