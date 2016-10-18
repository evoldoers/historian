#!/usr/bin/env node

var fs = require("fs")
var execSync = require("child_process").execSync
function exec (cmd) {
    console.log(cmd)
    execSync(cmd)
}

var pow = process.argv[2], dir = process.argv[3], makefileFromDir = process.argv[4], indelProb = process.argv[5], reps = process.argv[6], treeFile = process.argv[7]

var len = Math.pow (10, pow)
var maxIndelLen = Math.pow (2, parseInt(pow) + 2)

var paramFile = dir + "/params.isg"

exec ("bin/add-indelign-params.js " + len + " " + maxIndelLen + " " + indelProb + " " + treeFile + " >" + paramFile)
for (var rep = 1; rep <= reps; ++rep) {
    exec ("cd " + dir + "; make -f " + makefileFromDir + " dna" + rep + ".hist.json")
}
