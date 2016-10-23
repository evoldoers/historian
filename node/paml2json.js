#!/usr/bin/env node

var fs = require('fs'),
    extend = require('extend'),
    getopt = require('node-getopt')

var defaults = { insrate: .01,
                 delrate: .01,
                 insextprob: .66,
                 delextprob: .66,
                 alphabet: "arndcqeghilkmfpstwyv" }

var opt = getopt.create([
    ['a' , 'alphabet=STRING'  , 'alphabet (default='+defaults.alphabet+')'],
    ['i' , 'insrate=N'        , 'insertion rate (default='+defaults.insrate+')'],
    ['d' , 'delrate=N'        , 'insertion rate (default='+defaults.delrate+')'],
    ['x' , 'insextprob=N'     , 'insertion extension probability (default='+defaults.insextprob+')'],
    ['y' , 'delextprob=N'     , 'insertion extension probability (default='+defaults.delextprob+')'],
    ['h' , 'help'             , 'display this help message']
])              // create Getopt instance
.bindHelp()     // bind option 'help' to default action
.setHelp(
  "Usage: node " + process.argv[1].match(/(?:.*[\/\\])?(.*)$/)[1] + " <PAML exchangeability matrix>\n" +
  "\n" +
  "[[OPTIONS]]\n"
)
.parseSystem(); // parse command line


function inputError(err) {
    console.warn (err)
    process.exit()
}

opt.argv.length == 1 || inputError("Please specify a PAML exchangeability matrix")
var pamlFile = opt.argv[0]

var model = {}
Object.keys(defaults).forEach (function (arg) {
    model[arg] = (typeof(opt.options[arg]) === 'undefined'
                  ? defaults[arg]
                  : (typeof(defaults[arg] === 'string')
                     ? opt.options[arg]
                     : parseFloat(opt.options[arg])))
})
var alph = model.alphabet.split("")

function nonempty(str) { return /\S/.test(str) }

fs.readFile (pamlFile, 'utf8', (err, data) => {
    if (err) inputError(err)
    else {
        var rows = data.split('\n')
            .filter(nonempty)
            .map((line)=>{ return line.split(' ').filter(nonempty).map(parseFloat) })
        if (rows.length < alph.length)
            inputError ("Input file has " + rows.length + " nonempty rows, but alphabet has " + alph.length + " characters")
        while (rows[alph.length-1].length < alph.length && rows.length > alph.length) {
            rows[alph.length-1] = rows[alph.length-1].concat(rows[alph.length])
            rows.splice (alph.length, 1)
        }
        for (var n = 0; n < alph.length; ++n) {
            if (rows[n].length != n + 1)
                inputError ("Row #" + n + " of input file has " + rows[n].length + " columns; expected " + (n+1) + "\n" + rows[n])
        }
        var eqm = rows[alph.length - 1]
        model.rootprob = {}
        model.subrate = {}
        alph.forEach ((c) => { model.subrate[c] = {} })
        for (var i = 0; i < alph.length; ++i) {
            model.rootprob[alph[i]] = eqm[i]
            for (var j = 0; j < i; ++j) {
                var exch = rows[i-1][j]
                if (exch > 0) {
                    model.subrate[alph[i]][alph[j]] = exch * eqm[j]
                    model.subrate[alph[j]][alph[i]] = exch * eqm[i]
                }
            }
        }
        process.stdout.write (JSON.stringify (model, null, 2))
    }
})
