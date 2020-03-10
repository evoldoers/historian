#!/usr/bin/env node

var prog = process.argv[2]
var opts = process.argv.slice(3)

var hist = require ('../' + prog)
var fs = require ('fs')

hist.onRuntimeInitialized = () => {
  var nFiles = 0, filePrefix = 'FILE'
  opts = opts.map ((opt) => {
    if (opt.charAt(0) !== '-'
	&& fs.existsSync(opt)
	&& fs.lstatSync(opt).isFile()) {
      const fileBuffer = fs.readFileSync(opt)
      const filename = filePrefix + (++nFiles)
      hist.FS.writeFile (filename, new Uint8Array(fileBuffer))
      return filename
    } else
      return opt
  })

  hist.callMain (opts)
}
