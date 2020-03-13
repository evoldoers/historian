// This code goes at the start of the generated historian JS
Module.noInitialRun = true;
Module.noExitRuntime = true;

Module.print = function (x) {
  if (Module.stdout)
    Module.stdout.push (x)
  else
    console.log (x)
}

Module.printErr = function (x) {
  if (Module.stderr)
    Module.stderr.push (x)
  else
    console.log (x)
}

Module.runtimeInitPromise = new Promise ((resolve, reject) => {
  Module.onRuntimeInitialized = resolve
})

Module.runWithFiles = (args) => {
  return Module.runtimeInitPromise
    .then (() => {
      args = args || [];

      let nFiles = 0, filePrefix = 'FILE', outputs = []
      const wrappedArgs = args.map ((opt) => {
	if (typeof(opt) !== 'string') {
	  const filename = opt.filename || (filePrefix + (++nFiles))
          if (opt.input || opt.data) {
            const fileBuffer = Buffer.from (opt.data || '', 'utf-8')
	    Module.FS.writeFile (filename, new Uint8Array(fileBuffer))
          }
          if (opt.output)
            outputs.push (filename)
	  return filename
	} else
	  return opt
      })

      const oldStdout = Module.stdout
      const oldStderr = Module.stderr
      
      Module.stdout = []
      Module.stderr = []
      Module.callMain (wrappedArgs)
      const stdout = Module.stdout.join('')
      const stderr = Module.stderr.join('')
      Module.stdout = oldStdout
      Module.stderr = oldStderr

      let file = {}
      outputs.forEach ((filename) => {
        file[filename] = Module.FS.readFile (filename, { encoding: 'utf8'})
      })
      return { stdout, stderr, outputs, file }
    })
}

// Web worker
onmessage = function(e) {
  Module.runWithFiles(e.data).then ((result) => {
    postMessage (result);
  })
}
