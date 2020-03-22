# Historian

evol-historian is an emscripten-compiled version of [historian](https://github.com/evoldoers/historian).

~~~~
var historian = require ('historian')
historian.runWithFiles ('-h').then
 ((result) => console.log (result.stdout))
~~~~
