// based on https://github.com/mmalecki/ansispan
var ansi = {}
ansi.resetAnsi = function() { ansi.fg = 'white'; ansi.bg = 'black' }
ansi.resetAnsi()
ansi.foregroundColor = {
  '30': 'black',
  '31': 'red',
  '32': 'green',
  '33': 'yellow',
  '34': 'blue',
  '35': 'purple',
  '36': 'cyan',
  '37': 'white'
};
ansi.backgroundColor = {
  '40': 'black',
  '41': 'red',
  '42': 'green',
  '43': 'yellow',
  '44': 'blue',
  '45': 'purple',
  '46': 'cyan',
  '47': 'white'
};
ansi.ansi2html = function (str) {
  // https://stackoverflow.com/a/5499821
  var tagsToEscape = {
    '&': '&amp;',
    '<': '&lt;',
    '>': '&gt;'
  };

  function escapeTag(tag) {
    return tagsToEscape[tag] || tag;
  }

  function escapeAllTags(str) {
    return str.replace(/[&<>]/g, escapeTag);
  }

  function beginSpan() {
    return '<span style="color: ' + ansi.fg + '; background-color: ' + ansi.bg + '">';
  }

  function endSpan() {
    return '</span>'
  }

  str = escapeAllTags (str);
  const regex = new RegExp('(.*?)\033\\[([0-9]+)m', 'g');
  let span = beginSpan(), match, pos = 0;
  const ansi2html = this;
  while (match = regex.exec(str)) {
    span += match[1];
    pos = match.index + match[0].length;
    const code = match[2], fg = ansi.foregroundColor[code], bg = ansi.backgroundColor[code];
    let changed = false;
    if (code === '0') {
      ansi.resetAnsi();
      changed = true;
    } else if (fg) {
      ansi.fg = fg;
      changed = true;
    } else if (bg) {
      ansi.bg = bg;
      changed = true;
    }
    if (changed)
      span += endSpan() + beginSpan();
  }
  span += str.substr(pos) + endSpan();
  return span;
};

if (typeof module !== 'undefined' && module.exports) {
  module.exports = ansi;
}
