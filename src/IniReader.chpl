/* IniReader.chpl — tiny INI-style parameter-file reader.
 *
 * Syntax:  key = value   one per line; `#` or `;` start a comment;
 * `[section]` headers are allowed for organisation but ignored (keys
 * live in a flat namespace and use the same names as the command-line
 * flags).  See runtime_params.ini at the repository root.
 */
module IniReader {
  use IO, Map, FileSystem;
  import Cli;

  var table = new map(string, string);
  var iniPathUsed = "";      // the parameter file actually read, if any

  proc loadIni() {
    const explicit = !Cli.isUnset(Cli.paramsFile);
    const path = if explicit then Cli.paramsFile else "runtime_params.ini";
    var found = false;
    try! { found = exists(path); }
    if !found {
      if explicit then
        halt("parameter file not found: " + path);
      return;                       // optional default file is absent: fine
    }
    iniPathUsed = path;
    var text: string;
    try! {
      var f = open(path, ioMode.r);
      var rd = f.reader(locking=false);
      text = rd.readAll(string);
      rd.close();
    }
    for line0 in text.split("\n") {
      var line = line0;
      for mark in ("#", ";") {
        const cut = line.find(mark);
        if cut != -1 then line = try! line[..<cut];
      }
      line = line.strip();
      if line == "" || line.startsWith("[") then continue;
      const eq = line.find("=");
      if eq == -1 then
        halt("runtime params: cannot parse line: " + line0);
      const key = (try! line[..<eq]).strip();
      const val = (try! line[(eq+1)..]).strip();
      if key == "" then
        halt("runtime params: empty key in line: " + line0);
      table[key] = val;
    }
  }

  loadIni();   // executed at module initialization, before Params resolves

  proc hasKey(k: string) do return table.contains(k);

  proc getS(k: string): string do return try! table[k];

  proc getI(k: string): int {
    const s = getS(k);
    try {
      return s: int;
    } catch {
      halt("runtime params: '" + k + " = " + s + "' is not an integer");
    }
    return 0;
  }

  proc getR(k: string): real {
    const s = getS(k);
    try {
      return s: real;
    } catch {
      halt("runtime params: '" + k + " = " + s + "' is not a real number");
    }
    return 0.0;
  }
}
