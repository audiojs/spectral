// Writes pvoc.wat's WebAssembly bytes into pvoc.js: wat2wasm pvoc.wat -o pvoc.wasm && node scripts/wasm.js
import { readFileSync, writeFileSync, rmSync } from 'node:fs'
let here = new URL('../', import.meta.url), bytes = readFileSync(new URL('pvoc.wasm', here)), js = readFileSync(new URL('pvoc.js', here), 'utf8')
let rows = [], line = '  '
for (let b of bytes) { let s = b + ','; if (line.length + s.length > 120) { rows.push(line); line = '  ' } line += s }
rows.push(line)
let out = js.replace(/(const WASM = new Uint8Array\(\[)[^]*?(\]\))/, `$1\n${rows.join('\n')}\n$2`)
  .replace(/\/\/ pvoc\.wat, compiled by wat2wasm \([\d,]+ bytes\)/, `// pvoc.wat, compiled by wat2wasm (${bytes.length.toLocaleString('en')} bytes)`)
writeFileSync(new URL('pvoc.js', here), out)
rmSync(new URL('pvoc.wasm', here))
console.log(`pvoc.js: ${bytes.length} bytes`)
