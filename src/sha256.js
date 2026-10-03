/**
 * sha256.js — the one content-hash the results record uses, in ~70 lines of
 * dependency-free JavaScript.
 *
 * WHY THIS EXISTS INSTEAD OF `crypto.subtle.digest`
 * -------------------------------------------------
 * The results record has to carry a content hash of the input structure, and it
 * has to follow the convention `data/manifest.json` already established:
 *
 *     "4w52.pdb": "sha256:<64 hex chars>"
 *
 * (bench/require_inputs.js:59 `sha256File` computes exactly that, over the raw
 * file bytes, and refuses a bench whose input hash moved.)
 *
 * The WebCrypto options were both unusable for a record that must work offline,
 * from `file://`, and in plain Node:
 *   - `crypto.subtle` is async and only exists in a SECURE CONTEXT. `file://`
 *     pages are not one, so `sha256(4w52.pdb)` would be `null` for the user who
 *     opened the app by double-clicking index.html — i.e. silently missing for
 *     exactly the offline-first user this project targets.
 *   - `node:crypto` does not exist in a browser and importing it would make the
 *     record module non-importable under Node, which the round-trip test needs.
 * So the record would either carry no hash (unreproducible input identity) or
 * carry a DIFFERENT, home-made one that could not be compared against the
 * manifest. This module is the third option: a real SHA-256, synchronous,
 * identical in both runtimes.
 *
 * IT IS NOT A SECURITY PRIMITIVE AND IS NOT USED AS ONE
 * ---------------------------------------------------
 * Nothing here is constant-time and nothing here resists an adversary. It exists
 * to answer "are these two runs looking at the same bytes?", which is what the
 * manifest already uses it for.
 *
 * CORRECTNESS IS ASSERTED, NOT ASSERTED-BY-COMMENT
 * -----------------------------------------------
 * tests/test_results_record.js checks (a) the FIPS-180-4 vectors for "" and
 * "abc", (b) every one of the 64 round constants against a fixed table, and
 * (c) the digest of the repo's real 4w52.pdb against the value already recorded
 * in data/manifest.json — so a broken compression function cannot pass by
 * looking plausible.
 */

/** First 32 primes cubed: the fractional parts of the cube roots of the first
 *  64 primes, as specified by FIPS 180-4 §4.2.2. Doubled so JS bit-ops can
 *  stay in 32-bit integer land (>>> 0 keeps them unsigned). */
const K = new Uint32Array([
  0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
  0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
  0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
  0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
  0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
  0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
  0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
  0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]);

const rotr = (x, n) => (x >>> n) | (x << (32 - n));

/**
 * SHA-256 digest of a byte array, as 64 lowercase hex characters.
 * @param {Uint8Array} bytes input bytes (NOT a string — use `sha256Text`)
 * @returns {string} 64 hex chars, no `sha256:` prefix
 */
export function sha256Bytes(bytes) {
  if (!(bytes instanceof Uint8Array)) throw new Error("sha256Bytes: expected a Uint8Array.");
  const bitLen = bytes.length * 8;
  // Message + 0x80 + zero padding + 8-byte big-endian length, to a 64-byte multiple.
  const padded = new Uint8Array(((bytes.length + 9 + 63) >> 6) << 6);
  padded.set(bytes);
  padded[bytes.length] = 0x80;
  const dv = new DataView(padded.buffer);
  // Length in BITS as a 64-bit big-endian value. Split so the high word stays
  // exact past 2^32 bits (512 MiB) instead of silently wrapping.
  dv.setUint32(padded.length - 8, Math.floor(bitLen / 0x100000000), false);
  dv.setUint32(padded.length - 4, bitLen >>> 0, false);

  const H = new Uint32Array([
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
  ]);
  const w = new Uint32Array(64);
  for (let off = 0; off < padded.length; off += 64) {
    for (let i = 0; i < 16; i++) w[i] = dv.getUint32(off + 4 * i, false);
    for (let i = 16; i < 64; i++) {
      const s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >>> 3);
      const s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >>> 10);
      w[i] = (w[i - 16] + s0 + w[i - 7] + s1) >>> 0;
    }
    let [a, b, c, d, e, f, g, h] = H;
    for (let i = 0; i < 64; i++) {
      const S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
      const ch = (e & f) ^ (~e & g);
      const t1 = (h + S1 + ch + K[i] + w[i]) >>> 0;
      const S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
      const maj = (a & b) ^ (a & c) ^ (b & c);
      const t2 = (S0 + maj) >>> 0;
      h = g; g = f; f = e;
      e = (d + t1) >>> 0;
      d = c; c = b; b = a;
      a = (t1 + t2) >>> 0;
    }
    H[0] = (H[0] + a) >>> 0; H[1] = (H[1] + b) >>> 0;
    H[2] = (H[2] + c) >>> 0; H[3] = (H[3] + d) >>> 0;
    H[4] = (H[4] + e) >>> 0; H[5] = (H[5] + f) >>> 0;
    H[6] = (H[6] + g) >>> 0; H[7] = (H[7] + h) >>> 0;
  }
  let hex = "";
  for (let i = 0; i < 8; i++) hex += H[i].toString(16).padStart(8, "0");
  return hex;
}

/**
 * SHA-256 of a string's UTF-8 bytes, prefixed `sha256:` — the exact shape
 * `data/manifest.json` uses, so a results record and the manifest can be
 * compared with `===`.
 *
 * String, not bytes, on purpose: the browser hands this module the PDB *text* it
 * parsed, and that text is what the run used. Hashing the file on disk instead
 * would be a different (and, for an RCSB fetch, unknowable) byte sequence.
 * @param {string} text
 * @returns {string} "sha256:<64 hex>"
 */
export function sha256Text(text) {
  const s = typeof text === "string" ? text : "";
  return "sha256:" + sha256Bytes(new TextEncoder().encode(s));
}
