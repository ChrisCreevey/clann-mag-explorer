(function () {
  'use strict';

// Standard BLOSUM62 substitution matrix (Henikoff & Henikoff 1992),
// used by marker-genes.js for real-residue extension scoring (brief
// §Marker-gene identification module, step 3: "extension proceeds left
// and right using the actual BLOSUM62 matrix on the real (non-reduced)
// residues"). The reduced alphabet (reduced-alphabet.js) is only for the
// seeding lookup — scoring always happens on real residues, per the
// brief.

const ORDER = 'ARNDCQEGHILKMFPSTWYV';
// prettier-ignore
const ROWS = [
  [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0], // A
  [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3], // R
  [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3], // N
  [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3], // D
  [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1], // C
  [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2], // Q
  [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2], // E
  [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3], // G
  [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3], // H
  [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3], // I
  [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1], // L
  [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2], // K
  [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1], // M
  [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1], // F
  [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2], // P
  [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2], // S
  [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0], // T
  [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3], // W
  [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1], // Y
  [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4], // V
];

const X_SCORE = -1; // BLOSUM62 convention: X (any) scores roughly the matrix average
const STOP_SCORE = -100; // '*' or any other non-standard byte: force extension to halt here,
                          // rather than a special-cased boundary check (brief's "no special
                          // case needed" philosophy for sentinel-delimited segmentation)

// score[a*128+b] for charCodes a,b in [0,128). Built once.
const BLOSUM62_SCORE = new Int8Array(128 * 128).fill(STOP_SCORE);
for (let i = 0; i < 128; i++) {
  const code = 'X'.charCodeAt(0);
  BLOSUM62_SCORE[i * 128 + code] = X_SCORE;
  BLOSUM62_SCORE[code * 128 + i] = X_SCORE;
}
for (let i = 0; i < ORDER.length; i++) {
  for (let j = 0; j < ORDER.length; j++) {
    const a = ORDER.charCodeAt(i), b = ORDER.charCodeAt(j);
    BLOSUM62_SCORE[a * 128 + b] = ROWS[i][j];
  }
}

const exportsObj = { BLOSUM62_SCORE, ORDER };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.blosum62 = exportsObj;
}
})();
