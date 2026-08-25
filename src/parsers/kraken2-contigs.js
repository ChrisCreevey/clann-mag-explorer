(function () {
  'use strict';

// Per-contig taxonomic call (brief §Inputs: "Reuses the eDNA Explorer's
// Kraken2 report parser, run against the assembly's contigs rather than
// reads"). Scope decision, flagged rather than silently assumed: the
// eDNA Explorer's parser (breport.js, already ported — see
// docs/phase1-investigation.md §2/§Phase 4) reads Kraken2's *aggregate*
// .breport (one row per taxon in the whole-run hierarchy, read counts
// summed across every query) — it has no per-query field at all, so it
// structurally cannot answer "what did Kraken2 call THIS contig." What
// actually carries a per-query taxonomic call is Kraken2's other, simpler
// standard output format (`--output`, distinct from `--report`): one row
// per query, `C|U <TAB> seq_id <TAB> taxid[ (name)] <TAB> length <TAB>
// LCA mapping`. That's what this module parses. breport.js's *tree/LCA
// machinery* (taxonomy-tree.js) is still reused for the taxonomic-
// disagreement check itself, per the brief's actual intent — just fed
// from this per-contig format instead of a breport, since only this one
// carries a call per contig.

const TAXID_IN_NAME_RE = /\(taxid (\d+)\)/; // --use-names format: "Escherichia coli (taxid 562)"

function parseTaxIdField(field) {
  if (/^\d+$/.test(field)) return Number(field);
  const m = TAXID_IN_NAME_RE.exec(field);
  return m ? Number(m[1]) : null;
}

/**
 * @param {string} text
 * @returns {{contigId: string, classified: boolean, taxId: number|null}[]}
 */
function parseKraken2ContigCalls(text) {
  const calls = [];
  for (const rawLine of text.split(/\r?\n/)) {
    if (!rawLine.trim()) continue;
    const fields = rawLine.split('\t');
    if (fields.length < 4) continue;
    const [classFlag, contigId, taxIdField] = fields;
    if (classFlag !== 'C' && classFlag !== 'U') continue;
    calls.push({ contigId, classified: classFlag === 'C', taxId: parseTaxIdField(taxIdField) });
  }
  return calls;
}

/** Cheap content-based check for sniff.js. */
function looksLikeKraken2ContigOutput(text) {
  const lines = text.split(/\r?\n/).filter((l) => l.trim().length > 0).slice(0, 20);
  if (lines.length === 0) return false;
  return lines.every((line) => {
    const fields = line.split('\t');
    if (fields.length < 4) return false;
    if (fields[0] !== 'C' && fields[0] !== 'U') return false;
    return parseTaxIdField(fields[2]) !== null && !Number.isNaN(Number(fields[3]));
  });
}

const exportsObj = { parseKraken2ContigCalls, looksLikeKraken2ContigOutput };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.kraken2Contigs = exportsObj;
}
})();
