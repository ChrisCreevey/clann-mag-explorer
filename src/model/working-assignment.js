(function () {
  'use strict';

// Mutable session state for interactive reassignment (brief §Interactive
// reassignment): a single Map<contigId, binId> that starts from whichever
// bin information was loaded and gets edited as the student moves
// contigs around. Every reassignment action in the brief — move a
// contig, split a bin, merge two bins, spin a new bin out of a flagged
// cluster — reduces to the same primitive: reassign a set of contig IDs
// to a target bin ID (existing or brand new). Kept as pure
// state-transition functions (old Map in, new Map out) rather than
// mutating in place, so the UI layer can diff/undo/re-render off a
// simple reference-equality check if it ever needs to.

/**
 * The starting point for the editable assignment: the reconciled
 * majority vote when 2+ tools were loaded (brief's "high-confidence
 * core" plus best-guess placement for disputed contigs — a contig with
 * no majority at all, e.g. every tool split evenly, has no starting
 * bin and stays unassigned until the student places it), or the single
 * table's own bins when only one tool was loaded. No bin table at all
 * means nothing to reassign yet — an empty Map.
 * @param {Map<string,{contigId:string,binId:string}[]>|null} binTablesByTool
 * @param {{contigAgreement: {contigId:string, majorityMagId:string|null}[]}|null} reconciliation
 * @returns {Map<string,string>} contigId -> binId
 */
function deriveInitialAssignment(binTablesByTool, reconciliation) {
  const assignment = new Map();
  if (reconciliation) {
    for (const c of reconciliation.contigAgreement) {
      if (c.majorityMagId) assignment.set(c.contigId, c.majorityMagId);
    }
    return assignment;
  }
  if (binTablesByTool && binTablesByTool.size === 1) {
    const [[, table]] = binTablesByTool;
    for (const { contigId, binId } of table) assignment.set(contigId, binId);
  }
  return assignment;
}

/**
 * The one reassignment primitive everything else is built from: move
 * `contigIds` to `targetBinId`, leaving every other contig's assignment
 * untouched. Returns a new Map — the caller decides whether/how to track
 * that something changed (e.g. arming the beforeunload guard).
 */
function reassignContigs(assignment, contigIds, targetBinId) {
  const next = new Map(assignment);
  for (const contigId of contigIds) next.set(contigId, targetBinId);
  return next;
}

/** All contig IDs currently assigned to `binId`. */
function contigsInBin(assignment, binId) {
  const ids = [];
  for (const [contigId, b] of assignment) if (b === binId) ids.push(contigId);
  return ids;
}

/** Every distinct bin ID currently in use, sorted. */
function listBinIds(assignment) {
  return [...new Set(assignment.values())].sort();
}

/** A bin ID guaranteed not already in use, for "spin out a new bin" actions. */
function generateNewBinId(assignment, prefix = 'bin.new') {
  const existing = new Set(assignment.values());
  if (!existing.has(prefix)) return prefix;
  let i = 2;
  while (existing.has(`${prefix}.${i}`)) i++;
  return `${prefix}.${i}`;
}

/** Merge every contig currently in `sourceBinId` into `targetBinId`. */
function mergeBins(assignment, sourceBinId, targetBinId) {
  return reassignContigs(assignment, contigsInBin(assignment, sourceBinId), targetBinId);
}

/** Convert back to the {contigId, binId}[] shape bin-summary.js/bin-reconciliation.js expect. */
function assignmentToRows(assignment) {
  return [...assignment.entries()].map(([contigId, binId]) => ({ contigId, binId }));
}

const exportsObj = {
  deriveInitialAssignment, reassignContigs, contigsInBin, listBinIds,
  generateNewBinId, mergeBins, assignmentToRows,
};
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.workingAssignment = exportsObj;
}
})();
