(function () {
  'use strict';

// App shell: theme toggle and file-picker wiring only, for now. The load
// pipeline (sniffing, streaming FASTA parse, bin table loading) gets wired
// in as those modules land per the suggested build phases (brief §Suggested
// build phases). Structure follows clann-edna-explorer/src/app.js.

const THEME_KEY = 'clann-mag-explorer-theme';

function initTheme() {
  const saved = localStorage.getItem(THEME_KEY);
  if (saved === 'light' || saved === 'dark') {
    document.documentElement.setAttribute('data-theme', saved);
  }
  const btn = document.getElementById('themeBtn');
  btn.addEventListener('click', () => {
    const current = document.documentElement.getAttribute('data-theme');
    const prefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches;
    const isDark = current ? current === 'dark' : prefersDark;
    const next = isDark ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', next);
    localStorage.setItem(THEME_KEY, next);
  });
}

function initFilePicker() {
  const input = document.getElementById('folder-input');
  const openButtons = [document.getElementById('uploadBtn'), document.getElementById('emptyOpen')];
  openButtons.forEach((btn) => btn && btn.addEventListener('click', () => input.click()));
  input.addEventListener('change', () => {
    // Sniffing + load pipeline lands with src/parsers/sniff.js's real
    // implementation (Phase 1/2). For now, just acknowledge the selection.
    const err = document.getElementById('err');
    err.style.display = 'block';
    err.textContent = `${input.files.length} file(s) selected — load pipeline not implemented yet.`;
  });
}

document.addEventListener('DOMContentLoaded', () => {
  initTheme();
  initFilePicker();
});
})();
