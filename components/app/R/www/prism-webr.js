// components/app/R/www/prism-webr.js
//
// Client-side, sandboxed R execution for the PRISM board's chartbot.
// Runs LLM-generated ggplot2 code inside webR (R compiled to WASM) instead of
// evaluating it on the server, so arbitrary generated code never runs server-side.
//
// Assumes the "prism" Shiny module is mounted with module id "prism" (see
// components/app/R/ui.R), so every ns()-generated element id used here is the
// literal "prism-" prefix below.
import { WebR } from 'https://webr.r-wasm.org/latest/webr.mjs';

const ID_PREFIX = 'prism-';
const PACKAGES = ['ggplot2', 'ggrepel', 'ggprism', 'ggpubr', 'ggsignif'];

let webR = null;
let isReady = false;

function el(name) {
  return document.getElementById(ID_PREFIX + name);
}

function status(msg, color = '#64748b') {
  const node = el('webr-status');
  if (node) { node.textContent = msg; node.style.color = color; }
}

function setSendEnabled(enabled) {
  const btn = el('chartbot_send');
  if (btn) btn.disabled = !enabled;
}

function showPlot(svgContent) {
  const placeholder = el('plot-placeholder');
  const container = el('plot-container');
  const errorEl = el('plot-error');

  if (placeholder) placeholder.style.display = 'none';
  if (errorEl) errorEl.style.display = 'none';

  container.innerHTML = svgContent;

  const svg = container.querySelector('svg');
  if (svg) {
    svg.style.maxWidth = '100%';
    svg.style.height = 'auto';
    svg.removeAttribute('width');
    svg.removeAttribute('height');
  }
}

function showError(msg) {
  const node = el('plot-error');
  if (node) {
    node.textContent = 'Plot error: ' + msg;
    node.style.display = 'block';
  }
}

// Escape a string for safe embedding inside an R double-quoted string literal.
function escapeRString(s) {
  return s
    .replace(/\\/g, '\\\\')
    .replace(/"/g, '\\"')
    .replace(/\r\n/g, '\\n')
    .replace(/\n/g, '\\n')
    .replace(/\r/g, '\\n');
}

async function initWebR() {
  try {
    setSendEnabled(false);
    status('Initialising webR runtime…');

    webR = new WebR();
    await webR.init();

    status('Installing webR packages…');
    await webR.installPackages(PACKAGES, { quiet: true });
    await webR.evalR(
      'suppressPackageStartupMessages({' +
      PACKAGES.map((p) => `library(${p})`).join(';') +
      '})'
    );

    isReady = true;
    setSendEnabled(true);
    status('webR ready — ask for a plot!', '#059669');
  } catch (err) {
    status('webR error: ' + err.message, '#dc2626');
    console.error('[prism webR] init error:', err);
  }
}

async function executePlotCode(code, csv) {
  if (!isReady) {
    console.warn('[prism webR] not ready yet, ignoring executeCode');
    return;
  }

  status('Rendering plot…', '#64748b');

  const safeCode = code.replace(/\\/g, '\\\\');
  const safeCsv = escapeRString(csv || '');

  const wrapped = `
local({
  data <- read.csv(text = "${safeCsv}", row.names = 1)
  tmp <- tempfile(fileext = ".svg")
  grDevices::svg(tmp, width = 8, height = 5.5, onefile = FALSE)
  on.exit({ if (dev.cur() > 1) dev.off(); unlink(tmp) }, add = TRUE)
  tryCatch(
    {
      .plotbot_result <- withVisible({ ${safeCode} })
      if (.plotbot_result$visible) print(.plotbot_result$value)
    },
    error = function(e) stop(conditionMessage(e))
  )
  dev.off()
  paste(readLines(tmp, warn = FALSE), collapse = "\\n")
})`;

  try {
    const result = await webR.evalR(wrapped);
    const js = await result.toJs();
    const svg = js.values[0];

    showPlot(svg);
    status('Plot rendered!', '#059669');
  } catch (err) {
    showError(err.message);
    addChatMessage({ role: 'error', content: 'Plot error: ' + err.message });
    status('Plot error — see panel below.', '#d97706');
    console.error('[prism webR] plot error:', err);
  }
}

function addChatMessage({ role, content }) {
  const container = el('chat-messages');
  if (!container) return;

  const div = document.createElement('div');
  div.className = `prism-msg prism-msg-${role}`;
  div.textContent = content;

  container.appendChild(div);
  container.scrollTop = container.scrollHeight;
}

function bootstrap() {
  if (!window.Shiny || !window.Shiny.addCustomMessageHandler) {
    setTimeout(bootstrap, 150);
    return;
  }

  Shiny.addCustomMessageHandler('prism-executeCode', ({ code, csv }) => executePlotCode(code, csv));
  Shiny.addCustomMessageHandler('prism-addMessage', (data) => addChatMessage(data));

  initWebR();
}

if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', bootstrap);
} else {
  bootstrap();
}
