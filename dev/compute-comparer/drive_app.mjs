#!/usr/bin/env node
//
// Path A: drive the real Shiny upload flow to the point where it writes
// params.RData, then copy that file out before the app deletes it.
//
//   node drive_app.mjs <scenario-id> [--headed] [--url http://localhost:3838]
//
// See docs/app-automation.md for why we stop at params.RData and for the
// input-id reference this script encodes.

import { chromium } from 'playwright';
import { readFileSync, existsSync, mkdirSync, copyFileSync, readdirSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = dirname(fileURLToPath(import.meta.url));
const OPG = join(HERE, '..', '..');            // omicsplayground root (master worktree)
const USER_INPUT = join(OPG, 'data', 'USER_INPUT');

const argv = process.argv.slice(2);
const scenarioId = argv.find((a) => !a.startsWith('--'));
const headed = argv.includes('--headed');
const url = (argv.find((a) => a.startsWith('--url')) || '--url=http://localhost:3838').split('=')[1]
  || 'http://localhost:3838';

if (!scenarioId) {
  console.error('usage: node drive_app.mjs <scenario-id> [--headed] [--url=http://host:port]');
  process.exit(2);
}

// ---------------------------------------------------------------- scenario

const cfg = JSON.parse(readFileSync(join(HERE, 'scenarios.json'), 'utf8'));
const raw = cfg.scenarios.find((s) => s.id === scenarioId);
if (!raw) {
  console.error(`unknown scenario "${scenarioId}". known: ${cfg.scenarios.map((s) => s.id).join(', ')}`);
  process.exit(2);
}
// one level of merge is all the schema needs (defaults -> scenario, with qc/compute merged)
const sc = {
  ...cfg.defaults, ...raw,
  qc: { ...cfg.defaults.qc, ...(raw.qc || {}) },
  compute: { ...cfg.defaults.compute, ...(raw.compute || {}) },
};
const fx = cfg.fixtures[sc.fixture];
if (!fx) { console.error(`unknown fixture "${sc.fixture}"`); process.exit(2); }

const outDir = join(HERE, 'runs', scenarioId, 'A');
mkdirSync(outDir, { recursive: true });

const log = (...m) => console.log('[drive]', ...m);

// ---------------------------------------------------------------- helpers

/** Set a Shiny input server-side. Works for widgets with no practical DOM path
 *  (server-side selectize, sliders) and for anything inside a collapsed accordion. */
const setInput = (page, id, value) =>
  page.evaluate(([i, v]) => window.Shiny.setInputValue(i, v, { priority: 'event' }), [id, value]);

/** Wait until Shiny has no pending work.
 *  NB: do NOT look at `.recalculating` -- outputs on tabs that were never shown
 *  keep that class forever (33 of them sit there permanently on the welcome tab),
 *  so any wait including it never returns. `.shiny-busy` is the real signal. */
async function settle(page, quietMs = 1500, timeout = 300000) {
  const deadline = Date.now() + timeout;
  let quietSince = null;
  while (Date.now() < deadline) {
    const busy = await page.evaluate(() =>
      document.querySelectorAll('.shiny-busy').length > 0
      || document.documentElement.classList.contains('shiny-busy'));
    if (busy) { quietSince = null; } else if (quietSince === null) { quietSince = Date.now(); }
    if (quietSince !== null && Date.now() - quietSince >= quietMs) return;
    await page.waitForTimeout(250);
  }
  // Warn, don't throw. settle() is a best-effort quiet-wait; the hard gates are the
  // explicit waitForSelector calls. A busy app (big fixture, several sessions) can
  // legitimately stay busy longer than any fixed budget, and aborting the whole
  // capture over that loses a 10-minute run for no good reason.
  log(`settle() gave up after ${Math.round(timeout / 1000)}s of activity, continuing`);
}

/** Click, falling back to a direct dispatch.
 *
 *  Playwright's click waits for actionability, which includes *stability* -- the
 *  element must not move for two consecutive frames. Shiny re-renders these steps
 *  continuously after a CSV lands, so a perfectly visible, enabled button can fail
 *  actionability forever. dispatchEvent skips those checks and fires the handler
 *  directly, which is all these buttons need. */
async function clickOrDispatch(page, selector, label, timeout = 60000) {
  const loc = page.locator(selector).filter({ visible: true }).first();
  try {
    await loc.click({ timeout });
  } catch (e) {
    log(`${label} click intercepted, dispatching directly: ${e.message.split('\n')[0]}`);
    await clearOverlays(page);
    await settle(page);
    await page.locator(selector).first().dispatchEvent('click');
  }
}

/** Close anything covering the page: the chirp promo modal (ENABLE_CHIRP=TRUE in
 *  etc/OPTIONS pops one over the upload board) and any shinyalert. Deliberately
 *  spares wizard-modal-*, which is the wizard itself. */
async function clearOverlays(page) {
  const n = await page.evaluate(() => {
    let closed = 0;
    document.querySelectorAll('.modal.show').forEach((m) => {
      if (m.id && m.id.startsWith('wizard-modal-')) return;
      const inst = window.bootstrap?.Modal?.getInstance(m);
      if (inst) inst.hide();
      else m.querySelector('.btn-close, [data-bs-dismiss=modal]')?.click();
      closed++;
    });
    document.querySelectorAll('.swal-overlay--show-modal .swal-button').forEach((b) => {
      b.click(); closed++;
    });
    // a hidden modal can leave its backdrop behind, which still eats clicks
    if (!document.querySelector('.modal.show')) {
      document.querySelectorAll('.modal-backdrop').forEach((b) => b.remove());
    }
    return closed;
  });
  if (n > 0) { log(`dismissed ${n} overlay(s)`); await page.waitForTimeout(1000); }
  return n;
}

/** The wizard never locks (see docs/app-automation.md sec.4), so advance on our own
 *  readiness signal instead: wait for `probe` to exist, then click next.
 *  NB: state 'attached', not the default 'visible' -- Shiny's selectInput renders a
 *  selectize widget and hides the underlying <select>, so a bare select id is never
 *  "visible" no matter how long you wait. */
async function nextStep(page, probeSelector, label) {
  if (probeSelector) {
    log(`waiting for ${label}: ${probeSelector}`);
    await page.waitForSelector(probeSelector, { state: 'attached', timeout: 300000 });
  }
  await settle(page);
  await clearOverlays(page);
  // Generous timeout, and pick the visible button: after a CSV lands the step
  // re-renders, and the default 30s is not enough for the contrasts step on a
  // realistically sized fixture. Retry once after clearing overlays again.
  await clickOrDispatch(page, '.wizard-btn.btn.next', 'next', 60000);
  await settle(page);
  log(`advanced past ${label}`);
}

/** Set a select and verify it actually stuck, retrying if the server overwrote it.
 *
 *  Needed for bec_method: the batch-correction comparison calls
 *  `updateSelectInput(session, "bec_method", choices = methods)` with no `selected=`
 *  (upload_module_normalization.R:253-257), and Shiny then selects the FIRST choice --
 *  "ComBat". So any bec_method we set before that comparison finishes is silently
 *  reverted, and the app computes ComBat while the scenario asked for limma/NPM.
 *  (That reset is an app bug in its own right: a user's chosen BC method is discarded
 *  whenever the panel refreshes.) */
async function setInputSticky(page, id, value, tries = 6) {
  for (let i = 0; i < tries; i++) {
    // Drive the selectize instance, not Shiny.setInputValue. selectize keeps its own
    // state and does NOT write through to the hidden original <select>, so reading
    // `select.value` reports a stale value and a plain setInputValue leaves the
    // widget and the server disagreeing. `.selectize.setValue()` fires the real
    // change event, which is the path a user's click takes.
    const got = await page.evaluate(([i2, v]) => {
      const el = document.querySelector(`#${i2}`);
      if (!el) return null;
      if (el.selectize) { el.selectize.setValue(v, false); return el.selectize.getValue(); }
      el.value = v;
      el.dispatchEvent(new Event('change', { bubbles: true }));
      return el.value;
    }, [id, value]);
    await settle(page, 2500);
    const now = await page.evaluate((i2) => {
      const el = document.querySelector(`#${i2}`);
      return el?.selectize ? el.selectize.getValue() : el?.value;
    }, id);
    // Require the value to STAY put: the batch-correction comparison takes far
    // longer than a settle() and resets the select when it finally finishes, so a
    // single post-set check can pass and still be overwritten seconds later.
    if (now === String(value)) {
      let stable = true;
      for (let k = 0; k < 5; k++) {
        await page.waitForTimeout(3000);
        const v = await page.evaluate((i2) => {
          const el = document.querySelector(`#${i2}`);
          return el?.selectize ? el.selectize.getValue() : el?.value;
        }, id);
        if (v !== String(value)) { log(`${id} drifted to "${v}" after ${(k + 1) * 3}s`); stable = false; break; }
      }
      if (stable) return true;
      continue;
    }
    log(`${id} set to "${got}" but now reads "${now}", retrying (${i + 1}/${tries})`);
    await page.waitForTimeout(3000);
  }
  log(`WARNING: ${id} would not hold value "${value}"`);
  return false;
}

/** Upload a CSV. The step's UI is rendered server-side and can take well over
 *  setInputFiles' 30s default to appear (step 3 builds the comparison builder from
 *  the sample table first), so wait for the input explicitly with a long timeout. */
async function uploadFile(page, selector, path, label) {
  log(`uploading ${label}`);
  await page.waitForSelector(selector, { state: 'attached', timeout: 300000 });
  await page.setInputFiles(selector, path, { timeout: 120000 });
  await settle(page);
}

/** Snapshot the raw_* dirs that already exist, so we only pick up OUR run. */
const snapshotRawDirs = () =>
  new Set(existsSync(USER_INPUT) ? readdirSync(USER_INPUT).filter((d) => d.startsWith('raw_')) : []);

/** Poll for params.RData, re-triggering compute periodically.
 *
 *  The retry is not paranoia: the compute observer bails out (returning NULL, with no
 *  visible feedback) while `probetype()` is still "running" -- it is an async
 *  ExtendedTask kicked off by the counts upload. Re-sending `wizard_finished` re-fires
 *  the observer, so once detection lands the run proceeds. raw_dir() is already set by
 *  the uploads and is reused, so retries cannot produce a second directory. */
async function captureParams(before, retrigger, timeoutMs = 600000, retryMs = 15000) {
  const deadline = Date.now() + timeoutMs;
  let lastTrigger = 0;
  while (Date.now() < deadline) {
    if (existsSync(USER_INPUT)) {
      for (const d of readdirSync(USER_INPUT)) {
        if (!d.startsWith('raw_')) continue;
        const p = join(USER_INPUT, d, 'params.RData');
        if (existsSync(p)) return join(USER_INPUT, d);
      }
    }
    if (Date.now() - lastTrigger >= retryMs) {
      lastTrigger = Date.now();
      await retrigger().catch(() => {});
    }
    await new Promise((r) => setTimeout(r, 2000));
  }
  return null;
}

// ---------------------------------------------------------------- main

const browser = await chromium.launch({ headless: !headed });
const page = await browser.newPage({ viewport: { width: 1600, height: 1000 } });
page.on('console', (m) => { if (m.type() === 'error') log('page error:', m.text()); });

let rawDir = null;
try {
  log(`scenario=${scenarioId} fixture=${sc.fixture} (${fx.datatype}) url=${url}`);

  await page.goto(url, { waitUntil: 'domcontentloaded', timeout: 120000 });
  await page.waitForFunction(() => window.Shiny && window.Shiny.shinyapp
    && window.Shiny.shinyapp.$socket, { timeout: 180000 });
  await settle(page, 2000);
  log('app connected');

  // --- sign-in splash ---------------------------------------------------
  // Even with AUTHENTICATION=none the app shows splashLoginModal with a single
  // click-through button (AuthenticationModule.R:29-40). Matched by its label
  // rather than its id (ns("login_emailSubmit")) so the auth namespace can move.
  // Wait for it to appear at all -- on a cold app the modal can lag the socket by
  // several seconds -- then click via the resilient path.
  // The button id differs by branch (edgy: auth-login_submit_btn, master:
  // ns("login_emailSubmit")), so match on the label and fall back to a direct
  // dispatch: the modal container itself reports no layout, so Playwright's
  // actionability check on the button can never settle even though the button is
  // painted and clickable.
  await page.waitForSelector('.modal.show, #shiny-modal', { timeout: 60000 }).catch(() => {});
  const clicked = await page.evaluate(() => {
    const btn = [...document.querySelectorAll('button, a.btn, .btn')]
      .find((e) => /sure i am/i.test(e.textContent || ''));
    if (!btn) return null;
    btn.click();
    return btn.id || '(no id)';
  });
  if (clicked) {
    await settle(page, 2000);
    log(`passed sign-in splash (${clicked})`);
  } else {
    log('no sign-in splash present');
  }

  // --- navigate to the upload board (bigdash tab trigger) ---------------
  // Navigation differs by branch:
  //   master: bigdash bigTabItem("upload-tab") reached via a .tab-trigger that sits in
  //           a closed navbar dropdown -- not clickable with a real mouse event, but
  //           jQuery .trigger('click') fires bigdash's handler regardless.
  //   edgy:   the feat/bigdash-modules rework moved it to a bslib
  //           nav_panel_hidden("Upload") inside the "app-sidebar" navset. A hidden
  //           panel has no nav link to click, so drive the navset's Shiny input
  //           binding directly -- that is exactly what bslib::nav_select() sends.
  const nav = await page.evaluate(() => {
    const bigdash = window.$('.tab-trigger[data-target="upload-tab"]');
    if (bigdash.length) { bigdash.trigger('click'); return 'bigdash'; }
    const el = document.getElementById('app-sidebar');
    const binding = el && window.$(el).data('shiny-input-binding');
    if (binding && binding.receiveMessage) {
      binding.receiveMessage(el, { value: 'Upload' });
      return 'bslib';
    }
    return null;
  });
  log(`navigation: ${nav ?? 'NONE FOUND'}`);
  // Gate on a genuinely visible control, not on a selectize-hidden <select>.
  await page.waitForSelector('#upload-start_upload', { state: 'visible', timeout: 120000 });
  await settle(page);
  await clearOverlays(page);
  log('on upload board');

  // --- landing tab: datatype + organism --------------------------------
  await setInput(page, 'upload-selected_datatype', fx.datatype);
  await settle(page);
  await setInput(page, 'upload-selected_organism', fx.organism);
  await settle(page);

  await clearOverlays(page);
  await clickOrDispatch(page, '#upload-start_upload', 'start_upload', 60000);
  await page.waitForSelector('.wizard-step', { timeout: 60000 });
  await settle(page);
  log('wizard opened');

  // --- step 1: counts ---------------------------------------------------
  await uploadFile(page, '#upload-counts_preview-counts_csv', join(fx.dir, 'counts.csv'), 'counts.csv');
  await nextStep(page, null, 'step 1 (counts)');

  // --- step 2: samples --------------------------------------------------
  await uploadFile(page, '#upload-samples_preview-samples_csv', join(fx.dir, 'samples.csv'), 'samples.csv');
  await nextStep(page, null, 'step 2 (samples)');

  // --- step 3: contrasts ------------------------------------------------
  // show_comparison_builder is a reactiveVal initialised to TRUE
  // (upload_server.R:32), so this step opens on the INTERACTIVE builder and the
  // fileInputArea is not rendered at all. "Upload comparisons file"
  // (input$goUploadComparison, contrasts module:109) flips it to FALSE and reveals
  // the upload area. Without this the file input never appears.
  await clickOrDispatch(page, '#upload-contrasts_preview-goUploadComparison',
    'goUploadComparison', 60000);
  await settle(page);
  log('switched step 3 to upload mode');
  await uploadFile(page, '#upload-contrasts_preview-contrasts_csv', join(fx.dir, 'contrasts.csv'), 'contrasts.csv');
  await nextStep(page, null, 'step 3 (contrasts)');

  // --- step 4: QC / normalization --------------------------------------
  // Wait for the QC controls only AFTER advancing. Shiny suspends outputs that are
  // hidden (suspendWhenHidden defaults TRUE), and the QC uiOutput lives in a wizard
  // step that is display:none until you get there -- so it does not render at all
  // while we are still on step 3. Waiting for it before clicking next deadlocks.
  await page.waitForSelector('#upload-checkqc-normalization_method',
    { state: 'attached', timeout: 300000 });
  await settle(page);
  log('QC panel rendered');

  // Order matters: the enabling flag before the value it gates, so any server-side
  // update() on the dependent widget lands before we assert our value.
  const q = sc.qc;
  const qcOrder = [
    ['zero_as_na', q.zero_as_na],
    ['filtermissing', q.filtermissing],
    ['filterthreshold', q.filterthreshold],
    ['impute', q.impute],
    ['impute_method', q.impute_method],
    ['normalize', q.normalize],
    ['normalization_method', q.normalization_method],
    ['ref_gene', q.ref_gene],
    ['remove_outliers', q.remove_outliers],
    ['outlier_threshold', q.outlier_threshold],
    ['batchcorrect', q.batchcorrect],
    // bec_param BEFORE bec_method, deliberately. The batch-correction comparison
    // reactive depends on bec_param (and bec_full_features) but NOT on bec_method,
    // and it resets bec_method to its first choice, "ComBat"
    // (upload_module_normalization.R:251-257). Setting the method first and the
    // parameter second therefore throws the method away every time.
    ['bec_param', q.bec_param],
    ['bec_method', q.bec_method],
  ];
  for (const [id, v] of qcOrder) {
    if (v === null || v === undefined) continue;
    await setInput(page, `upload-checkqc-${id}`, v);
  }
  await settle(page, 2500);

  // bec_method must be re-asserted AFTER the batch-correction comparison has run,
  // because that comparison resets the select to its first choice (ComBat).
  if (q.batchcorrect && q.bec_method) {
    if (q.bec_param) { await setInput(page, 'upload-checkqc-bec_param', q.bec_param); await settle(page, 3000); }
    await setInputSticky(page, 'upload-checkqc-bec_method', q.bec_method);
  }
  log('QC settings applied');

  await nextStep(page, null, 'step 4 (QC)');

  // --- step 5: compute panel -------------------------------------------
  await page.waitForSelector('#upload-compute-selected_name', { state: 'attached', timeout: 120000 });
  await settle(page);

  const c = sc.compute;
  await setInput(page, 'upload-compute-selected_name', `cmp_${scenarioId}`);
  await setInput(page, 'upload-compute-selected_description', sc.description);
  await setInput(page, 'upload-compute-gene_methods', c.gene_methods);
  await setInput(page, 'upload-compute-gset_methods', c.gset_methods);
  await setInput(page, 'upload-compute-do_extra', c.do_extra);
  await setInput(page, 'upload-compute-extra_methods', c.extra_methods);
  await setInput(page, 'upload-compute-filter_methods', c.filter_methods);
  await setInput(page, 'upload-compute-dotimeseries', c.dotimeseries);
  await setInput(page, 'upload-compute-exclude_void', c.exclude_void);
  await settle(page, 2500);

  // Re-assert the QC settings: a server-side update() between step 4 and now can
  // silently revert a value set with Shiny.setInputValue. Gate 1 would catch it,
  // but re-asserting turns a confusing failure into a non-event.
  for (const [id, v] of qcOrder) {
    if (v === null || v === undefined) continue;
    await setInput(page, `upload-checkqc-${id}`, v);
  }
  await settle(page, 2000);
  log('compute settings applied');

  // --- trigger compute --------------------------------------------------
  const before = snapshotRawDirs();
  const trigger = async () => {
    await clearOverlays(page);
    // Re-assert every QC input immediately before firing. The batch-correction
    // comparison is slow and resets bec_method to ComBat whenever it finishes, which
    // can land AFTER an earlier verification passed -- so the only safe moment is
    // right before the compute observer reads compute_settings. settle() in between
    // lets the module's bc_method reactive propagate before the observer runs.
    for (const [id, v] of qcOrder) {
      if (v === null || v === undefined) continue;
      // Skip the batch-correction inputs. Re-asserting bec_param invalidates the
      // batch-correction comparison, which re-runs (slowly) and resets bec_method to
      // ComBat -- so this very re-assert was destroying the method it was meant to
      // protect. They are set once, verified stable, on step 4.
      if (id === 'bec_param' || id === 'bec_method') continue;
      await setInput(page, `upload-checkqc-${id}`, v);
    }
    await settle(page, 2000);
    await setInput(page, 'upload-upload_wizard', 'wizard_finished');
  };
  await trigger();
  log('compute triggered, polling for params.RData ...');

  rawDir = await captureParams(before, trigger);
  if (!rawDir) {
    // most likely cause: sv_required validation, or probetype detection failed
    const alert = await page.locator('.swal-title, .swal-text').allTextContents().catch(() => []);
    throw new Error(`params.RData never appeared. page alerts: ${JSON.stringify(alert)}`);
  }
  log(`captured from ${rawDir}`);

  for (const f of ['params.RData', 'CHECKS_OUTPUT', 'processx-error.log', 'processx-output.log']) {
    const src = join(rawDir, f);
    if (existsSync(src)) copyFileSync(src, join(outDir, f));
  }
  for (const f of ['counts.csv', 'samples.csv', 'contrasts.csv']) {
    copyFileSync(join(fx.dir, f), join(outDir, f));
  }
  log(`wrote ${outDir}`);
} catch (err) {
  log('FAILED:', err.message.split('\n')[0]);
  // Dump what was actually on screen -- guessing from a bare timeout is slow.
  const diag = await page.evaluate(() => ({
    openModals: [...document.querySelectorAll('.modal.show')].map((m) => m.id || '(no id)'),
    swal: [...document.querySelectorAll('.swal-modal')].map((s) => (s.textContent || '').trim().slice(0, 120)),
    swalButtons: [...document.querySelectorAll('.swal-button')].map((b) => b.textContent.trim()),
    backdrops: document.querySelectorAll('.modal-backdrop').length,
    nextButtons: [...document.querySelectorAll('.wizard-btn.btn.next')].map((b) => ({
      visible: !!b.offsetParent, disabled: b.disabled, cls: b.className,
    })),
    activeStep: document.querySelector('.wizard-step.active')?.dataset?.title || null,
  })).catch((e) => ({ diagError: e.message }));
  log('page state:', JSON.stringify(diag));
  await page.screenshot({ path: join(outDir, 'failure.png'), fullPage: true }).catch(() => {});
  await browser.close();
  process.exit(1);
}

await browser.close();
log('done');
