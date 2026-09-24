// Update the "pending received datasets" count badge in two places:
//  - the "Library" pill in the left sidebar (#app-sidebar)
//  - the "Shared datasets" inner tab title (#load-tabs)
// Called from loading_module_received.R via shinyjs::runjs("updateSharedBadges(n)").
function updateSharedBadges(n) {
  var targets = [
    document.querySelector('#app-sidebar a[data-value="Library"]'),
    document.querySelector('#load-tabs a[data-value="sharing_tab"]')
  ];
  targets.forEach(function (el) {
    if (!el) return;
    var old = el.querySelector('.shared-pending-badge');
    if (old) old.remove();
    if (n > 0) {
      var b = document.createElement('span');
      b.className = 'shared-pending-badge';
      b.textContent = n;
      el.appendChild(b);
    }
  });
}

// Red dot on Settings > "AI Features" while AI is off (show = true from
// appsettings_server.R). Once the user opens that tab it stays hidden in this
// browser (localStorage); turning AI on hides it everywhere (server side).
var AI_DOT_SEEN_KEY = 'opg-ai-settings-seen';
function aiSettingsSeen() {
  try { return localStorage.getItem(AI_DOT_SEEN_KEY) === '1'; } catch (e) { return false; }
}
function updateAiSettingsDot(show) {
  [
    document.querySelector('#app-sidebar a[data-value="Settings"]'),
    document.querySelector('#app_settings-tabs1 a[data-value="AI Features"]')
  ].forEach(function (el) {
    if (!el) return;
    var old = el.querySelector('.shared-pending-dot');
    if (old) old.remove();
    if (show && !aiSettingsSeen()) {
      var d = document.createElement('span');
      d.className = 'shared-pending-dot';
      el.appendChild(d);
    }
  });
}
document.addEventListener('click', function (e) {
  if (!e.target.closest('#app_settings-tabs1 a[data-value="AI Features"]')) return;
  try { localStorage.setItem(AI_DOT_SEEN_KEY, '1'); } catch (err) {}
  updateAiSettingsDot(false);
});
