// Update the "pending received datasets" indicators in the navbar:
//  - numeric badge on the "Sharing" sub-item
//  - small red dot on the parent "Datasets" dropdown toggle
// Called from loading_module_received.R via shinyjs::runjs("updateSharedBadges(n)").
function updateSharedBadges(n) {
  var sub = document.querySelector('a[data-target="sharing-tab"]');
  if (sub) {
    var oldSub = sub.querySelector('.shared-pending-badge');
    if (oldSub) oldSub.remove();
    if (n > 0) {
      var b = document.createElement('span');
      b.className = 'shared-pending-badge';
      b.textContent = n;
      sub.appendChild(b);
    }
  }
  document.querySelectorAll('.nav-link.dropdown-toggle').forEach(function (el) {
    if (el.textContent.trim().indexOf('Datasets') !== 0) return;
    var oldDot = el.querySelector('.shared-pending-dot');
    if (oldDot) oldDot.remove();
    if (n > 0) {
      var d = document.createElement('span');
      d.className = 'shared-pending-dot';
      el.appendChild(d);
    }
  });
}
