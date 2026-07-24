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
