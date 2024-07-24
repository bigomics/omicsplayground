function addTick(buttonId) {
  var button = document.getElementById(buttonId);
  var originalText = button.innerHTML;
  button.innerHTML = originalText + '<span class=\"tick\"> ✓</span>';
  button.classList.add('show-tick');
  setTimeout(function() {
    button.classList.remove('show-tick');
    setTimeout(function() {
      button.innerHTML = originalText;
    }, 500); // Match this duration with CSS transition duration
  }, 1000);
}