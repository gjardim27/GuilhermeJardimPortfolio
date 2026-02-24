// Dark mode toggle with localStorage persistence
(function () {
  var toggle = document.getElementById('themeToggle');
  var icon = toggle && toggle.querySelector('.toggle-icon');
  var root = document.documentElement;

  // Restore saved theme or default to light
  var saved = localStorage.getItem('theme');
  if (saved === 'dark' || saved === 'light') {
    root.setAttribute('data-theme', saved);
  }
  updateIcon();

  if (toggle) {
    toggle.addEventListener('click', function () {
      var current = root.getAttribute('data-theme');
      var next = current === 'dark' ? 'light' : 'dark';
      root.setAttribute('data-theme', next);
      localStorage.setItem('theme', next);
      updateIcon();
    });
  }

  function updateIcon() {
    if (!icon) return;
    var theme = root.getAttribute('data-theme');
    icon.textContent = theme === 'dark' ? '\u2600' : '\u263D';
  }
})();
