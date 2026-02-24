// set current year
document.getElementById("year").textContent = new Date().getFullYear();

// theme toggle
const root = document.documentElement;
const toggle = document.getElementById("themeToggle");

function setTheme(theme) {
  root.setAttribute("data-theme", theme);
  toggle.textContent = theme === "dark" ? "🌙" : "☀️";
  localStorage.setItem("theme", theme);
}

setTheme(localStorage.getItem("theme") || "dark");

toggle.addEventListener("click", () => {
  const next = root.getAttribute("data-theme") === "dark" ? "light" : "dark";
  setTheme(next);
});