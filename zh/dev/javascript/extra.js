// Copy from Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

// Block sitemap.xml fetches triggered by Weglot's hreflang tags detected by MkDocs Material
(() => {
    const EMPTY_SITEMAP = `<?xml version="1.0" encoding="UTF-8"?><urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9"></urlset>`;
  
    const originalFetch = window.fetch;
    window.fetch = function (url, options) {
      if (typeof url === "string" && url.includes("/sitemap.xml")) {
        return Promise.resolve(
          new Response(EMPTY_SITEMAP, { status: 200, headers: { "Content-Type": "application/xml" } }),
        );
      }
      return originalFetch.apply(this, arguments);
    };
  
    const originalXHROpen = XMLHttpRequest.prototype.open;
    XMLHttpRequest.prototype.open = function (method, url) {
      if (typeof url === "string" && url.includes("/sitemap.xml")) {
        this._blockRequest = true;
      }
      return originalXHROpen.apply(this, arguments);
    };
  
    const originalXHRSend = XMLHttpRequest.prototype.send;
    XMLHttpRequest.prototype.send = function () {
      if (this._blockRequest) {
        Object.defineProperty(this, "status", { value: 200 });
        Object.defineProperty(this, "responseText", { value: EMPTY_SITEMAP });
        Object.defineProperty(this, "response", { value: EMPTY_SITEMAP });
        Object.defineProperty(this, "responseXML", {
          value: new DOMParser().parseFromString(EMPTY_SITEMAP, "application/xml"),
        });
        this.dispatchEvent(new Event("load"));
        return;
      }
      return originalXHRSend.apply(this, arguments);
    };
  })();
  
  // Apply theme colors based on dark/light mode
  const applyTheme = (isDark) => {
    document.body.setAttribute("data-md-color-scheme", isDark ? "slate" : "default");
    document.body.setAttribute("data-md-color-primary", isDark ? "black" : "indigo");
  };
  
  // Sync widget theme with Material theme
  const syncWidgetTheme = () => {
    const isDark = document.body.getAttribute("data-md-color-scheme") === "slate";
    document.documentElement.setAttribute("data-theme", isDark ? "dark" : "light");
  };
  
  // Check and apply appropriate theme based on system/user preference
  const checkTheme = () => {
    const palette = JSON.parse(localStorage.getItem(".__palette") || "{}");
    if (palette.index === 0) {
      applyTheme(window.matchMedia("(prefers-color-scheme: dark)").matches);
      syncWidgetTheme();
    }
  };
  
  // Initialize theme handling on page load
  document.addEventListener("DOMContentLoaded", () => {
    checkTheme();
    syncWidgetTheme();
  
    // Watch for system theme changes
    window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", checkTheme);
  
    // Watch for theme toggle changes
    document.getElementById("__palette_1")?.addEventListener("change", (e) => {
      if (e.target.checked) setTimeout(checkTheme);
    });
  
    // Watch for Material theme changes and sync to widget
    new MutationObserver(syncWidgetTheme).observe(document.body, {
      attributes: true,
      attributeFilter: ["data-md-color-scheme"],
    });
  });
  
  
  (() => {
    function fixLanguageLinks() {
      const path = location.pathname; 
      const links = document.querySelectorAll(".md-select__link[hreflang]");
      if (!links.length) return;
  
      const suffix = location.search + location.hash;
  
      links.forEach((link) => {
        const lang = link.getAttribute("hreflang");
        let newPath = path;

        if (lang === "zh") {
          if (!path.includes("/CPhasing/zh/")) {
            newPath = path.replace("/CPhasing/", "/CPhasing/zh/");
          }
        } else {
          if (path.includes("/CPhasing/zh/")) {
            newPath = path.replace("/CPhasing/zh/", "/CPhasing/");
          }
        }
        link.href = newPath + suffix;
      });
    }
  
    fixLanguageLinks();
  
    if (typeof document$ !== "undefined") {
      document$.subscribe(() => setTimeout(fixLanguageLinks, 50));
    }
  })();