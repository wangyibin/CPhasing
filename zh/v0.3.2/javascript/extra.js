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
  
  
  // Dynamically fix Language & Version switchers to preserve active page path, queries, and hash anchors
  (() => {
    // Helper to evaluate version segment patterns (supports v0.1.0, latest, dev, stable etc.)
    function isVersionSegment(str) {
      if (!str) return false;
      return /^(v?\d+\.\d+\.\d+|latest|dev|stable)$/i.test(str);
    }

    // Parse paths as [language]/[version]/[page].
    // Chinese docs are deployed by mike with "--deploy-prefix zh".
    function parseMkDocsPath(path, rootPath) {
      const relative = path.substring(rootPath.length);
      const segments = relative.split("/").filter(Boolean);

      const language = segments[0] === "zh" ? segments.shift() : "";
      const version = isVersionSegment(segments[0]) ? segments.shift() : "";
      const pagePath = segments.join("/");

      return { version, isZh: language === "zh", pagePath };
    }

    function fixLanguageAndVersionLinks() {
      const path = location.pathname;
      const suffix = location.search + location.hash;

      let rootPath = "/";
      if (path.startsWith("/CPhasing/")) {
        rootPath = "/CPhasing/";
      }

      const currentInfo = parseMkDocsPath(path, rootPath);

      // 1. Fix Language Switcher links
      const langLinks = document.querySelectorAll(".md-select__link[hreflang]");
      langLinks.forEach((link) => {
        const lang = link.getAttribute("hreflang");
        const targetSegments = [];
        if (lang === "zh") {
          targetSegments.push("zh");
        }
        if (currentInfo.version) {
          targetSegments.push(currentInfo.version);
        }
        if (currentInfo.pagePath) {
          targetSegments.push(currentInfo.pagePath);
        }
        const targetPath = rootPath + targetSegments.join("/");
        link.href = targetPath + suffix;
      });

      // 2. Fix mike Version Switcher links
      const versionContainer = document.querySelector('[data-md-component="version"]');
      if (versionContainer) {
        const versionLinks = versionContainer.querySelectorAll("a");
        versionLinks.forEach((link) => {
          // Resolve target URL of the selected version option
          const testUrl = new URL(link.href, location.origin);
          const targetInfo = parseMkDocsPath(testUrl.pathname, rootPath);

          const targetSegments = [];
          if (currentInfo.isZh) {
            targetSegments.push("zh");
          }
          if (targetInfo.version) {
            targetSegments.push(targetInfo.version);
          }
          if (currentInfo.pagePath) {
            targetSegments.push(currentInfo.pagePath);
          }
          const targetPath = rootPath + targetSegments.join("/");
          link.href = targetPath + suffix;
        });
      }
    }

    // Run on initial load
    fixLanguageAndVersionLinks();

    // Re-bind when navigated via MkDocs instant loading
    if (typeof document$ !== "undefined") {
      document$.subscribe(() => {
        setTimeout(fixLanguageAndVersionLinks, 50);
      });
    }
  })();
