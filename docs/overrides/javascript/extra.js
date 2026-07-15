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
  
  // Preserve the current page, query parameters, and hash when switching
// between documentation languages and versions.
(() => {
  const VERSION_PATTERN =
    /^(?:v?\d+(?:\.\d+){1,3}(?:[-+][a-z0-9.-]+)?|latest|dev|stable)$/i;

  function isVersionSegment(segment) {
    return VERSION_PATTERN.test(segment || "");
  }

  function getRootPath(pathname = window.location.pathname) {
    return pathname.startsWith("/CPhasing/")
      ? "/CPhasing/"
      : "/";
  }

  /**
   * Parse a documentation URL.
   *
   * Examples:
   * /CPhasing/latest/installation/
   *   -> { version: "latest", language: "en", pagePath: "installation/" }
   *
   * /CPhasing/v0.3.1/zh/installation/
   *   -> { version: "v0.3.1", language: "zh", pagePath: "installation/" }
   */
  function parseDocumentationPath(pathname, rootPath) {
    let relativePath = pathname;

    if (relativePath.startsWith(rootPath)) {
      relativePath = relativePath.slice(rootPath.length);
    } else {
      relativePath = relativePath.replace(/^\/+/, "");
    }

    const segments = relativePath.split("/").filter(Boolean);

    let version = "";
    let language = "en";

    if (segments.length > 0 && isVersionSegment(segments[0])) {
      version = segments.shift();
    }

    if (segments[0] === "zh") {
      language = "zh";
      segments.shift();
    }

    const pagePath = segments.join("/");

    return {
      version,
      language,
      pagePath,
    };
  }

  function extractTargetVersion(link, rootPath) {
    const rawHref = link.getAttribute("href");

    if (!rawHref) {
      return "";
    }

    const targetUrl = new URL(rawHref, window.location.href);
    const targetInfo = parseDocumentationPath(
      targetUrl.pathname,
      rootPath,
    );

    if (targetInfo.version) {
      return targetInfo.version;
    }

    // Fallback: try to determine the version from the link text.
    const linkText = link.textContent.trim();

    if (isVersionSegment(linkText)) {
      return linkText;
    }

    return "";
  }

  function buildDocumentationUrl({
    rootPath,
    version,
    language,
    pagePath,
    suffix,
  }) {
    const segments = [];

    if (version) {
      segments.push(version);
    }

    if (language === "zh") {
      segments.push("zh");
    }

    if (pagePath) {
      segments.push(...pagePath.split("/").filter(Boolean));
    }

    let targetPath = rootPath + segments.join("/");

    // MkDocs normally uses directory-style URLs.
    if (!targetPath.endsWith("/")) {
      targetPath += "/";
    }

    targetPath = targetPath.replace(/\/{2,}/g, "/");

    return targetPath + suffix;
  }

  function fixLanguageLinks(currentInfo, rootPath, suffix) {
    const languageLinks =
      document.querySelectorAll(".md-select__link[hreflang]");

    languageLinks.forEach((link) => {
      const hreflang = link.getAttribute("hreflang");
      const targetLanguage =
        hreflang && hreflang.toLowerCase().startsWith("zh")
          ? "zh"
          : "en";

      link.href = buildDocumentationUrl({
        rootPath,
        version: currentInfo.version,
        language: targetLanguage,
        pagePath: currentInfo.pagePath,
        suffix,
      });
    });
  }

  function fixVersionLinks(currentInfo, rootPath, suffix) {
    const versionContainer = document.querySelector(
      '[data-md-component="version"]',
    );

    if (!versionContainer) {
      return;
    }

    const versionLinks = versionContainer.querySelectorAll("a[href]");

    versionLinks.forEach((link) => {
      const targetVersion = extractTargetVersion(link, rootPath);

      if (!targetVersion) {
        return;
      }

      link.href = buildDocumentationUrl({
        rootPath,
        version: targetVersion,
        language: currentInfo.language,
        pagePath: currentInfo.pagePath,
        suffix,
      });
    });
  }

  function fixLanguageAndVersionLinks() {
    const rootPath = getRootPath();
    const currentInfo = parseDocumentationPath(
      window.location.pathname,
      rootPath,
    );
    const suffix = window.location.search + window.location.hash;

    fixLanguageLinks(currentInfo, rootPath, suffix);
    fixVersionLinks(currentInfo, rootPath, suffix);
  }

  function scheduleLinkFix() {
    // Material may render selectors asynchronously.
    requestAnimationFrame(() => {
      fixLanguageAndVersionLinks();

      setTimeout(fixLanguageAndVersionLinks, 50);
      setTimeout(fixLanguageAndVersionLinks, 200);
    });
  }

  // Normal initial page load.
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", scheduleLinkFix);
  } else {
    scheduleLinkFix();
  }

  // MkDocs Material instant navigation.
  if (typeof document$ !== "undefined") {
    document$.subscribe(scheduleLinkFix);
  }

  // Version and language menus can be rendered or updated dynamically.
  const observer = new MutationObserver((mutations) => {
    const selectorChanged = mutations.some((mutation) =>
      Array.from(mutation.addedNodes).some(
        (node) =>
          node.nodeType === Node.ELEMENT_NODE &&
          (
            node.matches?.('[data-md-component="version"], .md-select') ||
            node.querySelector?.(
              '[data-md-component="version"], .md-select',
            )
          ),
      ),
    );

    if (selectorChanged) {
      scheduleLinkFix();
    }
  });

  observer.observe(document.documentElement, {
    childList: true,
    subtree: true,
  });
})();
  