document.addEventListener("DOMContentLoaded", async () => {
  const container = document.querySelector("#site-menu");

  if (!container) {
    console.error("Navigation container #site-menu was not found.");
    return;
  }

  try {
    const response = await fetch("components/menu.html");

    if (!response.ok) {
      throw new Error(
        `Could not load components/menu.html: HTTP ${response.status}`
      );
    }

    container.innerHTML = await response.text();

    const normalizePath = (pathname) => {
      let path = pathname.replace(/\/+$/, "");
      path = path.replace(/\/index\.html$/i, "");
      return path || "/";
    };

    const currentUrl = new URL(window.location.href);
    const currentPath = normalizePath(currentUrl.pathname);

    for (const link of container.querySelectorAll("a[href]")) {
      const linkUrl = new URL(link.href, document.baseURI);

      if (linkUrl.origin !== currentUrl.origin) {
        continue;
      }

      if (normalizePath(linkUrl.pathname) === currentPath) {
        link.setAttribute("aria-current", "page");
      }
    }
  } catch (error) {
    console.error("Navigation menu failed to load:", error);
    container.textContent = "Navigation menu could not be loaded.";
  }
});
