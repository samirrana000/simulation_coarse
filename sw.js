/**
 * sw.js — freshness mechanism for the app's ES-module graph.
 *
 * WHAT THIS REPLACES
 * ------------------
 * Until 2026-09-29 every module edge carried a hand-typed cache-bust
 * suffix: `import { ForceField } from "./forcefield.js?v=10"`. 121 literals
 * across 42 files, plus one on index.html's module tag. Two things were
 * wrong with that scheme:
 *
 *   1. It was a discipline, not a mechanism. A module edited without also
 *      editing its neighbours' literals kept the SAME url, so the browser
 *      served the cached bytes and the developer saw "my fix did nothing" —
 *      indistinguishable from a physics bug.
 *   2. Even when the version was bumped it was one edit per module edge, and
 *      the literal (`10`) had drifted from src/version.js (`1.1.0-fp7`).
 *
 * The version bump could not be a single edit without a build step: a query
 * string on the entry point does not propagate to static relative imports,
 * so each edge would still need its own literal. (Module specifiers are
 * resolved by the parser before any code runs, so a runtime loader cannot
 * rewrite them either.)
 *
 * WHAT THIS DOES INSTEAD
 * ----------------------
 * It removes the reason the query string existed. Cache freshness is a
 * property of the HTTP response, not of the URL, and this file takes
 * control of the cache instead of the URLs:
 *
 *   - every same-origin GET (the document, all 66 modules, the WGSL
 *     sources fetched at runtime) is re-fetched from the network with
 *     `cache: "no-store"`, which bypasses the HTTP cache outright — so the
 *     10%-of-age heuristic freshness that a bare `python3 -m http.server`
 *     response earns can never apply; and
 *   - this worker holds NO cache of its own. There is no `caches.open`,
 *     no `cache.match`, no offline copy. Its only response path is a
 *     network response it just received, so there is no stored copy of the
 *     app that can go stale. tests/test_cache_contract.js asserts that
 *     absence, because a future "let me add an offline mode" edit that
 *     reintroduced a cache read would silently restore the original bug.
 *
 * Consequences, stated so they are not a surprise:
 *   - Every reload re-downloads ~1 MB of modules from the dev server. That
 *     is the price of "what you see is what is on disk".
 *   - The app no longer works offline. Serving a stale graph is a worse
 *     failure for an interactive simulator than refusing to start.
 *   - Where service workers are unavailable (an insecure origin such as
 *     http://192.168.x.x, or file://) registration is skipped and the app
 *     falls back to the plain, query-free urls in index.html/src — correct,
 *     just no longer guaranteed fresh.
 *
 * Build-free by construction: this is a static file served verbatim, with
 * no imports, no bundler, and no dependency of any kind. Its own freshness
 * is the browser's job — index.html registers it with
 * `updateViaCache: "none"` and a query-free url, so the worker's bytes are
 * revalidated on every navigation and a changed sw.js installs itself.
 */

self.addEventListener("install", () => {
  // Nothing to pre-cache on purpose. Claim the page immediately so the
  // very first visit is covered too; index.html turns the resulting
  // controllerchange into exactly one reload (at most, per tab).
  self.skipWaiting();
});

self.addEventListener("activate", (event) => {
  // Drop any Cache Storage this scope accumulated from an older build of
  // this file. Nothing reads it (see above), so deleting is safe and keeps
  // a previous revision from lingering on a developer's disk.
  event.waitUntil(
    caches.keys()
      .then((keys) => Promise.all(keys.map((key) => caches.delete(key))))
      .then(() => self.clients.claim())
  );
});

self.addEventListener("fetch", (event) => {
  const request = event.request;
  if (request.method !== "GET") return;                    // never touch writes
  let url;
  try {
    url = new URL(request.url);
  } catch {
    return;                                                // not a URL: default
  }
  if (url.origin !== self.location.origin) return;          // cross-origin: default
  event.respondWith(fetch(request, { cache: "no-store" }));
});
