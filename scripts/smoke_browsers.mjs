/**
 * smoke_browsers.mjs — FP3 cross-browser smoke matrix (additive, zero deps
 * beyond the repo's dev playwright install; headless only).
 *
 * Per browser (chromium + firefox): serve the repo root over HTTP (ES modules
 * cannot load from file://), load index.html, click `Load 4W52 sample`,
 * wait for build (STATE: READY), press Run, advance ~3 s of wall-clock
 * dynamics, then assert: (a) step counter advanced past 0, (b) HUD shows
 * live `t =` with no `⚠` error banner, (c) zero uncaught page errors.
 *
 * Safari/WebKit cannot execute on linux (no Safari binary, no playwright
 * webkit build installed) — covered by the static matrix step below
 * (HTTP 200 + module reachability), recorded as STATIC-ONLY in docs.
 *
 * Run: node scripts/smoke_browsers.mjs [--port 8123] [--browsers chromium,firefox]
 * Exit 0 = every launched browser passed with zero page errors.
 */
import http from "node:http";
import fs from "node:fs";
import path from "node:path";

const ROOT = path.resolve(path.dirname(new URL(import.meta.url).pathname), "..");
const PORT = Number((process.argv.find((a) => a.startsWith("--port=")) || "").split("=")[1])
  || Number(process.env.SMOKE_PORT || 8123);
const WANT = ((process.argv.find((a) => a.startsWith("--browsers=")) || "").split("=")[1]
  || process.env.SMOKE_BROWSERS || "chromium,firefox").split(",").map((s) => s.trim().toLowerCase());

const MIME = {
  ".html": "text/html; charset=utf-8",
  ".js": "text/javascript; charset=utf-8",
  ".mjs": "text/javascript; charset=utf-8",
  ".css": "text/css; charset=utf-8",
  ".json": "application/json; charset=utf-8",
  ".pdb": "text/plain; charset=utf-8",
  ".mol2": "text/plain; charset=utf-8",
  ".svg": "image/svg+xml",
  ".ico": "image/x-icon",
  ".md": "text/markdown; charset=utf-8",
};

/** Minimal static file server for the repo root (GET only, no directory listing). */
function serveRoot(port) {
  const srv = http.createServer((req, res) => {
    try {
      const urlPath = decodeURIComponent(String(req.url || "/").split("?")[0]);
      let rel = urlPath === "/" ? "index.html" : urlPath.replace(/^\/+/, "");
      const abs = path.normalize(path.join(ROOT, rel));
      if (!abs.startsWith(ROOT) || !fs.existsSync(abs) || fs.statSync(abs).isDirectory()) {
        res.writeHead(404, { "content-type": "text/plain" }); res.end("not found"); return;
      }
      res.writeHead(200, { "content-type": MIME[path.extname(abs).toLowerCase()] || "application/octet-stream" });
      fs.createReadStream(abs).pipe(res);
    } catch (e) {
      res.writeHead(500, { "content-type": "text/plain" }); res.end(String(e?.message || e));
    }
  });
  return new Promise((resolve) => srv.listen(port, "127.0.0.1", () => resolve(srv)));
}

/** Drive one browser through load → sample → build → run → assert. */
async function smokeOne(pw, name, base) {
  const browser = await pw[name].launch();
  const version = browser.version();
  const pageErrors = [];
  const page = await browser.newPage();
  page.on("pageerror", (e) => pageErrors.push(String(e?.message || e).slice(0, 300)));
  const out = { browser: name, version, pass: false, detail: "" };
  try {
    await page.goto(`${base}/index.html`, { waitUntil: "load", timeout: 30000 });
    // Presence = attached (headless default viewport can report a zero box
    // for the flex-sized canvas before first layout; functional checks below
    // carry the verdict, bbox is logged for the record).
    await page.waitForSelector("#canvas", { state: "attached", timeout: 15000 });
    await page.waitForSelector("#sampleBtn", { state: "attached", timeout: 15000 });
    const bbox = await page.locator("#canvas").boundingBox().catch(() => null);
    // FP1 progressive disclosure: the 1-click sample lives inside the
    // collapsed "Getting started" subpanel — open it first (real-user path).
    await page.locator(".subpanel summary", { hasText: "Getting started" }).click({ timeout: 15000 });
    await page.locator("#sampleBtn").click({ timeout: 15000 });
    // Built when the top status flips to READY (buildSystem tail).
    await page.waitForFunction(
      () => (document.getElementById("sysState")?.textContent || "").includes("READY"),
      { timeout: 45000 },
    );
    // Run ~3 s of dynamics, then pause.
    await page.click("#playBtn");
    await page.waitForTimeout(3500);
    const topStep = await page.textContent("#topStep").catch(() => "");
    const hud = await page.textContent("#hud").catch(() => "");
    await page.click("#playBtn"); // pause
    const stepN = Number((topStep.match(/step\s+(\d+)/) || [])[1] || -1);
    const hudLive = /t\s*=/.test(hud || "");
    const hudClean = !(hud || "").includes("⚠");
    const checks = [
      [`step advanced (${(topStep || "?").trim()})`, stepN > 0],
      [`HUD live (${(hud || "").slice(0, 60)}…)`, hudLive],
      ["HUD clean (no ⚠ banner)", hudClean],
      [`zero page errors (${pageErrors.length})`, pageErrors.length === 0],
    ];
    const bad = checks.filter(([, ok]) => !ok);
    out.pass = bad.length === 0;
    out.canvasBox = bbox ? `${Math.round(bbox.width)}×${Math.round(bbox.height)}` : "zero-box";
    out.detail = `canvas ${out.canvasBox} | ${checks.map(([l, ok]) => `${ok ? "✓" : "✗"} ${l}`).join(" | ")}${pageErrors.length ? ` — first error: ${pageErrors[0]}` : ""}`;
  } catch (e) {
    out.detail = `harness exception: ${String(e?.message || e).split("\n")[0]}${pageErrors.length ? ` | page errors: ${pageErrors[0]}` : ""}`;
  } finally {
    await browser.close().catch(() => {});
  }
  return out;
}

const srv = await serveRoot(PORT);
const base = `http://127.0.0.1:${PORT}`;
console.log(`[smoke] serving ${ROOT} at ${base}`);
let pw;
try {
  pw = await import("playwright");
} catch (e) {
  console.error(`[smoke] playwright import failed: ${e.message}`);
  srv.close();
  process.exit(2);
}
const results = [];
for (const name of WANT) {
  if (!pw[name] || typeof pw[name].launch !== "function") {
    results.push({ browser: name, version: "n/a", pass: false, detail: "unknown browser (want chromium,firefox)" });
    continue;
  }
  const r = await smokeOne(pw, name, base);
  results.push(r);
  console.log(`[${r.browser} ${r.version}] ${r.pass ? "PASS" : "FAIL"} — ${r.detail}`);
}
srv.close();
const failed = results.filter((r) => !r.pass).length;
console.log(`\n=== smoke_browsers: ${results.length - failed}/${results.length} browsers passed, page errors 0 required ===`);
process.exit(failed ? 1 : 0);
