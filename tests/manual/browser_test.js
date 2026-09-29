/**
 * browser_test.js — Full Browser Automation & Feature Verification with Screenshots
 */

import { firefox } from "playwright";
import http from "http";
import fs from "fs";
import path from "path";

// Simple static file server
function startServer(port = 8799) {
  const mimeTypes = {
    ".html": "text/html",
    ".js": "application/javascript",
    ".css": "text/css",
    ".json": "application/json",
    ".pdb": "text/plain",
    ".mol2": "text/plain",
    ".svg": "image/svg+xml",
    ".ico": "image/x-icon",
  };

  const server = http.createServer((req, res) => {
    let reqPath = req.url.split("?")[0];
    if (reqPath === "/") reqPath = "/index.html";
    const filePath = path.join(process.cwd(), reqPath);

    fs.readFile(filePath, (err, data) => {
      if (err) {
        res.writeHead(404, { "Content-Type": "text/plain" });
        res.end("404 Not Found: " + reqPath);
        return;
      }
      const ext = path.extname(filePath);
      res.writeHead(200, {
        "Content-Type": mimeTypes[ext] || "application/octet-stream",
        "Access-Control-Allow-Origin": "*",
      });
      res.end(data);
    });
  });

  return new Promise((resolve) => {
    server.listen(port, "127.0.0.1", () => {
      console.log(`Test server running at http://127.0.0.1:${port}/`);
      resolve(server);
    });
  });
}

async function runBrowserTest() {
  const server = await startServer(8799);
  fs.mkdirSync("screenshots", { recursive: true });

  const browser = await firefox.launch({ headless: true });
  const context = await browser.newContext({ viewport: { width: 1280, height: 800 } });
  const page = await context.newPage();

  const consoleLogs = [];
  const errors = [];

  page.on("console", (msg) => {
    consoleLogs.push(`[${msg.type()}] ${msg.text()}`);
    if (msg.type() === "error") {
      console.error(`Browser Error: ${msg.text()}`);
    }
  });

  page.on("pageerror", (err) => {
    errors.push(err.message);
    console.error(`Page Error: ${err.message}`);
  });

  console.log("\n[1] Navigating to simulation app...");
  await page.goto("http://127.0.0.1:8799/index.html", { waitUntil: "networkidle" });
  await page.screenshot({ path: "screenshots/01_initial_load.png" });

  // -------------------------------------------------------------
  // Test Preset Loading: 4W52
  // -------------------------------------------------------------
  console.log("\n[2] Testing Preset Loading (4W52)...");
  await page.click("a[data-ex='4W52']");
  await page.waitForTimeout(1000);

  let hudText = await page.textContent("#hud");
  console.log("HUD after 4W52 load:", hudText);
  await page.screenshot({ path: "screenshots/02_4w52_loaded.png" });

  // -------------------------------------------------------------
  // Test Running Langevin Simulation
  // -------------------------------------------------------------
  console.log("\n[3] Testing Play / Simulation Step...");
  await page.click("#playBtn");
  await page.waitForTimeout(1500);
  hudText = await page.textContent("#hud");
  console.log("HUD after running simulation:", hudText);
  await page.screenshot({ path: "screenshots/03_simulation_running.png" });
  await page.click("#playBtn"); // pause

  // -------------------------------------------------------------
  // Test Heavy-Atom Mode
  // -------------------------------------------------------------
  console.log("\n[4] Testing All-Atom Heavy Mode...");
  await page.selectOption("#modelMode", "heavy");
  await page.click("#buildBtn");
  await page.waitForTimeout(1000);

  const selSummary = await page.textContent("#selSummary");
  console.log("Heavy Mode Summary:", selSummary);
  await page.click("#playBtn");
  await page.waitForTimeout(1500);
  hudText = await page.textContent("#hud");
  console.log("HUD in Heavy Mode:", hudText);
  await page.screenshot({ path: "screenshots/04_heavy_mode_running.png" });
  await page.click("#playBtn"); // pause

  // -------------------------------------------------------------
  // Test File Loader with 4w52.pdb & benzene.mol2
  // -------------------------------------------------------------
  console.log("\n[5] Testing File Loader with local PDB & MOL2...");
  await page.setInputFiles("#fileInput", "4w52.pdb");
  await page.waitForTimeout(1000);
  await page.setInputFiles("#mol2File", "benzene.mol2");
  await page.waitForTimeout(1000);
  await page.screenshot({ path: "screenshots/05_file_loader_mol2.png" });

  // -------------------------------------------------------------
  // Test Ligand Pocket Placement
  // NOTE (wiki P4/P5): "Ligand Placement" panel is collapsed by default
  // after the decluttering redesign — open it before clicking inside.
  // -------------------------------------------------------------
  console.log("\n[6] Testing Auto Pocket Placement...");
  const ligPanel = page.locator("details.panel", { has: page.locator("#placePocketBtn") });
  if (!(await ligPanel.getAttribute("open"))) await ligPanel.locator("summary").click();
  await page.locator("#placePocketBtn").scrollIntoViewIfNeeded();
  await page.click("#placePocketBtn");
  await page.waitForTimeout(800);
  const placeInfo = await page.textContent("#ligPlaceInfo");
  console.log("Placement info:", placeInfo);
  await page.screenshot({ path: "screenshots/06_pocket_placed.png" });

  // -------------------------------------------------------------
  // Test Settings Modal
  // -------------------------------------------------------------
  console.log("\n[7] Testing Settings Modal...");
  await page.click("#settingsBtn");
  await page.waitForTimeout(500);
  await page.screenshot({ path: "screenshots/07_settings_modal.png" });
  await page.click("#closeSettingsBtn");

  // -------------------------------------------------------------
  // Test Chemical Network Model Controls
  // NOTE (wiki P4/P5): network panel is collapsed by default — open it.
  // -------------------------------------------------------------
  console.log("\n[8] Testing Chemical Network Controls...");
  const netPanel = page.locator("details.panel#networkPanelDetails");
  if (!(await netPanel.getAttribute("open"))) await netPanel.locator("summary").first().click();
  await page.locator("#netStepBtn").scrollIntoViewIfNeeded();
  await page.click("#netStepBtn");
  await page.waitForTimeout(300);
  await page.locator("#netMultiStepBtn").scrollIntoViewIfNeeded();
  await page.click("#netMultiStepBtn");
  await page.waitForTimeout(1200);
  await page.screenshot({ path: "screenshots/08_chemical_network.png" });

  console.log("\n=== BROWSER TEST COMPLETE ===");
  console.log("Console Errors:", errors.length === 0 ? "NONE (Clean)" : errors);

  await browser.close();
  server.close();

  if (errors.length > 0) {
    console.error("Browser tests failed with errors:", errors);
    process.exit(1);
  } else {
    console.log("ALL BROWSER TESTS PASSED! Screenshots saved to ./screenshots/");
  }
}

runBrowserTest().catch((e) => {
  console.error("Test failure:", e);
  process.exit(1);
});
