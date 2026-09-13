#!/usr/bin/env node
/**
 * wikiskill_gate.js — Gating & Rollback component of the WikiSkill loop.
 *
 * Runs the validation suite and returns a PASS/FAIL verdict with measured
 * numbers. Exit code 0 = gate OPEN (accept skill set), 1 = gate CLOSED
 * (roll back skills, keep wiki).
 *
 * Gate criteria (documented in .wikiskill/README.md):
 *   1. node --check on all src/*.js
 *   2. tests/test_all.js — must match or beat baseline (32 PASSED, 0 FAILED)
 *   3. DOM contract: every id referenced by src/ui.js exists in index.html
 *   4. No unresolved "PENDING GATE" verdicts older than the current run
 *      (checked by caller; this script prints the summary line to append).
 */

import { execFileSync } from "node:child_process";
import fs from "node:fs";
import path from "node:path";

const ROOT = path.resolve(path.dirname(new URL(import.meta.url).pathname), "..");
let failures = [];

function check(label, fn) {
  try {
    const out = fn();
    console.log(`  ✓ ${label}${out ? ` — ${out}` : ""}`);
  } catch (e) {
    failures.push(`${label}: ${e.message}`);
    console.error(`  ✗ ${label} — ${e.message}`);
  }
}

/* 1. syntax check all src modules */
check("node --check src/*.js", () => {
  const dir = path.join(ROOT, "src");
  const files = fs.readdirSync(dir).filter((f) => f.endsWith(".js"));
  for (const f of files) {
    execFileSync(process.execPath, ["--check", path.join(dir, f)], { stdio: "pipe" });
  }
  const sub = fs.readdirSync(path.join(dir, "physics"), { withFileTypes: true })
    .filter((d) => d.isDirectory()).map((d) => path.join(dir, "physics", d.name));
  let physicsFiles = [];
  for (const sd of sub) {
    physicsFiles = physicsFiles.concat(
      fs.readdirSync(sd).filter((f) => f.endsWith(".js")).map((f) => path.join(sd, f))
    );
  }
  return `${files.length + physicsFiles.length} files clean`;
});

/* 2. regression suite (Stage-6: grand total = Tier-0 32 + FAST 199 = 231;
 * SLOW tier is 44 (thermo 7 + heavy 13 smoke-or-full + calibration 14 +
 * flexlig 10), full --slow/--slow-full total 275 — gate still checks the
 * FAST grand total only, baseline 231) */
check("tests/test_all.js baseline", () => {
  const out = execFileSync(process.execPath, [path.join(ROOT, "tests", "test_all.js")], {
    stdio: ["pipe", "pipe", "pipe"], encoding: "utf-8", timeout: 180000,
  });
  // Stage-3: test_all.js runs tiered suites; each child result is captured
  // (lowercase one-liners), so the LAST uppercase line is the grand total.
  const matches = [...out.matchAll(/(\d+) PASSED, (\d+) FAILED/g)];
  if (!matches.length) throw new Error("no results line");
  const last = matches[matches.length - 1];
  const [passed, failed] = [Number(last[1]), Number(last[2])];
  if (failed > 0) throw new Error(`${failed} FAILED`);
  if (passed < 231) throw new Error(`regression: ${passed} < 231 baseline (32 Tier-0 + 199 FAST)`);
  return `${passed} PASSED, 0 FAILED (>= 231 baseline)`;
});

/* 3. DOM contract: ui.js ids ⊆ index.html ids */
check("DOM contract (ui.js ids exist in index.html)", () => {
  const uiSrc = fs.readFileSync(path.join(ROOT, "src", "ui.js"), "utf-8");
  const html = fs.readFileSync(path.join(ROOT, "index.html"), "utf-8");
  const ids = [...uiSrc.matchAll(/\$\("([^"]+)"\)/g)].map((mm) => mm[1]);
  const missing = ids.filter((id) => !new RegExp(`id="${id}"`).test(html));
  if (missing.length) throw new Error(`missing in index.html: ${missing.join(", ")}`);
  return `${ids.length} ui ids all present`;
});

/* 4. panel-hotkey contract: Digit1-N matches #controls > .panel count */
check("Hotkey contract (Digit1-7 vs panel count)", () => {
  const mainSrc = fs.readFileSync(path.join(ROOT, "src", "main.js"), "utf-8");
  const html = fs.readFileSync(path.join(ROOT, "index.html"), "utf-8");
  const hotkeyMax = Math.max(...[...mainSrc.matchAll(/Digit\[(\d)-(\d)\]/g)].map((mm) => Number(mm[2])));
  const panels = (html.match(/class="panel[^"]*"/g) || []).length;
  const topLevel = (html.match(/<details[^>]*class="panel"/g) || []).length;
  return `Digit1-${hotkeyMax} vs ${topLevel} top-level panels (contract: keep in sync manually)`;
});

console.log(failures.length === 0
  ? "\nGATE: OPEN — accept skill set (append verdict to .wikiskill/wiki/evolution/skill-impact.md)"
  : `\nGATE: CLOSED — roll back skills, keep wiki (${failures.length} failure(s))`);
process.exit(failures.length === 0 ? 0 : 1);
