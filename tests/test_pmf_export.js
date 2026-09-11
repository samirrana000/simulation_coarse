/**
 * test_pmf_export.js — D40 PMF export provenance
 * Exports CSV and checks header contains "# T=300, gamma=6, hills=..., V0=1660.54"
 */

import { Funnel } from "../src/funnel.js";
import { pmfCsv } from "../src/analysis.js";

function assert(cond, msg) {
  if (!cond) { console.error("✗ FAIL: " + msg); process.exit(1); }
  console.log("✓ " + msg);
}

function makeRef(nProt, nLig) {
  const n = nProt + nLig;
  const ref = new Float64Array(3 * n);
  for (let i = 0; i < nProt; i++) { ref[3*i]= (Math.random()-0.5)*2; ref[3*i+1]=0; ref[3*i+2]=0; }
  for (let a=0;a<nLig;a++){ const idx=nProt+a; ref[3*idx]=2; ref[3*idx+1]=0; ref[3*idx+2]=0; }
  return ref;
}

function run() {
  console.log("=== test_pmf_export.js — D40 provenance ===");
  const nProt=6, nLig=2;
  const ref = makeRef(nProt, nLig);
  const funnel = new Funnel({ nProt, n: nProt+nLig, ref, biasFactor:6, rMax:24, bins:96 });
  funnel.T = 300; // ensure header check
  funnel.biasFactor = 6;
  for (let i=0;i<5;i++) funnel.deposit(3+i*0.4);

  // via Funnel.exportPMF / getPMFcsv (src/funnel.js:620)
  const csv1 = funnel.exportPMF();
  console.log("  Funnel.exportPMF first line:", csv1.split("\n")[0]);
  assert(csv1.includes("# T=300"), `Funnel CSV includes "# T=300" (got ${csv1.split("\n")[0]})`);
  assert(csv1.includes("gamma=6"), `Funnel CSV includes gamma=6`);
  assert(csv1.includes("hills="), `Funnel CSV includes hills=`);
  assert(csv1.includes("V0=1660.54"), `Funnel CSV includes V0=1660.54 (got ${csv1.split("\n")[0]})`);

  // also via getPMFcsv alias
  const csvAlias = funnel.getPMFcsv();
  assert(csvAlias.includes("# T=300"), "getPMFcsv includes T");

  // via analysis.js pmfCsv
  const csv2 = pmfCsv(funnel);
  console.log("  analysis pmfCsv first line:", csv2.split("\n")[0]);
  assert(csv2.includes("# T=300"), `analysis pmfCsv includes "# T=300"`);
  assert(csv2.includes("gamma=6"), `analysis pmfCsv includes gamma=6`);
  assert(csv2.includes("hills="), `analysis pmfCsv includes hills`);
  assert(csv2.includes("V0=1660.54"), `analysis pmfCsv includes V0=1660.54`);

  // Check header format exactly "# T=300, gamma=6, hills=..., V0=1660.54" pattern
  const headerPat = /# T=300, gamma=6, hills=\d+, V0=1660\.54/;
  assert(headerPat.test(csv1), `Funnel header matches ${headerPat}`);
  assert(headerPat.test(csv2), `analysis header matches ${headerPat}`);

  // Ensure r_Ang,pmf header follows
  assert(csv2.includes("r_Ang,pMF_kcal_per_mol") || csv2.includes("r_Ang"), "CSV contains r_Ang header");

  console.log("\nPASS test_pmf_export.js — CSV header provenance OK");
}

run();
