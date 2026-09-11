/**
 * test_funnel_grid.js — D35 Bias grid resolution vs sigma
 * Gaussian height w=0.4, sigma=0.3, delta=0.25 integrates to w*sqrt(2pi)*sigma within 1%
 */

function assert(cond, msg) {
  if (!cond) { console.error("✗ FAIL: " + msg); process.exit(1); }
  console.log("✓ " + msg);
}

function run() {
  console.log("=== test_funnel_grid.js — D35 grid resolution ===");
  const w = 0.4, sigma = 0.3, delta = 0.25;
  const rCenter = 12.0;
  const rMax = 24.0, bins = 96;
  // Verify delta matches rMax/bins
  assert(Math.abs(rMax / bins - delta) < 1e-9, `delta = rMax/bins = ${rMax/bins}`);

  // Integrate Gaussian over grid: sum w*exp(-(rk - rCenter)^2 /2 sigma^2) * delta
  const inv2sig2 = 1 / (2 * sigma * sigma);
  let sum = 0;
  for (let k = 0; k < bins; k++) {
    const rk = (k + 0.5) * delta;
    const dr = rk - rCenter;
    const g = w * Math.exp(-dr * dr * inv2sig2);
    sum += g * delta;
  }
  const expected = w * Math.sqrt(2 * Math.PI) * sigma;
  const relErr = Math.abs(sum - expected) / expected;
  console.log(`  w=${w}, sigma=${sigma}, delta=${delta}`);
  console.log(`  integral sum*Δr = ${sum.toFixed(6)}, expected w√2πσ = ${expected.toFixed(6)}, relErr=${(relErr*100).toFixed(3)}%`);
  assert(relErr < 0.01, `Gaussian integrates within 1% (relErr ${(relErr*100).toFixed(3)}% <1%)`);

  // Also test via Funnel grid deposit approximates integral
  console.log("  Also verifying Funnel bias grid integrates similarly...");
  // Use Funnel deposit of single hill and sum bias
  import("../src/funnel.js").then(({ Funnel }) => {
    const nProt = 5, nLig = 2;
    const n = nProt + nLig;
    const ref = new Float64Array(3 * n);
    for (let i = 0; i < nProt; i++) { ref[3*i]= (Math.random()-0.5)*2; ref[3*i+1]=0; ref[3*i+2]=0; }
    for (let a=0;a<nLig;a++){ const idx=nProt+a; ref[3*idx]=2; ref[3*idx+1]=0; ref[3*idx+2]=0; }
    const funnel = new Funnel({ nProt, n, ref, sigma, w0:w, biasFactor:1, rMax, bins });
    // deposit directly at rCenter with height w (biasFactor=1 => w=w0)
    funnel.deposit(rCenter);
    let sumBias = 0;
    for (let k=0;k<bins;k++) sumBias += funnel._bias[k] * delta;
    const relErr2 = Math.abs(sumBias - expected)/expected;
    console.log(`  Funnel deposit sum*Δr = ${sumBias.toFixed(6)}, relErr=${(relErr2*100).toFixed(3)}%`);
    assert(relErr2 < 0.01, `Funnel grid Gaussian within 1% (relErr ${(relErr2*100).toFixed(3)}%)`);
    console.log("\nPASS test_funnel_grid.js");
  }).catch(e => { console.error(e); process.exit(1); });
}

run();
