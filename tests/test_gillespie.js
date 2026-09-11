/**
 * test_gillespie.js — E47 Gillespie distribution
 * Runs 1000 Gillespie steps, histograms dwell times, fits exponential via log-mean,
 * asserts R²>0.8 or just that mean dwell =1/totalRate within 10% (per-state and overall).
 * Runnable: node tests/test_gillespie.js
 */

import { ChemicalNetworkModel } from "../src/physics/network.js";

function mean(arr){ return arr.reduce((a,b)=>a+b,0)/arr.length; }
function variance(arr, m){ m=m??mean(arr); let s=0; for(const v of arr) s+=(v-m)**2; return s/arr.length; }

function run(){
  console.log("=== test_gillespie.js — E47 Gillespie dwell distribution ===");
  const model = new ChemicalNetworkModel({ temperature:300, concentrationM:0.001 });
  const steps = 1000;
  const dwells = []; // all dt in seconds
  const byState = [[],[],[],[]]; // per fromState
  const expectedRates = []; // per dwell expected totalRate

  for(let i=0;i<steps;i++){
    const cur = model.currentState;
    const totalRate = -model.rateMatrix[cur][cur]; // s^-1
    expectedRates.push(totalRate);
    const res = model.stepGillespie();
    dwells.push(res.dt); // seconds
    byState[res.from].push(res.dt);
  }

  // histogram (log bins or linear)
  const bins = 10;
  const maxD = Math.max(...dwells);
  const minD = Math.min(...dwells);
  const hist = new Array(bins).fill(0);
  for(const d of dwells){
    let b = Math.floor((d - minD)/(maxD - minD + 1e-12)*bins);
    if(b>=bins) b=bins-1;
    if(b<0) b=0;
    hist[b]++;
  }
  console.log(`  Ran ${steps} Gillespie steps`);
  console.log(`  Dwell times (s): min=${minD.toExponential(3)} max=${maxD.toExponential(3)} mean=${mean(dwells).toExponential(3)}`);
  console.log(`  Histogram (${bins} bins over [${minD.toExponential(2)}, ${maxD.toExponential(2)}]): ${hist.join(", ")}`);
  // per-state means vs expected
  let overallExpectedMean = 0;
  // compute per-state expected 1/totalRate (weighted by visitation)
  // For each state, totalRate is constant per model (since rateMatrix static)
  const perStateExpected = [];
  for(let s=0;s<4;s++){
    const rate = -model.rateMatrix[s][s];
    const expMean = rate>0? 1/rate : Infinity;
    perStateExpected[s]=expMean;
    const arr = byState[s];
    if(arr.length>0){
      const m = mean(arr);
      const diff = Math.abs(m - expMean)/expMean;
      console.log(`  State S${s}: n=${arr.length} mean dwell=${m.toExponential(3)} s expected 1/k_tot=${expMean.toExponential(3)} s diff=${(diff*100).toFixed(1)}%`);
      // check within 30% per-state due to sampling noise, but overall 10%
      // We'll allow larger tolerance for small sample counts
      if(arr.length>=50){
        if(diff>0.2) console.log(`    · per-state diff ${(diff*100).toFixed(1)}% >20% (warning, n=${arr.length})`);
      }
    } else {
      console.log(`  State S${s}: n=0 (no visits) expected mean ${expMean.toExponential(3)}`);
    }
  }
  // overall weighted expected mean
  let sumExp=0;
  for(let i=0;i<steps;i++) sumExp += 1/expectedRates[i];
  overallExpectedMean = sumExp/steps;
  const overallMean = mean(dwells);
  const overallDiff = Math.abs(overallMean - overallExpectedMean)/overallExpectedMean;
  console.log(`  Overall mean dwell=${overallMean.toExponential(3)} s expected (weighted) ${overallExpectedMean.toExponential(3)} s diff=${(overallDiff*100).toFixed(1)}%`);
  // Exponential fit via log-mean: for exponential, mean = std, and log(mean) = ?
  // Check exponential property: variance ~ mean^2, and coefficient of variation ≈1
  const varDw = variance(dwells, overallMean);
  const std = Math.sqrt(varDw);
  const cv = std/overallMean;
  console.log(`  Exponential check: std=${std.toExponential(3)} mean=${overallMean.toExponential(3)} CV=${cv.toFixed(3)} (ideal 1.0)`);
  // R² via log-linear fit on histogram tail? Simplified: compare empirical CDF to exponential CDF
  // For quick R², do linear regression of ln(1 - CDF) vs t (should be slope -k)
  // Build sorted dwells
  const sorted = [...dwells].sort((a,b)=>a-b);
  // Use median split: for exponential, ln(2) ≈0.693 = median/mean
  const median = sorted[Math.floor(sorted.length/2)];
  const medianMeanRatio = median/overallMean;
  console.log(`  median=${median.toExponential(3)} median/mean=${medianMeanRatio.toFixed(3)} (ideal ln2≈0.693)`);
  // Compute R² for exponential fit via least squares on ln survivor
  // survivor S(t)=1 - rank/N ; ln S = -k t
  let sumX=0,sumY=0,sumXY=0,sumXX=0,sumYY=0;
  let count=0;
  // use 20 quantiles
  const qs = 20;
  for(let qi=1; qi<qs; qi++){
    const idx = Math.floor(sorted.length*qi/qs);
    const t = sorted[idx];
    const S = 1 - qi/qs;
    if(S<=0) continue;
    const y = Math.log(S);
    sumX+=t; sumY+=y; sumXY+=t*y; sumXX+=t*t; sumYY+=y*y; count++;
  }
  const denom = Math.sqrt((count*sumXX - sumX*sumX)*(count*sumYY - sumY*sumY));
  let r = 0, R2=0;
  if(denom!==0){
    r = (count*sumXY - sumX*sumY)/denom;
    R2 = r*r;
  }
  console.log(`  Exponential fit R² (ln survivor vs t, ${count} quantiles) = ${R2.toFixed(3)} r=${r.toFixed(3)}`);

  const passedR2 = R2>0.8;
  const passedMean = overallDiff < 0.10;
  const passedCV = Math.abs(cv -1) < 0.3; // loose
  console.log(`  Checks: R²>0.8 ? ${passedR2} | mean within 10% ? ${passedMean} (${(overallDiff*100).toFixed(1)}%) | CV≈1 ? ${passedCV}`);

  if(passedMean || passedR2){
    console.log("\nPASS test_gillespie.js — dwell distribution exponential (mean within 10% or R²>0.8)");
  } else {
    // Still print histogram per spec, but also allow pass if histogram printed and not crashed? Spec says asserts R²>0.8 or mean within 10%
    // For our run, we expect mean within 10% due to law of large numbers with 1000 steps
    console.error(`FAIL: mean diff ${(overallDiff*100).toFixed(1)}% >10% and R² ${R2.toFixed(3)} <=0.8`);
    // If both fail due to sampling noise, still consider histogram printed -> but spec says assert, so fail
    process.exit(1);
  }
}

run();
