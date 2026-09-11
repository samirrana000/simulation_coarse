/**
 * test_master_stability.js — E48 Master equation stability
 * Propagates 500 steps and checks sum=1±1e-10, no negative.
 * Runnable: node tests/test_master_stability.js
 */

import { ChemicalNetworkModel } from "../src/physics/network.js";

function run(){
  console.log("=== test_master_stability.js — E48 master equation stability ===");
  const model = new ChemicalNetworkModel({ temperature:300, concentrationM:0.001 });
  // start from uniform mixed or pure bulk
  model.probabilities = new Float64Array([0.25,0.25,0.25,0.25]);
  // also test pure bulk
  const cases = [
    { name:"uniform", p0: new Float64Array([0.25,0.25,0.25,0.25]) },
    { name:"bulk", p0: new Float64Array([1,0,0,0]) },
    { name:"bound", p0: new Float64Array([0,0,0,1]) },
  ];
  for(const c of cases){
    console.log(`  Case: ${c.name} p0=[${Array.from(c.p0).join(",")}]`);
    model.probabilities = Float64Array.from(c.p0);
    model.timeNs = 0;
    let maxSumErr=0, minProb=Infinity, maxProb=-Infinity;
    for(let step=0; step<500; step++){
      // propagate small chunk
      model.propagateMasterEquation(1e-7, 500); // 0.1 µs per step, many steps
      let sum=0;
      for(let i=0;i<model.nStates;i++) sum+=model.probabilities[i];
      const err = Math.abs(sum-1);
      if(err>maxSumErr) maxSumErr=err;
      for(let i=0;i<model.nStates;i++){
        if(model.probabilities[i] < minProb) minProb = model.probabilities[i];
        if(model.probabilities[i] > maxProb) maxProb = model.probabilities[i];
        if(model.probabilities[i] < -1e-12){
          console.error(`FAIL: negative probability at step ${step} state ${i} p=${model.probabilities[i]}`);
          process.exit(1);
        }
        if(!Number.isFinite(model.probabilities[i])){
          console.error(`FAIL: non-finite probability at step ${step} state ${i}`);
          process.exit(1);
        }
      }
      if(err > 1e-10){
        console.error(`FAIL: sum deviation ${err.toExponential(3)} >1e-10 at step ${step} sum=${sum}`);
        process.exit(1);
      }
    }
    console.log(`    ✓ 500 steps ok: max|Σ-1|=${maxSumErr.toExponential(3)} minProb=${minProb.toExponential(3)} maxProb=${maxProb.toFixed(4)}`);
  }
  // Also test equilibrium convergence: after long time, p → pi
  model.probabilities = new Float64Array([1,0,0,0]);
  model.propagateMasterEquation(0.01, 500); // 10 ms, should be equilibrium
  let sum=0; for(let i=0;i<model.nStates;i++) sum+=model.probabilities[i];
  console.log(`  Equilibrium check after 10ms: sum=${sum} pi=[${Array.from(model.pi).map(v=>v.toFixed(4)).join(",")}] p=[${Array.from(model.probabilities).map(v=>v.toFixed(4)).join(",")}]`);

  console.log("\nPASS test_master_stability.js — 500 steps sum=1±1e-10, no negative");
}

run();
