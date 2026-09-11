/**
 * test_live_tracking.js — E49 Live tracking link
 * Mocks classification vs true COM/contacts, asserts match.
 * In src/main.js where isLiveTrackingActive classification happens:
 *   const curMacroState = networkModel.classifyPose(parseFloat(cv) ||10, nc, ligRMSD);
 *   if(curMacroState !== networkModel.currentState){ ... update }
 * This test mimics that logic and verifies against true COM/contacts computed from positions.
 * Runnable: node tests/test_live_tracking.js
 */

import fs from "fs";
import path from "path";
import { parseCa, parseLigands, selectSystem } from "../src/pdb.js";
import { ForceField } from "../src/forcefield.js";
import { LangevinIntegrator } from "../src/integrator.js";
import { ChemicalNetworkModel } from "../src/physics/network.js";

function findFile(names){ for(const p of names) if(fs.existsSync(p)) return p; throw new Error(`Cannot find ${names[0]}`); }

function computeFeatures(pos, ff, pocketCOM, ligStart, nLig, nProt){
  // true COM/contacts vs classification inputs
  let lx=0,ly=0,lz=0;
  for(let a=0;a<nLig;a++){const c=3*(ligStart+a); lx+=pos[c]; ly+=pos[c+1]; lz+=pos[c+2];}
  lx/=nLig; ly/=nLig; lz/=nLig;
  const cv = Math.sqrt((lx-pocketCOM[0])**2 + (ly-pocketCOM[1])**2 + (lz-pocketCOM[2])**2);
  let nc=0;
  for(let i=0;i<nProt;i++){
    const xi=pos[3*i],yi=pos[3*i+1],zi=pos[3*i+2];
    for(let a=0;a<nLig;a++){
      const c=3*(ligStart+a);
      const dx=pos[c]-xi, dy=pos[c+1]-yi, dz=pos[c+2]-zi;
      if(dx*dx+dy*dy+dz*dz <30.25) nc++;
    }
  }
  let s=0;
  for(let a=0;a<nLig;a++){const c=3*(ligStart+a); const dx=pos[c]-ff.ref[c], dy=pos[c+1]-ff.ref[c+1], dz=pos[c+2]-ff.ref[c+2]; s+=dx*dx+dy*dy+dz*dz;}
  const ligRMSD = Math.sqrt(s/nLig);
  return { cv, nc, ligRMSD, ligCOM:[lx,ly,lz] };
}

function computePocketCOM(ref, nProt, ligStart, nLig){
  let lx=0,ly=0,lz=0;
  for(let a=0;a<nLig;a++){const c=3*(ligStart+a); lx+=ref[c]; ly+=ref[c+1]; lz+=ref[c+2];}
  lx/=nLig; ly/=nLig; lz/=nLig;
  const pocketIdx=[];
  for(let i=0;i<nProt;i++){const c=3*i; const dx=ref[c]-lx, dy=ref[c+1]-ly, dz=ref[c+2]-lz; if(Math.sqrt(dx*dx+dy*dy+dz*dz)<=8.0) pocketIdx.push(i);}
  if(pocketIdx.length===0){
    const byDist=[];
    for(let i=0;i<nProt;i++){const c=3*i; const dx=ref[c]-lx,dy=ref[c+1]-ly,dz=ref[c+2]-lz; byDist.push([dx*dx+dy*dy+dz*dz,i]);}
    byDist.sort((a,b)=>a[0]-b[0]); for(let k=0;k<Math.min(6,nProt);k++) pocketIdx.push(byDist[k][1]);
  }
  let px=0,py=0,pz=0;
  for(const i of pocketIdx){const c=3*i; px+=ref[c]; py+=ref[c+1]; pz+=ref[c+2];}
  px/=pocketIdx.length; py/=pocketIdx.length; pz/=pocketIdx.length;
  return [px,py,pz];
}

// Mock of main.js live tracking block
function liveTrackingUpdate(networkModel, cv, nc, ligRMSD, isLiveTrackingActive) {
  if (!isLiveTrackingActive) return { updated:false, cur:networkModel.currentState };
  const curMacroState = networkModel.classifyPose(parseFloat(cv) || 10, nc, ligRMSD);
  if (curMacroState !== networkModel.currentState) {
    networkModel.currentState = curMacroState;
    networkModel.probabilities.fill(0);
    networkModel.probabilities[curMacroState]=1.0;
    return { updated:true, cur:curMacroState, prev: networkModel.currentState };
  }
  return { updated:false, cur: curMacroState };
}

function run(){
  console.log("=== test_live_tracking.js — E49 live tracking link ===");
  const pdbPath = findFile(["4w52.pdb","./4w52.pdb","../4w52.pdb", path.join(path.dirname(new URL(import.meta.url).pathname), "../4w52.pdb")]);
  const pdbText = fs.readFileSync(pdbPath,"utf-8");
  const parsed = parseCa(pdbText);
  const ligands = parseLigands(pdbText);
  if(ligands.length===0){ console.log("  No ligands found — using synthetic classification only"); }
  const sel = selectSystem(parsed,{});
  // Use first ligand set if multiple
  const ligSet = ligands.length? ligands.slice(0,1) : [];
  const ff = new ForceField(sel, { rc:10, gamma:2.0, binding:{on:true,holo:true} }, ligSet);
  const nProt = ff.nProt, nLig = ff.nLigAtoms, ligStart = nProt;
  if(nLig===0){
    // synthetic test without traj
    console.log("  No ligand atoms — synthetic mock classification");
    const model = new ChemicalNetworkModel({temperature:300});
    const cases = [
      { cv:4.0, nc:15, rmsd:1.2, expect:3 },
      { cv:7.0, nc:6, rmsd:3.0, expect:2 },
      { cv:12.0, nc:2, rmsd:5.0, expect:1 },
      { cv:25.0, nc:0, rmsd:10.0, expect:0 },
    ];
    for(const t of cases){
      const classified = model.classifyPose(parseFloat(t.cv)||10, t.nc, t.rmsd);
      const mock = liveTrackingUpdate(model, t.cv, t.nc, t.rmsd, true);
      if(classified!==t.expect){ console.error(`FAIL: classify ${t.cv}/${t.nc}/${t.rmsd} => ${classified} expect ${t.expect}`); process.exit(1); }
      if(mock.cur!==t.expect){ console.error(`FAIL: liveTrackingUpdate cur ${mock.cur} != expect ${t.expect}`); process.exit(1); }
      console.log(`  ✓ cv=${t.cv} nc=${t.nc} rmsd=${t.rmsd} => S${classified} liveTracking match`);
      // reset
      model.currentState = 99; // force change detection next
    }
    console.log("\nPASS test_live_tracking.js — mock classification vs true thresholds match (synthetic)");
    return;
  }
  const pocketCOM = computePocketCOM(ff.ref, nProt, ligStart, nLig);
  const integ = new LangevinIntegrator(ff.ref, ff, 110);
  integ.setTemperature(300); integ.setFriction(8);
  const model = new ChemicalNetworkModel({temperature:300});

  let mismatches=0;
  let total=50;
  let isLiveTrackingActive = true;
  for(let frame=0; frame<total; frame++){
    integ.step();
    const { cv, nc, ligRMSD } = computeFeatures(integ.pos, ff, pocketCOM, ligStart, nLig, nProt);
    const trueState = model.classifyPose(cv, nc, ligRMSD);
    // Simulate main.js block: parseFloat(cv)||10 handling
    const mockInputCV = parseFloat(cv.toFixed(2)) || 10;
    const liveState = model.classifyPose(mockInputCV, nc, ligRMSD);
    if(trueState!==liveState){
      mismatches++;
      console.error(`  ✗ frame ${frame}: cv=${cv.toFixed(2)} (mock ${mockInputCV}) nc=${nc} rmsd=${ligRMSD.toFixed(2)} => true S${trueState} vs mock S${liveState}`);
    }
    // also test liveTrackingUpdate logic updates model correctly
    const before = model.currentState;
    const res = liveTrackingUpdate(model, cv, nc, ligRMSD, isLiveTrackingActive);
    if(res.cur!==trueState){ console.error(`FAIL: liveTrackingUpdate cur ${res.cur} != true ${trueState}`); process.exit(1); }
    // verify probabilities updated if changed
    if(res.updated){
      if(model.probabilities[trueState]!==1.0){ console.error(`FAIL: probabilities not set to 1 at ${trueState}`); process.exit(1); }
    }
  }
  console.log(`  Ran ${total} frames of Langevin with live tracking mock`);
  console.log(`  Mismatches true vs mock (parseFloat handling): ${mismatches}/${total}`);
  if(mismatches>0){ console.error(`FAIL: ${mismatches} mismatches in live tracking`); process.exit(1); }
  console.log(`  ✓ All ${total} frames: classification vs true COM/contacts match, liveTrackingActive logic correct`);

  // Test toggle off: when isLiveTrackingActive false, no update
  isLiveTrackingActive=false;
  model.currentState=0; model.probabilities.fill(0); model.probabilities[0]=1;
  const { cv, nc, ligRMSD } = computeFeatures(integ.pos, ff, pocketCOM, ligStart, nLig, nProt);
  const beforeState = model.currentState;
  liveTrackingUpdate(model, 4.0, 15, 1.0, isLiveTrackingActive);
  if(model.currentState!==beforeState){ console.error(`FAIL: live tracking inactive should not update`); process.exit(1); }
  console.log(`  ✓ isLiveTrackingActive=false correctly disables updates`);

  console.log("\nPASS test_live_tracking.js — live tracking classification matches true COM/contacts");
}

run();
