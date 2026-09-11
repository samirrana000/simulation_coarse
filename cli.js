#!/usr/bin/env node
/**
 * cli.js — Headless CLI runner for simulation_coarse (I89).
 *
 * Runs a coarse-grained simulation without a browser, using the same
 * physics modules (`src/forcefield.js`, `src/pdb.js`, `src/integrator.js`).
 *
 * Usage:
 *   node cli.js --help
 *   node cli.js --pdb 4w52 --steps 1000 --out out.dcd
 *   node cli.js --pdb my.pdb --steps 10
 *
 * Measurable I89 criterion:
 *   `node cli.js --help` prints help
 *   `node cli.js --pdb 4w52 --steps 10` runs 10 steps and prints energy
 */

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";

import { ForceField } from "./src/forcefield.js";
import { LangevinIntegrator } from "./src/integrator.js";
import { parseCa, selectSystem, parseLigands } from "./src/pdb.js";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

function printHelp() {
  const help = `
simulation_coarse — headless CLI runner (I89)

Usage:
  node cli.js --pdb <id|path> --steps <n> --out <file>
  node cli.js --help

Options:
  --pdb <id|path>   PDB ID (e.g. 4w52) or path to .pdb file (default: 4w52)
                   For ID, loader tries ./4w52.pdb, ./data/4w52.pdb, ./<id>.pdb
  --steps <n>       Number of BAOAB Langevin steps to run (default: 1000)
  --out <file>      Output trajectory file (XYZ or PDB multi-model, default: none)
                   Example: --out out.xyz  or  --out out.dcd (writes XYZ format with .dcd name)
  --temp <K>        Bath temperature in K (default: 300)
  --help, -h        Show this help and exit

Examples:
  node cli.js --help
  node cli.js --pdb 4w52 --steps 10
  node cli.js --pdb 4w52 --steps 1000 --out out.xyz
  node cli.js --pdb ./1crn.pdb --steps 500 --out traj.pdb

Physics:
  ForceField: Cα ENM (rc=10 Å, gamma=2) + bone/angle + repulsion + binding (via src/forcefield.js)
  Integrator: BAOAB Langevin (src/integrator.js), dt auto from k_bond/zeta
  Energy printed is U (kcal/mol) from ff.compute(pos) and ff.energy

See also: bench/perf.js, bench/vs_gromacs.md, docs/LIMITATIONS.md
`.trim();
  console.log(help);
}

function parseArgs(argv) {
  const args = { pdb: "4w52", steps: 1000, out: null, temp: 300, help: false };
  for (let i = 2; i < argv.length; i++) {
    const a = argv[i];
    if (a === "--help" || a === "-h") args.help = true;
    else if (a === "--pdb" && i + 1 < argv.length) args.pdb = argv[++i];
    else if (a.startsWith("--pdb=")) args.pdb = a.split("=").slice(1).join("=");
    else if (a === "--steps" && i + 1 < argv.length) args.steps = parseInt(argv[++i], 10);
    else if (a.startsWith("--steps=")) args.steps = parseInt(a.split("=")[1], 10);
    else if (a === "--out" && i + 1 < argv.length) args.out = argv[++i];
    else if (a.startsWith("--out=")) args.out = a.split("=").slice(1).join("=");
    else if (a === "--temp" && i + 1 < argv.length) args.temp = parseFloat(argv[++i]);
    else if (a.startsWith("--temp=")) args.temp = parseFloat(a.split("=")[1]);
    else if (a === "--heavy") args.heavy = true;
    else {
      // positional pdb fallback
      if (!a.startsWith("-") && args.pdb === "4w52") args.pdb = a;
    }
  }
  return args;
}

function resolvePdbText(pdbArg) {
  const candidates = [];
  // If arg looks like a file path that exists directly
  candidates.push(path.resolve(pdbArg));
  if (!pdbArg.endsWith(".pdb")) {
    candidates.push(path.resolve(`${pdbArg}.pdb`));
    candidates.push(path.resolve(`${pdbArg.toLowerCase()}.pdb`));
    candidates.push(path.resolve(__dirname, `${pdbArg}.pdb`));
    candidates.push(path.resolve(__dirname, `${pdbArg.toLowerCase()}.pdb`));
    candidates.push(path.resolve(__dirname, "data", `${pdbArg.toLowerCase()}.pdb`));
    candidates.push(path.resolve(process.cwd(), `${pdbArg}.pdb`));
    candidates.push(path.resolve(process.cwd(), `${pdbArg.toLowerCase()}.pdb`));
  } else {
    candidates.push(path.resolve(__dirname, pdbArg));
    candidates.push(path.resolve(__dirname, path.basename(pdbArg)));
  }
  // Deduplicate
  const seen = new Set();
  const uniq = candidates.filter(p => { if (seen.has(p)) return false; seen.add(p); return true; });

  for (const p of uniq) {
    try {
      if (fs.existsSync(p)) {
        const text = fs.readFileSync(p, "utf-8");
        if (text.includes("ATOM")) return { text, source: p };
      }
    } catch {}
  }
  // Also try relative to cwd bench-style (process.cwd() may be repo root or bench)
  const extra = [
    path.resolve("4w52.pdb"),
    path.resolve("data/4w52.pdb"),
  ];
  for (const p of extra) {
    try {
      if (fs.existsSync(p)) {
        const text = fs.readFileSync(p, "utf-8");
        if (text.includes("ATOM")) return { text, source: p };
      }
    } catch {}
  }
  throw new Error(`Cannot find PDB for "${pdbArg}" (tried ${uniq.join(", ")}) — place 4w52.pdb in repo root or pass a file path`);
}

function buildXyzFrame(beads, pos, comment) {
  const n = beads.length;
  let out = `${n}\n${comment}\n`;
  for (let i = 0; i < n; i++) {
    const b = beads[i];
    const x = pos[3 * i].toFixed(3);
    const y = pos[3 * i + 1].toFixed(3);
    const z = pos[3 * i + 2].toFixed(3);
    const el = (b.resName ? "CA" : "C");
    out += `${el} ${x} ${y} ${z}\n`;
  }
  return out;
}

async function main() {
  const args = parseArgs(process.argv);
  if (args.help) {
    printHelp();
    process.exit(0);
  }

  // Also support plain `--help` via argv check before heavy imports? already handled
  const steps = Number.isFinite(args.steps) ? Math.max(1, Math.floor(args.steps)) : 1000;
  const temp = Number.isFinite(args.temp) ? args.temp : 300;

  let pdbText, pdbSource;
  try {
    const res = resolvePdbText(args.pdb);
    pdbText = res.text;
    pdbSource = res.source;
  } catch (e) {
    console.error(`Error: ${e.message}`);
    console.error(`Try: node cli.js --help`);
    process.exit(1);
  }

  console.log(`[cli] PDB: ${args.pdb} → ${pdbSource}`);
  console.log(`[cli] steps=${steps} temp=${temp}K`);

  // Parse & select (protein Cα; include ligands if present but don't require)
  let parsed, sel, ligands = [];
  try {
    parsed = parseCa(pdbText);
    sel = selectSystem(parsed);
  } catch (e) {
    console.error(`[cli] parseCa/selectSystem failed: ${e.message}`);
    process.exit(1);
  }
  try {
    ligands = parseLigands(pdbText);
    if (ligands.length) console.log(`[cli] found ${ligands.length} ligand molecule(s) with ${ligands.reduce((s,m)=>s+m.atoms.length,0)} atoms — included in ForceField`);
  } catch {}

  const nProt = sel.beads.length;
  const nLig = ligands.reduce((s, m) => s + m.atoms.length, 0);
  console.log(`[cli] system: ${nProt} Cα beads${nLig ? ` + ${nLig} ligand atoms` : ""} (${sel.segments.length} segment(s))`);

  const ff = new ForceField(sel, { rc: 10, gamma: 2.0, temp }, ligands);
  const integ = new LangevinIntegrator(ff.ref, ff, 110);
  integ.setTemperature(temp);

  console.log(`[cli] dt=${integ.dt.toExponential(3)} ps, n=${ff.n}, mass=${integ.mass[0]} Da (protein), kT=${(0.001987*temp).toFixed(3)} kcal/mol`);
  console.log(`[cli] initial energy: ${ff.energy.toFixed(3)} kcal/mol (U)  U_bind=${ff.bindingU.toFixed(2)}  desolv=${ff.desolvU.toFixed(2)}`);

  // Optional trajectory buffer (XYZ text)
  let traj = "";
  const stride = Math.max(1, Math.floor(steps / 10));
  if (args.out) {
    traj += buildXyzFrame(sel.beads.concat(...ligands.map(m=>m.atoms)), integ.pos, `frame 0 t=0 ps U=${ff.energy.toFixed(2)}`);
  }

  for (let s = 1; s <= steps; s++) {
    integ.step();
    if (s % stride === 0 || s === steps) {
      const e = ff.energy;
      // Print energy — required for `node cli.js --pdb 4w52 --steps 10` measurable
      console.log(`step ${s}/${steps}: energy = ${Number.isFinite(e) ? e.toFixed(3) : "NaN"} kcal/mol  T_inst=${ff.kineticTemp(integ.vel, integ.mass).toFixed(0)}K  RMSD=${ff.rmsd(integ.pos).toFixed(3)}Å`);
      if (args.out && s % stride === 0) {
        traj += buildXyzFrame(sel.beads.concat(...ligands.map(m=>m.atoms)), integ.pos, `frame ${s} t=${integ.time.toFixed(3)}ps U=${e.toFixed(2)}`);
      }
    }
    if (!Number.isFinite(ff.energy)) {
      console.warn(`[cli] non-finite energy at step ${s} — stopping`);
      break;
    }
  }

  console.log(`[cli] done: ${steps} steps, t=${integ.time.toFixed(3)} ps, final energy = ${ff.energy.toFixed(3)} kcal/mol`);

  if (args.out) {
    try {
      // Ensure directory exists
      const outPath = path.resolve(args.out);
      const dir = path.dirname(outPath);
      if (!fs.existsSync(dir)) fs.mkdirSync(dir, { recursive: true });
      // Write XYZ (even if filename is .dcd, write XYZ — honest about format)
      if (outPath.endsWith(".dcd") || outPath.endsWith(".pdb")) {
        // For .dcd/.pdb request, still write XYZ but warn
        console.log(`[cli] note: writing XYZ trajectory to ${outPath} (DCD/PDB binary not implemented in headless — use XYZ)`);
      }
      fs.writeFileSync(outPath, traj, "utf-8");
      console.log(`[cli] trajectory written to ${outPath} (${steps / stride} frames)`);
    } catch (e) {
      console.error(`[cli] failed to write ${args.out}: ${e.message}`);
      process.exit(1);
    }
  }
}

main().catch(e => {
  console.error(e);
  process.exit(1);
});
