/**
 * PmfTrajectoryScrubber.tsx — Blueprint component (not wired to index.html yet)
 *
 * Links a 2D free-energy landscape (distance-to-pocket vs contacts/RMSD)
 * to a frame-buffered trajectory viewer with zero main-thread freezing.
 *
 * Design contract vs current app:
 * - Replaces stub scrub in src/main.js:687 + <input> missing in index.html:213-239
 * - Consumes src/recorder.js frames via async TrajStore (chunked, worker-fed)
 * - PMF binning + k-medoid clustering run in a Blob Web Worker (no jank)
 * - Render path is rAF-only + FPS throttle + step-skip + loop + bookmarks
 *
 * Drop-in: `import { PmfTrajectoryScrubber } from './PmfTrajectoryScrubber'`
 * and pass `store` + `onSeek(frameIdx, coords)` which calls
 * `viewer.render(await store.getFrame(i))`.
 */

import React, {
  useCallback,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";

// ---------------------------------------------------------------------------
// Types (mirror src/recorder.js + src/analysis.js concepts)
// ---------------------------------------------------------------------------

export interface TrajStore {
  readonly numFrames: number;
  readonly numAtoms: number;
  /** Must be O(1) chunk fetch; may hit IndexedDB — hence async. */
  getFrame(i: number): Promise<Float32Array>;
  getTimePs(i: number): number;
  /** Precomputed CVs per frame: [distToPocket, nContacts]. */
  getCV(i: number): [number, number];
}

export interface PmfBin {
  freeEnergy: Float32Array; // nX*nY, kT units, min=0
  counts: Uint32Array;
  /** frame indices per bin (capped, e.g. 8/bin for memory) */
  members: number[][];
  nX: number;
  nY: number;
  xRange: [number, number]; // distance Å
  yRange: [number, number]; // contacts
  centroids: { x: number; y: number; frame: number; label: string }[];
}

interface Props {
  store: TrajStore;
  /** Called on every committed seek; parent renders 3D frame. */
  onSeek: (frame: number, coords: Float32Array) => void;
  kT?: number; // kcal/mol, default 0.596 at 300K
}

// ---------------------------------------------------------------------------
// Worker source: PMF histogram + lightweight k-means on (x,y), off-thread.
// ---------------------------------------------------------------------------

const WORKER_SRC = `
onmessage = (e) => {
  const { xs, ys, nX, nY, xRange, yRange, kT, k } = e.data;
  const n = xs.length;
  const counts = new Uint32Array(nX * nY);
  const members = Array.from({ length: nX * nY }, () => []);
  const dx = (xRange[1] - xRange[0]) / nX, dy = (yRange[1] - yRange[0]) / nY;
  for (let i = 0; i < n; i++) {
    let bx = Math.floor((xs[i] - xRange[0]) / dx), by = Math.floor((ys[i] - yRange[0]) / dy);
    bx = Math.max(0, Math.min(nX - 1, bx)); by = Math.max(0, Math.min(nY - 1, by));
    const b = by * nX + bx;
    counts[b]++;
    if (members[b].length < 8) members[b].push(i);
  }
  let cmax = 1; for (let c of counts) if (c > cmax) cmax = c;
  const fe = new Float32Array(nX * nY);
  for (let b = 0; b < fe.length; b++) fe[b] = counts[b] > 0 ? -kT * Math.log(counts[b] / cmax) : Infinity;
  // k-means on (x,y) normalized
  const cx = Array.from({length:k},(_,j)=>xRange[0]+(j+0.5)/k*(xRange[1]-xRange[0]));
  const cy = Array.from({length:k},()=> (yRange[0]+yRange[1])/2);
  const assign = new Int32Array(n);
  for (let it = 0; it < 20; it++) {
    const sx = new Float64Array(k), sy = new Float64Array(k), sn = new Float64Array(k);
    for (let i = 0; i < n; i++) {
      let best = 0, bd = Infinity;
      for (let j = 0; j < k; j++) { const d=(xs[i]-cx[j])**2+(ys[i]-cy[j])**2; if(d<bd){bd=d;best=j;} }
      assign[i]=best; sx[best]+=xs[i]; sy[best]+=ys[i]; sn[best]++;
    }
    for (let j = 0; j < k; j++) if (sn[j]>0){cx[j]=sx[j]/sn[j];cy[j]=sy[j]/sn[j];}
  }
  // medoid per cluster (nearest frame to centroid)
  const labels = ["Bound","Intermediate","Encounter","Unbound"];
  const centroids = cx.map((x,j)=>{
    let bf=0,bd=Infinity;
    for (let i=0;i<n;i++) if(assign[i]===j){const d=(xs[i]-x)**2+(ys[i]-cy[j])**2; if(d<bd){bd=d;bf=i;}}
    return { x, y: cy[j], frame: bf, label: labels[j] ?? ("S"+j) };
  });
  onmessage && postMessage({ freeEnergy: fe.buffer, counts: counts.buffer, members, centroids, nX, nY }, [fe.buffer, counts.buffer]);
};
`;

// Inferno-like LUT without deps
function colorForG(g: number, gmax: number): string {
  if (!isFinite(g)) return "#0f172a";
  const t = Math.max(0, Math.min(1, 1 - g / Math.max(1e-6, gmax)));
  // dark → purple → orange → yellow
  const stops: [number, [number, number, number]][] = [
    [0, [15, 23, 42]], [0.4, [76, 29, 149]], [0.7, [234, 88, 12]], [1, [250, 204, 21]],
  ];
  for (let s = 1; s < stops.length; s++) {
    if (t <= stops[s][0]) {
      const [t0, c0] = stops[s - 1], [t1, c1] = stops[s];
      const u = (t - t0) / (t1 - t0);
      const c = c0.map((v, i) => Math.round(v + (c1[i] - v) * u));
      return `rgb(${c[0]},${c[1]},${c[2]})`;
    }
  }
  return "#facc15";
}

export function PmfTrajectoryScrubber({ store, onSeek, kT = 0.596 }: Props) {
  const [pmf, setPmf] = useState<PmfBin | null>(null);
  const [frame, setFrame] = useState(0);
  const [playing, setPlaying] = useState(false);
  const [fps, setFps] = useState(30);
  const [skip, setSkip] = useState(1);
  const [loop, setLoop] = useState(true);
  const [bookmarks, setBookmarks] = useState<{ frame: number; label: string }[]>([]);
  const heatRef = useRef<HTMLCanvasElement>(null);
  const rafRef = useRef<number>(0);
  const lastDrawRef = useRef(0);
  const frameRef = useRef(0);
  frameRef.current = frame;

  // ---- 1. Off-thread PMF build (runs once per store) -----------------------
  useEffect(() => {
    const n = store.numFrames;
    if (n < 2) return;
    const xs = new Float64Array(n), ys = new Float64Array(n);
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (let i = 0; i < n; i++) {
      const [x, y] = store.getCV(i);
      xs[i] = x; ys[i] = y;
      if (x < x0) x0 = x; if (x > x1) x1 = x;
      if (y < y0) y0 = y; if (y > y1) y1 = y;
    }
    if (x1 - x0 < 1e-6) { x0 -= 1; x1 += 1; }
    if (y1 - y0 < 1e-6) { y0 -= 1; y1 += 1; }
    const blob = new Blob([WORKER_SRC], { type: "text/javascript" });
    const w = new Worker(URL.createObjectURL(blob));
    w.onmessage = (e) => {
      const d = e.data;
      setPmf({
        freeEnergy: new Float32Array(d.freeEnergy),
        counts: new Uint32Array(d.counts),
        members: d.members as number[][],
        nX: d.nX, nY: d.nY,
        xRange: [x0, x1], yRange: [y0, y1],
        centroids: d.centroids,
      });
      w.terminate();
    };
    w.postMessage({ xs, ys, nX: 30, nY: 24, xRange: [x0, x1], yRange: [y0, y1], kT, k: 4 });
    return () => w.terminate();
  }, [store, kT]);

  // ---- 2. Heatmap paint (memoized, main thread only draws rects) ------------
  const gmax = useMemo(() => {
    if (!pmf) return 1;
    let m = 0;
    for (const g of pmf.freeEnergy) if (isFinite(g) && g > m) m = g;
    return m || 1;
  }, [pmf]);

  useEffect(() => {
    const cv = heatRef.current;
    if (!cv || !pmf) return;
    const ctx = cv.getContext("2d");
    if (!ctx) return;
    const W = (cv.width = 300), H = (cv.height = 200);
    const cw = W / pmf.nX, ch = H / pmf.nY;
    for (let by = 0; by < pmf.nY; by++)
      for (let bx = 0; bx < pmf.nX; bx++) {
        ctx.fillStyle = colorForG(pmf.freeEnergy[by * pmf.nX + bx], gmax);
        ctx.fillRect(bx * cw, H - (by + 1) * ch, cw + 0.5, ch + 0.5);
      }
    // centroids
    ctx.font = "10px system-ui";
    for (const c of pmf.centroids) {
      const px = ((c.x - pmf.xRange[0]) / (pmf.xRange[1] - pmf.xRange[0])) * W;
      const py = H - ((c.y - pmf.yRange[0]) / (pmf.yRange[1] - pmf.yRange[0])) * H;
      ctx.strokeStyle = "#fff"; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.arc(px, py, 6, 0, 7); ctx.stroke();
      ctx.fillStyle = "#fff"; ctx.fillText(c.label, px + 8, py + 3);
    }
    // current-frame crosshair
    const [cx, cy] = store.getCV(frameRef.current);
    const px = ((cx - pmf.xRange[0]) / (pmf.xRange[1] - pmf.xRange[0])) * W;
    const py = H - ((cy - pmf.yRange[0]) / (pmf.yRange[1] - pmf.yRange[0])) * H;
    ctx.fillStyle = "#38bdf8";
    ctx.fillRect(px - 1, 0, 2, H); ctx.fillRect(0, py - 1, W, 2);
  }, [pmf, gmax, store, frame]);

  // ---- 3. Seek: PMF click → nearest member frame → 3D snap ------------------
  const seek = useCallback(
    async (i: number) => {
      const n = store.numFrames;
      let j = ((Math.round(i) % n) + n) % n;
      setFrame(j);
      // Chunked async fetch keeps UI free; parent renders (double-buffered).
      const coords = await store.getFrame(j);
      onSeek(j, coords);
    },
    [store, onSeek]
  );

  const onHeatClick = useCallback(
    (e: React.MouseEvent<HTMLCanvasElement>) => {
      if (!pmf) return;
      const r = (e.target as HTMLCanvasElement).getBoundingClientRect();
      const bx = Math.floor(((e.clientX - r.left) / r.width) * pmf.nX);
      const by = Math.floor((1 - (e.clientY - r.top) / r.height) * pmf.nY);
      const b = Math.max(0, Math.min(pmf.nX * pmf.nY - 1, by * pmf.nX + bx));
      const cand = pmf.members[b];
      if (cand?.length) void seek(cand[0]);
    },
    [pmf, seek]
  );

  // ---- 4. Transport: rAF + FPS throttle + step-skip + loop -------------------
  useEffect(() => {
    if (!playing) return;
    let dead = false;
    const tick = async (t: number) => {
      if (dead) return;
      const interval = 1000 / Math.max(1, fps);
      if (t - lastDrawRef.current >= interval) {
        lastDrawRef.current = t;
        let nxt = frameRef.current + Math.max(1, Math.floor(skip));
        if (nxt >= store.numFrames) {
          if (!loop) { setPlaying(false); return; }
          nxt = 0;
        }
        await seek(nxt);
      }
      rafRef.current = requestAnimationFrame(tick);
    };
    rafRef.current = requestAnimationFrame(tick);
    return () => { dead = true; cancelAnimationFrame(rafRef.current); };
  }, [playing, fps, skip, loop, store, seek]);

  // Keyboard: Space play, B bookmark, arrows step
  useEffect(() => {
    const h = (e: KeyboardEvent) => {
      if (/INPUT|SELECT|TEXTAREA/.test((document.activeElement?.tagName ?? ""))) return;
      if (e.code === "Space") { e.preventDefault(); setPlaying((p) => !p); }
      else if (e.key === "b" || e.key === "B")
        setBookmarks((b) => [...b, { frame: frameRef.current, label: `F${frameRef.current}` }].slice(-12));
      else if (e.key === "ArrowRight") void seek(frameRef.current + skip);
      else if (e.key === "ArrowLeft") void seek(frameRef.current - skip);
    };
    window.addEventListener("keydown", h);
    return () => window.removeEventListener("keydown", h);
  }, [seek, skip]);

  return (
    <div style={{ display: "grid", gap: 8, font: "12px system-ui" }}>
      <canvas
        ref={heatRef}
        onClick={onHeatClick}
        title="Click a basin to snap the 3D view to that frame"
        style={{ width: "100%", borderRadius: 8, border: "1px solid #334155", cursor: "crosshair" }}
      />
      <div style={{ display: "flex", gap: 6, flexWrap: "wrap" }}>
        {pmf?.centroids.map((c) => (
          <button key={c.label} onClick={() => void seek(c.frame)} title={`Snap to ${c.label} medoid (frame ${c.frame})`}>
            ○ {c.label}
          </button>
        ))}
      </div>
      <div style={{ display: "flex", gap: 8, alignItems: "center", background: "#0d121d", padding: "6px 10px", borderRadius: 8 }}>
        <button onClick={() => setPlaying((p) => !p)}>{playing ? "⏸ Pause" : "▶ Play"}</button>
        <input
          type="range" min={0} max={Math.max(0, store.numFrames - 1)} value={frame}
          onChange={(e) => void seek(Number(e.target.value))} style={{ flex: 1 }} aria-label="Frame scrubber"
        />
        <span>{frame} / {store.numFrames - 1} · {store.getTimePs(frame).toFixed(1)} ps</span>
      </div>
      <div style={{ display: "flex", gap: 10, alignItems: "center", flexWrap: "wrap" }}>
        <label>FPS <input type="range" min={1} max={60} value={fps} onChange={(e) => setFps(Number(e.target.value))} /> {fps}</label>
        <label>Skip <input type="number" min={1} max={20} value={skip} onChange={(e) => setSkip(Number(e.target.value))} style={{ width: 48 }} /></label>
        <label><input type="checkbox" checked={loop} onChange={(e) => setLoop(e.target.checked)} /> Loop</label>
        <button onClick={() => setBookmarks((b) => [...b, { frame, label: `F${frame}` }].slice(-12))}>+ Bookmark (B)</button>
      </div>
      <div style={{ display: "flex", gap: 4, flexWrap: "wrap" }}>
        {bookmarks.map((b, i) => (
          <button key={i} onClick={() => void seek(b.frame)} title={`Frame ${b.frame}`}>◆ {b.label} ✕</button>
        ))}
      </div>
    </div>
  );
}
