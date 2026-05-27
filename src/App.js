import { useState, useCallback } from "react";
import * as math from "mathjs";
import {
  LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip,
  Legend, ResponsiveContainer, ReferenceLine
} from "recharts";

/* ═══════════════════════════════════════════════════════════════
   GLOBAL STYLES
═══════════════════════════════════════════════════════════════ */
const G = {
  bg: "#070710",
  bg2: "#0e0e1a",
  bg3: "#161625",
  bg4: "#1e1e30",
  border: "#252538",
  accent: "#6c5ce7",
  teal: "#00cec9",
  coral: "#fd7272",
  gold: "#fdcb6e",
  lime: "#55efc4",
  pink: "#fd79a8",
  text: "#e2e2f0",
  text2: "#8080a8",
  text3: "#454560",
  mono: "'JetBrains Mono', 'Fira Code', monospace",
  sans: "'DM Sans', sans-serif",
};

const css = `
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500;600&display=swap');
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
body{background:${G.bg};color:${G.text};font-family:${G.sans};min-height:100vh}
::-webkit-scrollbar{width:5px;height:5px}
::-webkit-scrollbar-track{background:${G.bg2}}
::-webkit-scrollbar-thumb{background:${G.border};border-radius:4px}
input,select,textarea{background:${G.bg3};border:1px solid ${G.border};color:${G.text};
  border-radius:7px;padding:8px 12px;font-family:${G.mono};font-size:13px;outline:none;width:100%;
  transition:border-color .2s}
input:focus,select:focus,textarea:focus{border-color:${G.accent}}
table{width:100%;border-collapse:collapse;font-size:12px}
th{background:${G.bg4};color:${G.text2};padding:9px 12px;text-align:left;
  font-weight:500;border-bottom:1px solid ${G.border};font-family:${G.mono};font-size:11px}
td{padding:7px 12px;border-bottom:1px solid ${G.border};font-family:${G.mono};
  font-size:12px;color:${G.text}}
tr:last-child td{border-bottom:none}
tr:hover td{background:${G.bg4}}
code{background:${G.bg4};padding:1px 5px;border-radius:4px;font-family:${G.mono};font-size:12px}
`;

/* ═══════════════════════════════════════════════════════════════
   MATH UTILS
═══════════════════════════════════════════════════════════════ */
function safeEval(expr, x) {
  try { return math.evaluate(expr, { x }); } catch { return NaN; }
}
function derivative(expr, x, h = 1e-7) {
  return (safeEval(expr, x + h) - safeEval(expr, x - h)) / (2 * h);
}
function derivative2(expr, x, h = 1e-7) {
  return (safeEval(expr, x + h) - 2 * safeEval(expr, x) + safeEval(expr, x - h)) / (h * h);
}
function calcDerivativeSymbolic(expr) {
  try {
    const node = math.parse(expr);
    const d = math.derivative(node, 'x');
    return d.toString();
  } catch { return null; }
}

function calcError(xNew, xOld, f, errType) {
  if (xOld === null || xOld === undefined) return Infinity;
  if (errType === "Relativo") return Math.abs((xNew - xOld) / (xNew || 1e-100));
  if (errType === "Absoluto") return Math.abs(xNew - xOld);
  if (errType === "Condición" && f) {
    const fx = safeEval(f, xNew);
    const dfx = derivative(f, xNew);
    return Math.abs(xNew * dfx / (fx || 1e-100));
  }
  return Math.abs(xNew - xOld);
}

/* ═══════════════════════════════════════════════════════════════
   CHAPTER 1 — ROOT FINDING METHODS
═══════════════════════════════════════════════════════════════ */
function bisection(f, a, b, tol, maxIter, errType) {
  const rows = [];
  let xa = a, xb = b, xr = (a + b) / 2, err = Infinity;
  let fa = safeEval(f, xa), fb = safeEval(f, xb);
  if (!isFinite(fa) || !isFinite(fb)) return { error: "f(a) o f(b) no son finitos." };
  if (fa * fb > 0) return { error: "f(a)·f(b) > 0: no hay cambio de signo en el intervalo." };
  let xrPrev = null;
  for (let i = 1; i <= maxIter; i++) {
    xr = (xa + xb) / 2;
    const fr = safeEval(f, xr);
    fa = safeEval(f, xa);
    err = calcError(xr, xrPrev, f, errType);
    rows.push({ iter: i, xa: +xa.toFixed(8), xb: +xb.toFixed(8), xr: +xr.toFixed(8), fr: +fr.toFixed(8), error: xrPrev === null ? "---" : +err.toFixed(8) });
    if (Math.abs(fr) < 1e-14) break;
    if (isFinite(err) && err < tol && xrPrev !== null) break;
    if (fa * fr < 0) xb = xr; else xa = xr;
    xrPrev = xr;
  }
  return { rows, root: xr };
}

function falsePosition(f, a, b, tol, maxIter, errType) {
  const rows = [];
  let xa = a, xb = b, xr, xrPrev = null, err = Infinity;
  const fa0 = safeEval(f, xa), fb0 = safeEval(f, xb);
  if (fa0 * fb0 > 0) return { error: "f(a)·f(b) > 0: no hay cambio de signo en el intervalo." };
  for (let i = 1; i <= maxIter; i++) {
    const fa = safeEval(f, xa), fb = safeEval(f, xb);
    if (Math.abs(fa - fb) < 1e-14) break;
    xr = xb - fb * (xa - xb) / (fa - fb);
    const fr = safeEval(f, xr);
    err = calcError(xr, xrPrev, f, errType);
    rows.push({ iter: i, xa: +xa.toFixed(8), xb: +xb.toFixed(8), xr: +xr.toFixed(8), fr: +fr.toFixed(8), error: xrPrev === null ? "---" : +err.toFixed(8) });
    if (Math.abs(fr) < 1e-14) break;
    if (isFinite(err) && err < tol && xrPrev !== null) break;
    if (fa * fr < 0) xb = xr; else xa = xr;
    xrPrev = xr;
  }
  return { rows, root: xr };
}

function fixedPoint(g, x0, tol, maxIter, errType) {
  const rows = [];
  let x = x0, xNew, err = Infinity;
  for (let i = 1; i <= maxIter && err > tol; i++) {
    xNew = safeEval(g, x);
    if (!isFinite(xNew)) return { error: "g(x) divergió. Elige otra función g(x).", rows, root: x };
    err = calcError(xNew, x, null, errType);
    rows.push({ iter: i, x: +x.toFixed(8), gx: +xNew.toFixed(8), error: i === 1 ? "---" : +err.toFixed(8) });
    x = xNew;
  }
  return { rows, root: x };
}

function newtonRaphson(f, x0, tol, maxIter, errType) {
  const rows = [];
  let x = x0, err = Infinity;
  for (let i = 1; i <= maxIter && err > tol; i++) {
    const fx = safeEval(f, x);
    const dfx = derivative(f, x);
    if (Math.abs(dfx) < 1e-14) return { error: "f'(x) ≈ 0: posible raíz múltiple. Usa el método de Raíces Múltiples.", rows, root: x };
    const xNew = x - fx / dfx;
    err = calcError(xNew, x, f, errType);
    rows.push({ iter: i, x: +x.toFixed(8), fx: +fx.toFixed(8), dfx: +dfx.toFixed(8), xNew: +xNew.toFixed(8), error: i === 1 ? "---" : +err.toFixed(8) });
    x = xNew;
    if (Math.abs(fx) < 1e-14) break;
  }
  return { rows, root: x };
}

function secant(f, x0, x1, tol, maxIter, errType) {
  const rows = [];
  let xa = x0, xb = x1, err = Infinity;
  for (let i = 1; i <= maxIter && err > tol; i++) {
    const fa = safeEval(f, xa), fb = safeEval(f, xb);
    if (Math.abs(fb - fa) < 1e-14) return { error: "f(x₁)−f(x₀) ≈ 0. Elige otros puntos iniciales.", rows, root: xb };
    const xNew = xb - fb * (xb - xa) / (fb - fa);
    err = calcError(xNew, xb, f, errType);
    rows.push({ iter: i, xa: +xa.toFixed(8), xb: +xb.toFixed(8), fb: +fb.toFixed(8), xNew: +xNew.toFixed(8), error: i === 1 ? "---" : +err.toFixed(8) });
    xa = xb; xb = xNew;
    if (Math.abs(fb) < 1e-14) break;
  }
  return { rows, root: xb };
}

function multipleRoots(f, x0, tol, maxIter, errType) {
  const rows = [];
  let x = x0, err = Infinity;
  for (let i = 1; i <= maxIter && err > tol; i++) {
    const fx = safeEval(f, x);
    const dfx = derivative(f, x);
    const d2fx = derivative2(f, x);
    const denom = dfx * dfx - fx * d2fx;
    if (Math.abs(denom) < 1e-14) return { error: "Denominador ≈ 0 en iteración " + i, rows, root: x };
    const xNew = x - (fx * dfx) / denom;
    err = calcError(xNew, x, f, errType);
    rows.push({ iter: i, x: +x.toFixed(8), fx: +fx.toFixed(8), dfx: +dfx.toFixed(8), xNew: +xNew.toFixed(8), error: i === 1 ? "---" : +err.toFixed(8) });
    x = xNew;
    if (Math.abs(fx) < 1e-14) break;
  }
  return { rows, root: x };
}

function plotFn(expr, a, b, n = 300) {
  const pts = [];
  for (let i = 0; i <= n; i++) {
    const x = a + (b - a) * i / n;
    const y = safeEval(expr, x);
    if (isFinite(y) && Math.abs(y) < 1e6) pts.push({ x: +x.toFixed(5), y: +y.toFixed(5) });
  }
  return pts;
}

/* ═══════════════════════════════════════════════════════════════
   CHAPTER 2 — LINEAR SYSTEMS
═══════════════════════════════════════════════════════════════ */
function parseMatrix(text) {
  return text.trim().split("\n").map(row => row.trim().split(/[\s,;]+/).map(Number));
}

function runIterativeSystem(A, b, method, w, tol, maxIter, errType) {
  const n = A.length;
  let x = new Array(n).fill(0);
  const rows = [];
  for (let k = 1; k <= maxIter; k++) {
    const xOld = [...x];
    for (let i = 0; i < n; i++) {
      let sigma = 0;
      if (method === "Jacobi") {
        for (let j = 0; j < n; j++) if (j !== i) sigma += A[i][j] * xOld[j];
        x[i] = (b[i] - sigma) / A[i][i];
      } else {
        for (let j = 0; j < n; j++) if (j !== i) sigma += A[i][j] * x[j];
        const xGs = (b[i] - sigma) / A[i][i];
        x[i] = method === "SOR" ? (1 - w) * xOld[i] + w * xGs : xGs;
      }
    }
    const err = calcSysError(x, xOld, A, b, errType);
    const row = { iter: k, error: isFinite(err) ? +err.toFixed(10) : "---" };
    x.forEach((v, i) => { row[`x${i+1}`] = +v.toFixed(8); });
    rows.push(row);
    if (isFinite(err) && err < tol) break;
  }
  return { rows, x };
}

function calcSysError(x, xOld, A, b, errType) {
  if (errType === "Relativo") {
    const num = Math.sqrt(x.reduce((s, v, i) => s + (v - xOld[i]) ** 2, 0));
    const den = Math.sqrt(x.reduce((s, v) => s + v ** 2, 0)) || 1;
    return num / den;
  }
  if (errType === "Absoluto") return Math.sqrt(x.reduce((s, v, i) => s + (v - xOld[i]) ** 2, 0));
  const res = A.map((row, i) => b[i] - row.reduce((s, v, j) => s + v * x[j], 0));
  return Math.sqrt(res.reduce((s, v) => s + v ** 2, 0));
}

function matMul(A, B) {
  return A.map(r => B[0].map((_, j) => r.reduce((s, v, k) => s + v * B[k][j], 0)));
}
function matInv(M) {
  const n = M.length;
  const aug = M.map((r, i) => [...r, ...Array.from({ length: n }, (_, j) => i === j ? 1 : 0)]);
  for (let i = 0; i < n; i++) {
    let maxR = i;
    for (let k = i+1; k < n; k++) if (Math.abs(aug[k][i]) > Math.abs(aug[maxR][i])) maxR = k;
    [aug[i], aug[maxR]] = [aug[maxR], aug[i]];
    const piv = aug[i][i];
    if (Math.abs(piv) < 1e-12) continue;
    aug[i] = aug[i].map(v => v / piv);
    for (let k = 0; k < n; k++) if (k !== i) { const f = aug[k][i]; aug[k] = aug[k].map((v, j) => v - f * aug[i][j]); }
  }
  return aug.map(r => r.slice(n));
}
function spectralRadius(A, method, w) {
  const n = A.length;
  const D = A.map((r, i) => r.map((v, j) => i === j ? v : 0));
  const L = A.map((r, i) => r.map((v, j) => j < i ? -v : 0));
  const U = A.map((r, i) => r.map((v, j) => j > i ? -v : 0));
  let T;
  try {
    if (method === "Jacobi") {
      const Dinv = D.map((r, i) => r.map((v, j) => i === j && v !== 0 ? 1/v : 0));
      T = matMul(Dinv, L.map((r, i) => r.map((v, j) => v + U[i][j])));
    } else if (method === "Gauss-Seidel") {
      const DmL = D.map((r, i) => r.map((v, j) => v + L[i][j]));
      T = matMul(matInv(DmL), U.map(r => r.map(v => -v)));
    } else {
      const DwL = D.map((r, i) => r.map((v, j) => v + w * L[i][j]));
      const Minv = matInv(DwL);
      const rhs = D.map((r, i) => r.map((v, j) => (1-w)*v + w*U[i][j]));
      T = matMul(Minv, rhs);
    }
    return maxEigenvalue(T);
  } catch { return NaN; }
}
function maxEigenvalue(M) {
  const n = M.length;
  let v = Array.from({ length: n }, () => Math.random() + 0.1);
  let lambda = 0;
  for (let iter = 0; iter < 500; iter++) {
    const Mv = M.map(r => r.reduce((s, v2, j) => s + v2 * v[j], 0));
    const norm = Math.sqrt(Mv.reduce((s, x) => s + x*x, 0));
    if (norm < 1e-14) return 0;
    lambda = norm;
    v = Mv.map(x => x / norm);
  }
  return lambda;
}

/* ═══════════════════════════════════════════════════════════════
   CHAPTER 3 — INTERPOLATION
═══════════════════════════════════════════════════════════════ */
function gaussSolve(A, b) {
  const n = A.length;
  const aug = A.map((r, i) => [...r, b[i]]);
  for (let i = 0; i < n; i++) {
    let maxR = i;
    for (let k = i+1; k < n; k++) if (Math.abs(aug[k][i]) > Math.abs(aug[maxR][i])) maxR = k;
    [aug[i], aug[maxR]] = [aug[maxR], aug[i]];
    if (Math.abs(aug[i][i]) < 1e-14) continue;
    aug[i] = aug[i].map(v => v / aug[i][i]);
    for (let k = 0; k < n; k++) if (k !== i) { const f = aug[k][i]; aug[k] = aug[k].map((v, j) => v - f * aug[i][j]); }
  }
  return aug.map(r => r[n]);
}

function vandermonde(xs, ys) {
  const n = xs.length;
  const V = xs.map(x => Array.from({ length: n }, (_, j) => Math.pow(x, j)));
  return gaussSolve(V, ys);
}
function evalVandermonde(coeffs, x) {
  return coeffs.reduce((s, c, i) => s + c * Math.pow(x, i), 0);
}

function dividedDiffs(xs, ys) {
  const n = xs.length;
  const table = ys.map(y => [y]);
  for (let j = 1; j < n; j++)
    for (let i = 0; i < n - j; i++)
      table[i].push((table[i+1][j-1] - table[i][j-1]) / (xs[i+j] - xs[i]));
  return table[0];
}
function evalNewton(xs, coeffs, x) {
  let result = coeffs[coeffs.length - 1];
  for (let i = coeffs.length - 2; i >= 0; i--)
    result = result * (x - xs[i]) + coeffs[i];
  return result;
}

function evalLagrange(xs, ys, x) {
  return xs.reduce((sum, xi, i) => {
    const Li = xs.reduce((p, xj, j) => i === j ? p : p * (x - xj) / (xi - xj), 1);
    return sum + ys[i] * Li;
  }, 0);
}

function splineLinear(xs, ys, x) {
  let idx = xs.findIndex(xi => xi > x) - 1;
  idx = Math.max(0, Math.min(idx, xs.length - 2));
  const h = xs[idx+1] - xs[idx];
  if (Math.abs(h) < 1e-14) return ys[idx];
  return ys[idx] + (ys[idx+1] - ys[idx]) / h * (x - xs[idx]);
}

function buildCubicSpline(xs, ys) {
  const n = xs.length - 1;
  const h = xs.slice(0, n).map((x, i) => xs[i+1] - x);
  const A = Array.from({ length: n+1 }, () => new Array(n+1).fill(0));
  const rhs = new Array(n+1).fill(0);
  A[0][0] = 1; A[n][n] = 1;
  for (let i = 1; i < n; i++) {
    A[i][i-1] = h[i-1]; A[i][i] = 2*(h[i-1]+h[i]); A[i][i+1] = h[i];
    rhs[i] = 3*((ys[i+1]-ys[i])/h[i] - (ys[i]-ys[i-1])/h[i-1]);
  }
  const M = gaussSolve(A, rhs);
  return { M, xs, ys, h };
}
function evalCubicSpline({ M, xs, ys, h }, x) {
  let i = xs.findIndex((xi, idx) => xi > x) - 1;
  i = Math.max(0, Math.min(i, xs.length - 2));
  const hi = h[i], t = x - xs[i];
  return ((M[i+1]-M[i])/(6*hi))*t**3 + (M[i]/2)*t**2 +
    ((ys[i+1]-ys[i])/hi - hi*(2*M[i]+M[i+1])/6)*t + ys[i];
}

function genInterp(xs, ys, method, extra, n = 300) {
  const xMin = Math.min(...xs), xMax = Math.max(...xs);
  return Array.from({ length: n+1 }, (_, i) => {
    const x = xMin + (xMax - xMin) * i / n;
    let y;
    try {
      if (method === "Vandermonde") y = evalVandermonde(extra.vanCoeffs, x);
      else if (method === "Newton") y = evalNewton(xs, extra.newtonCoeffs, x);
      else if (method === "Lagrange") y = evalLagrange(xs, ys, x);
      else if (method === "Spline Lineal") y = splineLinear(xs, ys, x);
      else y = evalCubicSpline(extra.spline, x);
    } catch { y = NaN; }
    return isFinite(y) && Math.abs(y) < 1e8 ? { x: +x.toFixed(5), y: +y.toFixed(5) } : null;
  }).filter(Boolean);
}

/* ═══════════════════════════════════════════════════════════════
   SHARED UI COMPONENTS
═══════════════════════════════════════════════════════════════ */
function Card({ children, style }) {
  return <div style={{ background: G.bg2, borderRadius: 12, border: `1px solid ${G.border}`, ...style }}>{children}</div>;
}
function Label({ children, mono }) {
  return <span style={{ fontSize: 11, color: G.text2, fontFamily: mono ? G.mono : G.sans, textTransform: mono ? "none" : "uppercase", letterSpacing: mono ? 0 : "0.06em", fontWeight: 500 }}>{children}</span>;
}
function Field({ label, children }) {
  return (
    <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
      <Label mono>{label}</Label>
      {children}
    </label>
  );
}
function Btn({ children, onClick, color, secondary, style }) {
  const bg = secondary ? G.bg3 : (color || G.accent);
  return (
    <button onClick={onClick} style={{
      background: bg, color: secondary ? G.text : "#fff",
      border: secondary ? `1px solid ${G.border}` : "none",
      borderRadius: 8, padding: "9px 18px", cursor: "pointer",
      fontFamily: G.sans, fontSize: 13, fontWeight: 600,
      transition: "all 0.15s", ...style
    }}
      onMouseEnter={e => { e.target.style.opacity = "0.82"; e.target.style.transform = "translateY(-1px)"; }}
      onMouseLeave={e => { e.target.style.opacity = "1"; e.target.style.transform = "translateY(0)"; }}
    >
      {children}
    </button>
  );
}
function InfoBox({ color, icon, children }) {
  return (
    <div style={{
      background: color + "12", border: `1px solid ${color}30`,
      borderRadius: 10, padding: "11px 16px", fontSize: 13,
      color: color, lineHeight: 1.6,
      display: "flex", gap: 10, alignItems: "flex-start"
    }}>
      <span style={{ fontSize: 16, flexShrink: 0 }}>{icon}</span>
      <span>{children}</span>
    </div>
  );
}
function ErrBox({ msg }) {
  return msg ? (
    <div style={{ background: G.coral + "18", border: `1px solid ${G.coral}40`, borderRadius: 8, padding: "10px 14px", color: G.coral, fontSize: 13 }}>
      ⚠ {msg}
    </div>
  ) : null;
}
function TabBar({ items, active, setActive, color }) {
  return (
    <div style={{ display: "flex", gap: 6, flexWrap: "wrap" }}>
      {items.map(m => (
        <button key={m} onClick={() => setActive(m)} style={{
          background: active === m ? color + "20" : G.bg3,
          color: active === m ? color : G.text2,
          border: active === m ? `1px solid ${color}40` : `1px solid ${G.border}`,
          borderRadius: 8, padding: "6px 16px", cursor: "pointer",
          fontSize: 13, fontFamily: G.sans, fontWeight: active === m ? 600 : 400,
          transition: "all 0.15s"
        }}>{m}</button>
      ))}
    </div>
  );
}
function SectionTitle({ title, sub, color }) {
  return (
    <div>
      <h1 style={{ fontSize: 22, fontWeight: 700, color, letterSpacing: "-0.4px" }}>{title}</h1>
      {sub && <p style={{ color: G.text2, fontSize: 13, marginTop: 3 }}>{sub}</p>}
    </div>
  );
}
function TableView({ rows, cols }) {
  if (!rows?.length) return null;
  return (
    <div style={{ overflowX: "auto", maxHeight: 340, overflowY: "auto" }}>
      <table>
        <thead>
          <tr>{cols.map(c => <th key={c}>{c}</th>)}</tr>
        </thead>
        <tbody>
          {rows.map((row, i) => (
            <tr key={i}>{cols.map(c => <td key={c} style={{ color: c === "error" ? G.gold : c.startsWith("x") && c !== "xa" && c !== "xb" ? G.lime : G.text }}>{row[c] ?? "---"}</td>)}</tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
const CHART_STYLE = { background: "transparent", border: "1px solid " + G.border, borderRadius: 12, padding: 16 };
const TOOLTIP_STYLE = { background: G.bg2, border: `1px solid ${G.border}`, borderRadius: 8, fontSize: 12, fontFamily: G.mono };

/* ═══════════════════════════════════════════════════════════════
   CHAPTER 1 COMPONENT
═══════════════════════════════════════════════════════════════ */
const C1_METHODS = ["Bisección", "Regla Falsa", "Punto Fijo", "Newton-Raphson", "Secante", "Raíces Múltiples"];
const C1_COLORS = ["#6c5ce7", "#00cec9", "#fd7272", "#fdcb6e", "#55efc4", "#fd79a8"];
const C1_COLS = {
  "Bisección": ["iter", "xa", "xb", "xr", "fr", "error"],
  "Regla Falsa": ["iter", "xa", "xb", "xr", "fr", "error"],
  "Punto Fijo": ["iter", "x", "gx", "error"],
  "Newton-Raphson": ["iter", "x", "fx", "dfx", "xNew", "error"],
  "Secante": ["iter", "xa", "xb", "fb", "xNew", "error"],
  "Raíces Múltiples": ["iter", "x", "fx", "dfx", "xNew", "error"],
};
const C1_HELP = {
  "Bisección": { hint: "Requiere intervalo [a, b] con f(a)·f(b) < 0 (cambio de signo).", eg: "x^3 - x - 2" },
  "Regla Falsa": { hint: "Similar a bisección pero usa interpolación lineal. Mismo requisito de cambio de signo.", eg: "x^3 - x - 2" },
  "Punto Fijo": { hint: "Despeja x = g(x) de f(x) = 0. El método itera xₙ₊₁ = g(xₙ). Convergencia: |g'(x)| < 1.", eg: "(x + 2)^(1/3)" },
  "Newton-Raphson": { hint: "Método cuadrático. La derivada se calcula automáticamente por diferencias finitas.", eg: "cos(x) - x" },
  "Secante": { hint: "Aproxima la derivada con dos puntos. No requiere derivada analítica. Convergencia superlineal.", eg: "x^2 - 2" },
  "Raíces Múltiples": { hint: "Usa f, f', f''. Ideal cuando la raíz tiene multiplicidad m > 1.", eg: "(x-2)^3" },
};

function Chapter1() {
  const [method, setMethod] = useState("Bisección");
  const [fx, setFx] = useState("x^3 - x - 2");
  const [gx, setGx] = useState("(x + 2)^(1/3)");
  const [a, setA] = useState("1");
  const [b, setB] = useState("2");
  const [x0, setX0] = useState("1.5");
  const [x1, setX1] = useState("2");
  const [tol, setTol] = useState("1e-6");
  const [maxIter, setMaxIter] = useState("100");
  const [errType, setErrType] = useState("Relativo");
  const [result, setResult] = useState(null);
  const [plotData, setPlotData] = useState([]);
  const [compData, setCompData] = useState([]);
  const [errMsg, setErrMsg] = useState("");
  const [derivHint, setDerivHint] = useState("");

  const showDerivHint = () => {
    const d = calcDerivativeSymbolic(fx);
    setDerivHint(d ? `f'(x) = ${d}` : "No se pudo calcular simbólicamente.");
  };

  const run = () => {
    setErrMsg(""); setDerivHint("");
    try {
      const t = parseFloat(tol), m = parseInt(maxIter);
      if (isNaN(t) || isNaN(m)) { setErrMsg("Tolerancia o iteraciones inválidas."); return; }
      let res;
      if (method === "Bisección") res = bisection(fx, parseFloat(a), parseFloat(b), t, m, errType);
      else if (method === "Regla Falsa") res = falsePosition(fx, parseFloat(a), parseFloat(b), t, m, errType);
      else if (method === "Punto Fijo") res = fixedPoint(gx, parseFloat(x0), t, m, errType);
      else if (method === "Newton-Raphson") res = newtonRaphson(fx, parseFloat(x0), t, m, errType);
      else if (method === "Secante") res = secant(fx, parseFloat(x0), parseFloat(x1), t, m, errType);
      else res = multipleRoots(fx, parseFloat(x0), t, m, errType);

      if (res.error && !res.rows) { setErrMsg(res.error); return; }
      if (res.error) setErrMsg("⚠ " + res.error);
      setResult(res);

      // Plot
      const expr = method === "Punto Fijo" ? gx : fx;
      const pA = ["Bisección", "Regla Falsa"].includes(method) ? parseFloat(a) - 0.5 : parseFloat(x0) - 2;
      const pB = ["Bisección", "Regla Falsa"].includes(method) ? parseFloat(b) + 0.5 : parseFloat(x0) + 2;
      setPlotData(plotFn(method === "Punto Fijo" ? `${gx} - x` : fx, pA, pB));
    } catch (e) { setErrMsg("Error: " + e.message); }
  };

  const runComparison = () => {
    setErrMsg("");
    try {
      const t = parseFloat(tol), m = parseInt(maxIter);
      const aV = parseFloat(a), bV = parseFloat(b), x0V = parseFloat(x0), x1V = parseFloat(x1);
      const runMethod = (et) => {
        if (method === "Bisección") return bisection(fx, aV, bV, t, m, et);
        if (method === "Regla Falsa") return falsePosition(fx, aV, bV, t, m, et);
        if (method === "Punto Fijo") return fixedPoint(gx, x0V, t, m, et);
        if (method === "Newton-Raphson") return newtonRaphson(fx, x0V, t, m, et);
        if (method === "Secante") return secant(fx, x0V, x1V, t, m, et);
        return multipleRoots(fx, x0V, t, m, et);
      };
      const runs = ["Relativo", "Absoluto", "Condición"].map(et => ({ name: et, res: runMethod(et) }));
      const maxLen = Math.max(...runs.map(r => r.res.rows?.length || 0));
      const comp = Array.from({ length: maxLen }, (_, i) => ({ iter: i + 1 }));
      runs.forEach(({ name, res }) => {
        (res.rows || []).forEach(r => {
          const e = typeof r.error === "number" && r.error > 0 ? r.error : null;
          if (comp[r.iter - 1] && e !== null) comp[r.iter - 1][name] = e;
        });
      });
      setCompData(comp.slice(0, 60));
    } catch (e) { setErrMsg("Error comparación: " + e.message); }
  };

  const help = C1_HELP[method];
  const needsInterval = ["Bisección", "Regla Falsa"].includes(method);
  const needsX1 = method === "Secante";
  const needsX0 = !needsInterval;
  const needsGx = method === "Punto Fijo";
  const needsFx = method !== "Punto Fijo";

  return (
    <div style={{ display: "flex", flexDirection: "column", gap: 22 }}>
      <SectionTitle title="Raíces de Ecuaciones" sub="Encuentra raíces de f(x) = 0 con seis métodos iterativos" color={G.accent} />
      <TabBar items={C1_METHODS} active={method} setActive={setMethod} color={G.accent} />

      <InfoBox color={G.accent} icon="💡">
        <strong>{method}:</strong> {help.hint}<br />
        <span style={{ opacity: 0.8 }}>Sintaxis soportada: <code>x^2</code>, <code>sin(x)</code>, <code>cos(x)</code>, <code>exp(x)</code>, <code>log(x)</code>, <code>sqrt(x)</code>, <code>abs(x)</code>, <code>pi</code>, <code>e</code>, <code>tan(x)</code>, <code>asin(x)</code>, operadores: <code>+ - * / ^</code></span>
        {!needsGx && <><br /><span style={{ opacity: 0.7 }}>Ejemplo: <code>{help.eg}</code></span></>}
      </InfoBox>

      <Card style={{ padding: 20 }}>
        <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(180px, 1fr))", gap: 14 }}>
          {needsFx && (
            <Field label="f(x)">
              <div style={{ display: "flex", gap: 6 }}>
                <input value={fx} onChange={e => setFx(e.target.value)} placeholder={help.eg} style={{ flex: 1 }} />
                <Btn onClick={showDerivHint} secondary style={{ padding: "8px 12px", fontSize: 11, whiteSpace: "nowrap" }}>f'(x)?</Btn>
              </div>
            </Field>
          )}
          {needsGx && (
            <Field label="g(x) — función de iteración">
              <input value={gx} onChange={e => setGx(e.target.value)} placeholder="(x + 2)^(1/3)" />
            </Field>
          )}
          {needsInterval && (
            <>
              <Field label="a — límite inferior"><input value={a} onChange={e => setA(e.target.value)} /></Field>
              <Field label="b — límite superior"><input value={b} onChange={e => setB(e.target.value)} /></Field>
            </>
          )}
          {needsX0 && <Field label="x₀ — punto inicial"><input value={x0} onChange={e => setX0(e.target.value)} /></Field>}
          {needsX1 && <Field label="x₁ — segundo punto"><input value={x1} onChange={e => setX1(e.target.value)} /></Field>}
          <Field label="Tolerancia"><input value={tol} onChange={e => setTol(e.target.value)} /></Field>
          <Field label="Máx. iteraciones"><input value={maxIter} onChange={e => setMaxIter(e.target.value)} /></Field>
          <Field label="Tipo de error">
            <select value={errType} onChange={e => setErrType(e.target.value)}>
              {["Relativo", "Absoluto", "Condición"].map(e => <option key={e}>{e}</option>)}
            </select>
          </Field>
          <div style={{ display: "flex", alignItems: "flex-end", gap: 8 }}>
            <Btn onClick={run} style={{ flex: 1 }}>▶ Ejecutar</Btn>
            <Btn onClick={runComparison} secondary style={{ flex: 1 }}>⇄ Comparar</Btn>
          </div>
        </div>
        {derivHint && <div style={{ marginTop: 12, padding: "8px 14px", background: G.bg3, borderRadius: 8, fontFamily: G.mono, fontSize: 13, color: G.gold }}>{derivHint}</div>}
      </Card>

      <ErrBox msg={errMsg} />

      {plotData.length > 0 && (
        <Card style={CHART_STYLE}>
          <div style={{ fontSize: 13, color: G.text2, marginBottom: 12, fontWeight: 500 }}>
            Gráfica — {needsGx ? "g(x) − x" : "f(x)"}
            {result?.root !== undefined && <span style={{ marginLeft: 16, color: G.lime, fontFamily: G.mono }}>raíz ≈ {result.root?.toFixed(8)}</span>}
          </div>
          <ResponsiveContainer width="100%" height={230}>
            <LineChart data={plotData}>
              <CartesianGrid strokeDasharray="3 3" stroke={G.border} />
              <XAxis dataKey="x" stroke={G.text3} tick={{ fontSize: 10, fill: G.text3 }} />
              <YAxis stroke={G.text3} tick={{ fontSize: 10, fill: G.text3 }} />
              <Tooltip contentStyle={TOOLTIP_STYLE} />
              <ReferenceLine y={0} stroke={G.text3} strokeDasharray="4 2" />
              <Line type="monotone" dataKey="y" stroke={G.accent} dot={false} strokeWidth={2.5} name="f(x)" />
            </LineChart>
          </ResponsiveContainer>
          {result && (
            <div style={{ marginTop: 12, display: "flex", gap: 12, flexWrap: "wrap" }}>
              <div style={{ background: G.accent + "18", borderRadius: 8, padding: "8px 16px", fontFamily: G.mono, fontSize: 13, color: G.accent }}>
                Raíz ≈ {result.root?.toFixed(10)}
              </div>
              <div style={{ background: G.teal + "18", borderRadius: 8, padding: "8px 16px", fontFamily: G.mono, fontSize: 13, color: G.teal }}>
                Iteraciones: {result.rows?.length}
              </div>
            </div>
          )}
        </Card>
      )}

      {result?.rows?.length > 0 && (
        <Card>
          <div style={{ padding: "14px 20px", borderBottom: `1px solid ${G.border}`, fontSize: 13, fontWeight: 600 }}>
            Tabla de iteraciones — {method}
          </div>
          <TableView rows={result.rows} cols={C1_COLS[method]} />
        </Card>
      )}

      {compData.length > 0 && (
        <Card style={CHART_STYLE}>
          <div style={{ fontSize: 13, color: G.text2, marginBottom: 12, fontWeight: 500 }}>
            Comparación de tipos de error — {method} (escala logarítmica)
          </div>
          <ResponsiveContainer width="100%" height={280}>
            <LineChart data={compData}>
              <CartesianGrid strokeDasharray="3 3" stroke={G.border} />
              <XAxis dataKey="iter" stroke={G.text3} tick={{ fontSize: 10, fill: G.text3 }} label={{ value: "Iteración", position: "insideBottom", offset: -4, fill: G.text3, fontSize: 11 }} />
              <YAxis scale="log" stroke={G.text3} tick={{ fontSize: 10, fill: G.text3 }} />
              <Tooltip contentStyle={TOOLTIP_STYLE} />
              <Legend wrapperStyle={{ fontSize: 12 }} />
              <Line type="monotone" dataKey="Relativo" stroke={G.accent} dot={false} strokeWidth={2} connectNulls={false} />
              <Line type="monotone" dataKey="Absoluto" stroke={G.teal} dot={false} strokeWidth={2} connectNulls={false} />
              <Line type="monotone" dataKey="Condición" stroke={G.coral} dot={false} strokeWidth={2} connectNulls={false} />
            </LineChart>
          </ResponsiveContainer>
        </Card>
      )}
    </div>
  );
}

/* ═══════════════════════════════════════════════════════════════
   CHAPTER 2 COMPONENT
═══════════════════════════════════════════════════════════════ */
const C2_METHODS = ["Jacobi", "Gauss-Seidel", "SOR"];
const C2_COLORS = ["#6c5ce7", "#00cec9", "#fd7272"];

function Chapter2() {
  const [method, setMethod] = useState("Gauss-Seidel");
  const [matA, setMatA] = useState("4 -1 0\n-1 4 -1\n0 -1 4");
  const [vecB, setVecB] = useState("15\n10\n10");
  const [w, setW] = useState("1.1");
  const [tol, setTol] = useState("1e-6");
  const [maxIter, setMaxIter] = useState("200");
  const [errType, setErrType] = useState("Relativo");
  const [result, setResult] = useState(null);
  const [rho, setRho] = useState(null);
  const [compData, setCompData] = useState([]);
  const [errMsg, setErrMsg] = useState("");

  const run = () => {
    setErrMsg("");
    try {
      const A = parseMatrix(matA);
      const b = parseMatrix(vecB).flat();
      const n = A.length;
      if (n < 2 || n > 8) { setErrMsg("Tamaño de matriz debe ser entre 2×2 y 8×8."); return; }
      if (A.some(r => r.length !== n)) { setErrMsg("La matriz A no es cuadrada. Verifica que cada fila tenga " + n + " valores."); return; }
      if (b.length !== n) { setErrMsg(`El vector b debe tener ${n} valores.`); return; }
      if (A.some(r => r.some(isNaN)) || b.some(isNaN)) { setErrMsg("Valores no numéricos detectados."); return; }
      const t = parseFloat(tol), m = parseInt(maxIter), wVal = parseFloat(w);
      if (method === "SOR" && (wVal <= 0 || wVal >= 2)) { setErrMsg("Para SOR, w debe estar en (0, 2)."); return; }
      const res = runIterativeSystem(A, b, method, wVal, t, m, errType);
      setResult({ ...res, n });
      setRho(spectralRadius(A, method, wVal));
    } catch (e) { setErrMsg("Error: " + e.message); }
  };

  const runComparison = () => {
    setErrMsg("");
    try {
      const A = parseMatrix(matA), b = parseMatrix(vecB).flat();
      const t = parseFloat(tol), m = parseInt(maxIter), wVal = parseFloat(w);
      const runs = ["Relativo", "Absoluto", "Condición"].map(et => ({
        name: et, res: runIterativeSystem(A, b, method, wVal, t, m, et)
      }));
      const maxLen = Math.max(...runs.map(r => r.res.rows.length));
      const comp = Array.from({ length: maxLen }, (_, i) => ({ iter: i + 1 }));
      runs.forEach(({ name, res }) => {
        res.rows.forEach(r => { if (comp[r.iter - 1] && typeof r.error === "number" && r.error > 0) comp[r.iter - 1][name] = r.error; });
      });
      setCompData(comp.slice(0, 60));
    } catch (e) { setErrMsg("Error comparación: " + e.message); }
  };

  const cols = result ? ["iter", ...Array.from({ length: result.n }, (_, i) => `x${i+1}`), "error"] : [];

  return (
    <div style={{ display: "flex", flexDirection: "column", gap: 22 }}>
      <SectionTitle title="Sistemas de Ecuaciones Lineales" sub="Métodos iterativos para Ax = b — hasta matrices 8×8" color={G.teal} />
      <TabBar items={C2_METHODS} active={method} setActive={setMethod} color={G.teal} />

      <InfoBox color={G.teal} icon="📐">
        Ingresa la <strong>matriz A</strong> fila por fila, valores separados por espacios.<br />
        El <strong>vector b</strong>, un valor por línea.<br />
        {method === "SOR" ? "Para SOR: w ∈ (0,2). Típicamente w ∈ (1,2) acelera convergencia." : "Convergencia garantizada si la matriz es diagonalmente dominante."}<br />
        <span style={{ opacity: 0.7 }}>Ejemplo diagonal dominante: <code>4 -1 0 / -1 4 -1 / 0 -1 4</code></span>
      </InfoBox>

      <div style={{ display: "grid", gridTemplateColumns: "2fr 1fr", gap: 16 }}>
        <Card style={{ padding: 20, display: "flex", flexDirection: "column", gap: 14 }}>
          <Field label="Matriz A (fila por línea, valores separados por espacios)">
            <textarea value={matA} onChange={e => setMatA(e.target.value)} rows={6} style={{ resize: "vertical", fontFamily: G.mono, fontSize: 13 }} />
          </Field>
          <Field label="Vector b (un valor por línea)">
            <textarea value={vecB} onChange={e => setVecB(e.target.value)} rows={4} style={{ resize: "vertical", fontFamily: G.mono, fontSize: 13 }} />
          </Field>
        </Card>
        <Card style={{ padding: 20, display: "flex", flexDirection: "column", gap: 12 }}>
          {method === "SOR" && (
            <Field label="w (factor de relajación)">
              <input value={w} onChange={e => setW(e.target.value)} />
            </Field>
          )}
          <Field label="Tolerancia"><input value={tol} onChange={e => setTol(e.target.value)} /></Field>
          <Field label="Máx. iteraciones"><input value={maxIter} onChange={e => setMaxIter(e.target.value)} /></Field>
          <Field label="Tipo de error">
            <select value={errType} onChange={e => setErrType(e.target.value)}>
              {["Relativo", "Absoluto", "Condición"].map(e => <option key={e}>{e}</option>)}
            </select>
          </Field>
          <div style={{ marginTop: "auto", display: "flex", flexDirection: "column", gap: 8 }}>
            <Btn onClick={run} color={G.teal} style={{ color: "#000" }}>▶ Ejecutar</Btn>
            <Btn onClick={runComparison} secondary>⇄ Comparar errores</Btn>
          </div>
        </Card>
      </div>

      <ErrBox msg={errMsg} />

      {rho !== null && (
        <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 16 }}>
          <Card style={{ padding: 18, background: rho < 1 ? G.teal + "12" : G.coral + "12", border: `1px solid ${rho < 1 ? G.teal : G.coral}30` }}>
            <div style={{ fontSize: 11, color: G.text2, marginBottom: 4, textTransform: "uppercase", letterSpacing: "0.06em" }}>Radio Espectral ρ(T)</div>
            <div style={{ fontSize: 32, fontWeight: 700, fontFamily: G.mono, color: rho < 1 ? G.teal : G.coral }}>{rho.toFixed(6)}</div>
            <div style={{ fontSize: 13, marginTop: 10, color: rho < 1 ? G.teal : G.coral, fontWeight: 500 }}>
              {rho < 1 ? "✓ CONVERGE — ρ(T) < 1" : "✗ Puede NO CONVERGER — ρ(T) ≥ 1"}
            </div>
            <div style={{ fontSize: 12, color: G.text2, marginTop: 6 }}>
              {rho < 1 ? `Tasa de convergencia ≈ ${(-Math.log10(rho)).toFixed(3)} dígitos/iteración` : "Considera un método directo o reformulación de la matriz."}
            </div>
          </Card>
          {result && (
            <Card style={{ padding: 18 }}>
              <div style={{ fontSize: 11, color: G.text2, marginBottom: 10, textTransform: "uppercase", letterSpacing: "0.06em" }}>Solución — {result.rows.length} iteraciones</div>
              {result.x.map((v, i) => (
                <div key={i} style={{ fontFamily: G.mono, fontSize: 13, color: G.lime, marginBottom: 5 }}>
                  x<sub>{i+1}</sub> = <strong>{v.toFixed(10)}</strong>
                </div>
              ))}
            </Card>
          )}
        </div>
      )}

      {result?.rows?.length > 0 && (
        <Card>
          <div style={{ padding: "14px 20px", borderBottom: `1px solid ${G.border}`, fontSize: 13, fontWeight: 600 }}>
            Iteraciones — {method}
          </div>
          <TableView rows={result.rows} cols={cols} />
        </Card>
      )}

      {compData.length > 0 && (
        <Card style={CHART_STYLE}>
          <div style={{ fontSize: 13, color: G.text2, marginBottom: 12, fontWeight: 500 }}>
            Comparación de tipos de error — {method} (escala log)
          </div>
          <ResponsiveContainer width="100%" height={280}>
            <LineChart data={compData}>
              <CartesianGrid strokeDasharray="3 3" stroke={G.border} />
              <XAxis dataKey="iter" stroke={G.text3} tick={{ fontSize: 10, fill: G.text3 }} />
              <YAxis scale="log" stroke={G.text3} tick={{ fontSize: 10, fill: G.text3 }} />
              <Tooltip contentStyle={TOOLTIP_STYLE} />
              <Legend wrapperStyle={{ fontSize: 12 }} />
              <Line type="monotone" dataKey="Relativo" stroke={G.accent} dot={false} strokeWidth={2} connectNulls={false} />
              <Line type="monotone" dataKey="Absoluto" stroke={G.teal} dot={false} strokeWidth={2} connectNulls={false} />
              <Line type="monotone" dataKey="Condición" stroke={G.coral} dot={false} strokeWidth={2} connectNulls={false} />
            </LineChart>
          </ResponsiveContainer>
        </Card>
      )}
    </div>
  );
}

/* ═══════════════════════════════════════════════════════════════
   CHAPTER 3 COMPONENT
═══════════════════════════════════════════════════════════════ */
const C3_METHODS = ["Vandermonde", "Newton", "Lagrange", "Spline Lineal", "Spline Cúbico"];
const VAL_PCTS = [10, 20, 30, 40];

function Chapter3() {
  const [method, setMethod] = useState("Vandermonde");
  const [pointsText, setPointsText] = useState("1, 2\n2, 3\n3, 5\n4, 4\n5, 7\n6, 8");
  const [evalX, setEvalX] = useState("3.5");
  const [result, setResult] = useState(null);
  const [plotData, setPlotData] = useState([]);
  const [polyStr, setPolyStr] = useState("");
  const [valErrors, setValErrors] = useState([]);
  const [errMsg, setErrMsg] = useState("");

  const parsePoints = () => {
    const lines = pointsText.trim().split("\n")
      .filter(l => l.trim())
      .map(l => l.trim().split(/[\s,;]+/).map(Number));
    const xs = lines.map(l => l[0]), ys = lines.map(l => l[1]);
    if (xs.some(isNaN) || ys.some(isNaN)) throw new Error("Formato inválido. Usa: x, y (una por línea).");
    if (xs.length < 2) throw new Error("Se necesitan al menos 2 puntos.");
    if (xs.length > 10) throw new Error("Máximo 10 puntos.");
    const unique = new Set(xs);
    if (unique.size !== xs.length) throw new Error("Los valores de x deben ser únicos.");
    // Sort by x
    const pts = xs.map((x, i) => ({ x, y: ys[i] })).sort((a, b) => a.x - b.x);
    return { xs: pts.map(p => p.x), ys: pts.map(p => p.y) };
  };

  const run = () => {
    setErrMsg("");
    try {
      const { xs, ys } = parsePoints();
      const xEval = parseFloat(evalX);
      let yEval, extra = {}, poly;

      if (method === "Vandermonde") {
        extra.vanCoeffs = vandermonde(xs, ys);
        yEval = evalVandermonde(extra.vanCoeffs, xEval);
        poly = "P(x) = " + extra.vanCoeffs.map((c, i) => {
          const s = i === 0 ? c.toFixed(5) : `${c >= 0 ? "+" : ""}${c.toFixed(5)}·x^${i}`;
          return s;
        }).join(" ");
      } else if (method === "Newton") {
        extra.newtonCoeffs = dividedDiffs(xs, ys);
        yEval = evalNewton(xs, extra.newtonCoeffs, xEval);
        poly = "Coeficientes de diferencias divididas: [" + extra.newtonCoeffs.map(c => c.toFixed(5)).join(", ") + "]";
      } else if (method === "Lagrange") {
        yEval = evalLagrange(xs, ys, xEval);
        poly = `L(${xEval}) = ${yEval.toFixed(8)} (Forma de Lagrange, ${xs.length - 1}º grado)`;
      } else if (method === "Spline Lineal") {
        yEval = splineLinear(xs, ys, xEval);
        poly = `Spline lineal por tramos (${xs.length - 1} segmentos)`;
      } else {
        extra.spline = buildCubicSpline(xs, ys);
        yEval = evalCubicSpline(extra.spline, xEval);
        poly = `Spline cúbico natural (${xs.length - 1} polinomios cúbicos)`;
      }

      setPolyStr(poly);
      setResult({ yEval, xs, ys });
      setPlotData(genInterp(xs, ys, method, extra));
      computeValidation(xs, ys);
    } catch (e) { setErrMsg(e.message); }
  };

  const computeValidation = (xs, ys) => {
    const errors = VAL_PCTS.map(pct => {
      const n = xs.length;
      const nVal = Math.max(1, Math.round(n * pct / 100));
      const step = Math.floor(n / nVal);
      const idxVal = Array.from({ length: nVal }, (_, i) => Math.min(i * step, n - 1));
      const idxTrain = Array.from({ length: n }, (_, i) => i).filter(i => !idxVal.includes(i));
      if (idxTrain.length < 2) return { pct, error: null, label: `${pct}%` };
      const xTr = idxTrain.map(i => xs[i]), yTr = idxTrain.map(i => ys[i]);
      const xVal = idxVal.map(i => xs[i]), yVal = idxVal.map(i => ys[i]);
      try {
        const coeffs = vandermonde(xTr, yTr);
        const errs = xVal.map((x, i) => Math.abs(evalVandermonde(coeffs, x) - yVal[i]));
        const avg = errs.reduce((a, b) => a + b, 0) / errs.length;
        return { pct, error: +avg.toFixed(6), label: `${pct}%` };
      } catch { return { pct, error: null, label: `${pct}%` }; }
    });
    setValErrors(errors);
  };

  return (
    <div style={{ display: "flex", flexDirection: "column", gap: 22 }}>
      <SectionTitle title="Interpolación" sub="Ajusta un polinomio a un conjunto de puntos (hasta 10)" color={G.coral} />
      <TabBar items={C3_METHODS} active={method} setActive={setMethod} color={G.coral} />

      <InfoBox color={G.coral} icon="📍">
        Ingresa hasta <strong>10 puntos</strong> en formato <code>x, y</code> — uno por línea.<br />
        Los puntos se ordenan automáticamente por x. Los valores de x deben ser únicos.<br />
        <span style={{ opacity: 0.8 }}>Ejemplo: <code>1, 2.5</code> → <code>2, 4.1</code> → <code>3, 3.8</code></span><br />
        {method === "Punto Fijo" ? "" : <span style={{ opacity: 0.7 }}>La gráfica de validación usa Vandermonde sobre el subconjunto de entrenamiento.</span>}
      </InfoBox>

      <div style={{ display: "grid", gridTemplateColumns: "2fr 1fr", gap: 16 }}>
        <Field label={`Puntos (x, y — uno por línea, máx 10)`}>
          <textarea value={pointsText} onChange={e => setPointsText(e.target.value)} rows={9}
            style={{ resize: "vertical", fontFamily: G.mono, fontSize: 13 }} />
        </Field>
        <Card style={{ padding: 20, display: "flex", flexDirection: "column", gap: 14 }}>
          <Field label="Evaluar en x =">
            <input value={evalX} onChange={e => setEvalX(e.target.value)} />
          </Field>
          <div style={{ marginTop: "auto" }}>
            <Btn onClick={run} color={G.coral} style={{ width: "100%" }}>▶ Interpolar</Btn>
          </div>
        </Card>
      </div>

      <ErrBox msg={errMsg} />

      {polyStr && (
        <Card style={{ padding: 18, background: G.coral + "10", border: `1px solid ${G.coral}28` }}>
          <div style={{ fontSize: 11, color: G.text2, marginBottom: 6, textTransform: "uppercase", letterSpacing: "0.06em" }}>Polinomio — {method}</div>
          <div style={{ fontFamily: G.mono, fontSize: 12, color: G.coral + "dd", wordBreak: "break-all", lineHeight: 1.7 }}>{polyStr}</div>
          {result && (
            <div style={{ marginTop: 12, padding: "8px 14px", background: G.bg3, borderRadius: 8, fontFamily: G.mono, fontSize: 14, color: G.coral }}>
              P({evalX}) = <strong>{result.yEval.toFixed(10)}</strong>
            </div>
          )}
        </Card>
      )}

      {plotData.length > 0 && (
        <Card style={CHART_STYLE}>
          <div style={{ fontSize: 13, color: G.text2, marginBottom: 12, fontWeight: 500 }}>Gráfica de interpolación — {method}</div>
          <ResponsiveContainer width="100%" height={260}>
            <LineChart data={plotData}>
              <CartesianGrid strokeDasharray="3 3" stroke={G.border} />
              <XAxis dataKey="x" stroke={G.text3} tick={{ fontSize: 10, fill: G.text3 }} />
              <YAxis stroke={G.text3} tick={{ fontSize: 10, fill: G.text3 }} />
              <Tooltip contentStyle={TOOLTIP_STYLE} />
              <Line type="monotone" dataKey="y" stroke={G.coral} dot={false} strokeWidth={2.5} name="P(x)" />
            </LineChart>
          </ResponsiveContainer>
          {result && (
            <div style={{ marginTop: 12, display: "flex", gap: 8, flexWrap: "wrap" }}>
              {result.xs.map((x, i) => (
                <div key={i} style={{ background: G.bg3, borderRadius: 6, padding: "4px 10px", fontSize: 12, fontFamily: G.mono, color: G.coral + "cc" }}>
                  ({x}, {result.ys[i]})
                </div>
              ))}
            </div>
          )}
        </Card>
      )}

      {valErrors.length > 0 && (
        <Card style={CHART_STYLE}>
          <div style={{ fontSize: 13, color: G.text2, marginBottom: 12, fontWeight: 500 }}>
            Error de validación según porcentaje reservado (Vandermonde)
          </div>
          <ResponsiveContainer width="100%" height={200}>
            <LineChart data={valErrors}>
              <CartesianGrid strokeDasharray="3 3" stroke={G.border} />
              <XAxis dataKey="label" stroke={G.text3} tick={{ fontSize: 11, fill: G.text3 }} />
              <YAxis stroke={G.text3} tick={{ fontSize: 11, fill: G.text3 }} />
              <Tooltip contentStyle={TOOLTIP_STYLE} formatter={v => [v, "Error promedio abs."]} />
              <Line type="monotone" dataKey="error" stroke={G.coral} strokeWidth={2.5} dot={{ fill: G.coral, r: 5 }} />
            </LineChart>
          </ResponsiveContainer>
          <div style={{ display: "grid", gridTemplateColumns: "repeat(4,1fr)", gap: 10, marginTop: 14 }}>
            {valErrors.map(({ pct, error }) => (
              <div key={pct} style={{ background: G.bg3, borderRadius: 8, padding: "10px 12px", textAlign: "center" }}>
                <div style={{ fontSize: 11, color: G.text2 }}>{pct}% validación</div>
                <div style={{ fontFamily: G.mono, fontSize: 14, color: error === null ? G.text3 : G.coral, marginTop: 4 }}>
                  {error === null ? "N/A" : error}
                </div>
              </div>
            ))}
          </div>
        </Card>
      )}
    </div>
  );
}

/* ═══════════════════════════════════════════════════════════════
   APP ROOT
═══════════════════════════════════════════════════════════════ */
const CHAPTERS = [
  { id: 1, label: "Cap. 1", title: "Raíces", fullTitle: "Raíces de Ecuaciones", icon: "√x", color: G.accent },
  { id: 2, label: "Cap. 2", title: "Sistemas", fullTitle: "Sistemas Lineales", icon: "Ax=b", color: G.teal },
  { id: 3, label: "Cap. 3", title: "Interpolación", fullTitle: "Interpolación", icon: "P(x)", color: G.coral },
];

export default function App() {
  const [active, setActive] = useState(1);
  const ch = CHAPTERS.find(c => c.id === active);

  return (
    <>
      <style>{css}</style>
      <div style={{ minHeight: "100vh", display: "flex", flexDirection: "column", background: G.bg }}>
        {/* Header */}
        <header style={{
          background: G.bg2, borderBottom: `1px solid ${G.border}`,
          padding: "0 2rem", display: "flex", alignItems: "center",
          gap: "1.5rem", height: 58, position: "sticky", top: 0, zIndex: 100
        }}>
          <div style={{ display: "flex", alignItems: "center", gap: 10 }}>
            <div style={{
              width: 34, height: 34, borderRadius: 9,
              background: `linear-gradient(135deg, ${G.accent}, ${G.teal})`,
              display: "flex", alignItems: "center", justifyContent: "center",
              fontSize: 15, fontWeight: 700, fontFamily: G.mono, color: "#fff"
            }}>∑</div>
            <span style={{ fontWeight: 600, fontSize: 15, letterSpacing: "-0.3px" }}>
              Métodos Numéricos
            </span>
          </div>
          <nav style={{ display: "flex", gap: 4 }}>
            {CHAPTERS.map(ch => (
              <button key={ch.id} onClick={() => setActive(ch.id)} style={{
                background: active === ch.id ? ch.color + "1a" : "transparent",
                color: active === ch.id ? ch.color : G.text2,
                border: active === ch.id ? `1px solid ${ch.color}35` : "1px solid transparent",
                borderRadius: 7, padding: "5px 14px", cursor: "pointer",
                fontSize: 13, fontWeight: active === ch.id ? 600 : 400,
                transition: "all 0.15s", fontFamily: G.sans, display: "flex", alignItems: "center", gap: 7
              }}>
                <span style={{ fontFamily: G.mono, fontSize: 11, opacity: 0.7 }}>{ch.label}</span>
                {ch.fullTitle}
              </button>
            ))}
          </nav>
          <div style={{ marginLeft: "auto", fontSize: 11, color: G.text3, fontFamily: G.mono }}>
            {ch.icon}
          </div>
        </header>

        {/* Content */}
        <main style={{ flex: 1, padding: "28px 32px", maxWidth: 1100, margin: "0 auto", width: "100%" }}>
          {active === 1 && <Chapter1 />}
          {active === 2 && <Chapter2 />}
          {active === 3 && <Chapter3 />}
        </main>
      </div>
    </>
  );
}