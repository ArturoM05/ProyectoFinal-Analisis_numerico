import { useState } from "react";
import { BarChart, Bar, Cell, LabelList, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from "recharts";

const METHODS = ["Jacobi", "Gauss-Seidel", "SOR"];
const ERROR_TYPES = ["Relativo", "Absoluto", "Condición"];
const COLORS = ["#7c6af7", "#4ecdc4", "#ff6b6b"];

function parseMatrix(text) {
  return text.trim().split("\n").map(row => row.trim().split(/[\s,]+/).map(Number));
}

function parseVector(text) {
  return text.trim().split(/[\s,]+/).map(Number);
}

function runIterative(A, b, method, w, x0, tol, maxIter, errType) {
  const n = A.length;
  let x = x0 ? [...x0] : new Array(n).fill(0);
  const rows = [];
  rows.push({ iter: 0, ...Object.fromEntries(x.map((v, i) => [`x${i+1}`, +v.toFixed(6)])), error: "---" });

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
    const err = calcError(x, xOld, A, b, errType);
    rows.push({ iter: k, ...Object.fromEntries(x.map((v, i) => [`x${i+1}`, +v.toFixed(6)])), error: isFinite(err) ? +err.toFixed(8) : "---" });
    if (isFinite(err) && err < tol) break;
  }
  return { rows, x };
}

function calcError(x, xOld, A, b, errType) {
  if (errType === "Relativo") {
    const maxX = Math.max(...x.map(v => Math.abs(v)));
    const delta = Math.max(...x.map((v, i) => Math.abs(v - xOld[i])));
    return maxX === 0 ? delta : delta / maxX;
  }
  if (errType === "Absoluto") {
    return Math.max(...x.map((v, i) => Math.abs(v - xOld[i])));
  }
  // Condición: residual
  const res = A.map((row, i) => b[i] - row.reduce((s, v, j) => s + v * x[j], 0));
  return Math.sqrt(res.reduce((s, v) => s + v ** 2, 0));
}

function spectralRadius(A, method, w) {
  const n = A.length;
  const D = A.map((r, i) => r.map((v, j) => i === j ? v : 0));
  const L = A.map((r, i) => r.map((v, j) => j < i ? -v : 0));
  const U = A.map((r, i) => r.map((v, j) => j > i ? -v : 0));

  // T matrix
  let T;
  if (method === "Jacobi") {
    const Dinv = D.map((r, i) => r.map((v, j) => i === j ? 1 / v : 0));
    T = matMul(Dinv, matAdd(L, U));
  } else if (method === "Gauss-Seidel") {
    const M = matAdd(D, L.map(r => r.map(v => -v)));
    const Minv = matInv(M);
    T = matMul(Minv, U.map(r => r.map(v => -v)));
  } else {
    const DwL = matAdd(D, L.map(r => r.map(v => -w * v)));
    const Minv = matInv(DwL);
    const rhs = matAdd(D.map(r => r.map(v => (1-w)*v)), U.map(r => r.map(v => -w*v)));
    T = matMul(Minv, rhs);
  }
  return maxEigenvalue(T);
}

function matAdd(A, B) { return A.map((r, i) => r.map((v, j) => v + B[i][j])); }
function matMul(A, B) {
  const n = A.length;
  return A.map(r => B[0].map((_, j) => r.reduce((s, v, k) => s + v * B[k][j], 0)));
}
function matInv(M) {
  const n = M.length;
  const aug = M.map((r, i) => [...r, ...Array.from({ length: n }, (_, j) => i === j ? 1 : 0)]);
  for (let i = 0; i < n; i++) {
    let pivot = aug[i][i];
    if (Math.abs(pivot) < 1e-12) { for (let k = i+1; k < n; k++) if (Math.abs(aug[k][i]) > Math.abs(pivot)) { [aug[i], aug[k]] = [aug[k], aug[i]]; pivot = aug[i][i]; break; } }
    aug[i] = aug[i].map(v => v / pivot);
    for (let k = 0; k < n; k++) if (k !== i) { const f = aug[k][i]; aug[k] = aug[k].map((v, j) => v - f * aug[i][j]); }
  }
  return aug.map(r => r.slice(n));
}
function maxEigenvalue(M) {
  const n = M.length;
  let v = Array.from({ length: n }, () => Math.random());
  let lambda = 0;
  for (let iter = 0; iter < 1000; iter++) {
    const Mv = M.map(r => r.reduce((s, v2, j) => s + v2 * v[j], 0));
    const norm = Math.sqrt(Mv.reduce((s, x) => s + x*x, 0));
    if (norm < 1e-14) return 0;
    lambda = norm;
    v = Mv.map(x => x / norm);
  }
  return lambda;
}

const DEFAULT_A = "4 -1 0\n-1 4 -1\n0 -1 4";
const DEFAULT_B = "15\n10\n10";

export default function Chapter2() {
  const [method, setMethod] = useState("Gauss-Seidel");
  const [size, setSize] = useState(3);
  const [matA, setMatA] = useState(DEFAULT_A);
  const [vecB, setVecB] = useState(DEFAULT_B);
  const [x0, setX0] = useState("0 0 0");
  const [w, setW] = useState("1.1");
  const [tol, setTol] = useState("1e-6");
  const [maxIter, setMaxIter] = useState("100");
  const [errType, setErrType] = useState("Absoluto");
  const [result, setResult] = useState(null);
  const [rho, setRho] = useState(null);
  const [compData, setCompData] = useState([]);
  const [error, setError] = useState("");

  const run = () => {
    setError("");
    try {
      const A = parseMatrix(matA);
      const b = parseMatrix(vecB).flat();
      const x0Vec = parseVector(x0);
      const n = A.length;
      if (A.some(r => r.length !== n) || b.length !== n) throw new Error("Dimensiones incorrectas");
      if (x0Vec.length !== n || x0Vec.some(isNaN)) throw new Error("x0 debe tener " + n + " valores numéricos");
      const t = parseFloat(tol), m = parseInt(maxIter), wVal = parseFloat(w);
      const res = runIterative(A, b, method, wVal, x0Vec, t, m, errType);
      setResult({ ...res, n });
      const rhoVal = spectralRadius(A, method, wVal);
      setRho(rhoVal);
    } catch (e) { setError("Error: " + e.message); }
  };

  const runComparison = () => {
    try {
      const A = parseMatrix(matA);
      const b = parseMatrix(vecB).flat();
      const t = parseFloat(tol), m = parseInt(maxIter), wVal = parseFloat(w);
      const x0Vec = parseVector(x0);
      const methods = [
        { name: "Jacobi", res: runIterative(A, b, "Jacobi", 1, x0Vec, t, m, errType) },
        { name: "Gauss-Seidel", res: runIterative(A, b, "Gauss-Seidel", 1, x0Vec, t, m, errType) },
        { name: "SOR", res: runIterative(A, b, "SOR", wVal, x0Vec, t, m, errType) },
      ];
      const comp = methods.map(({ name, res }) => {
        const last = res.rows[res.rows.length - 1];
        const finalError = last && typeof last.error === "number" && isFinite(last.error) ? last.error : null;
        return { name, error: finalError };
      });
      console.log("Chapter2 comparison data:", comp);
      setCompData(comp);
    } catch (e) { setError("Error comparación: " + e.message); }
  };

  const cols = result ? ["iter", ...Array.from({ length: result.n }, (_, i) => `x${i+1}`), "error"] : [];

  return (
    <div style={{ display: "flex", flexDirection: "column", gap: 24 }}>
      <div>
        <h1 style={{ fontSize: 24, fontWeight: 600, color: "#4ecdc4" }}>Sistemas de Ecuaciones Lineales</h1>
        <p style={{ color: "var(--text2)", fontSize: 14, marginTop: 4 }}>Métodos iterativos para resolver Ax = b</p>
      </div>

      {/* Method selector */}
      <div style={{ display: "flex", gap: 8 }}>
        {METHODS.map(m => (
          <button key={m} onClick={() => setMethod(m)} style={{
            background: method === m ? "#4ecdc422" : "var(--bg3)",
            color: method === m ? "#4ecdc4" : "var(--text2)",
            border: method === m ? "1px solid #4ecdc444" : "1px solid var(--border)",
            borderRadius: 8, padding: "6px 20px", cursor: "pointer",
            fontSize: 13, fontFamily: "var(--sans)"
          }}>{m}</button>
        ))}
      </div>

      {/* Help */}
      <div style={{ background: "#4ecdc411", border: "1px solid #4ecdc433", borderRadius: 10, padding: "12px 16px", fontSize: 13, color: "#7eddd8" }}>
        💡 Ingresa la matriz A fila por fila (valores separados por espacios) y el vector b, uno por línea.
        Tamaño máximo: 8×8.
        <br />Para el método de Jacobi usa <strong>Error Absoluto</strong> si quieres el "mayor delta" entre iteraciones como en la tabla del profe.
        {method === "SOR" && " Para SOR se recomienda w ∈ (1, 2) para acelerar convergencia."}
      </div>

      {/* Inputs */}
      <div style={{ display: "grid", gridTemplateColumns: "2fr 1fr", gap: 16 }}>
        <div style={{ background: "var(--bg2)", borderRadius: 12, padding: 20, border: "1px solid var(--border)", display: "flex", flexDirection: "column", gap: 16 }}>
          <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
            <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>Matriz A (una fila por línea)</span>
            <textarea value={matA} onChange={e => setMatA(e.target.value)} rows={5} style={{ resize: "vertical", fontFamily: "var(--mono)", fontSize: 13 }} />
          </label>
          <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
            <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>Vector b (un valor por línea)</span>
            <textarea value={vecB} onChange={e => setVecB(e.target.value)} rows={3} style={{ resize: "vertical", fontFamily: "var(--mono)", fontSize: 13 }} />
          </label>
          <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
            <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>Vector inicial x₀</span>
            <input value={x0} onChange={e => setX0(e.target.value)} placeholder="2 2 2 2" />
          </label>
        </div>
        <div style={{ background: "var(--bg2)", borderRadius: 12, padding: 20, border: "1px solid var(--border)", display: "flex", flexDirection: "column", gap: 12 }}>
          {method === "SOR" && (
            <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
              <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>w (relajación)</span>
              <input value={w} onChange={e => setW(e.target.value)} />
            </label>
          )}
          <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
            <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>Tolerancia</span>
            <input value={tol} onChange={e => setTol(e.target.value)} />
          </label>
          <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
            <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>Máx. iteraciones</span>
            <input value={maxIter} onChange={e => setMaxIter(e.target.value)} />
          </label>
          <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
            <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>Tipo de error</span>
            <select value={errType} onChange={e => setErrType(e.target.value)}>
              {ERROR_TYPES.map(e => <option key={e}>{e}</option>)}
            </select>
          </label>
          <button onClick={run} style={{ marginTop: "auto", background: "#4ecdc4", color: "#000" }}>▶ Ejecutar</button>
          <button onClick={runComparison} className="secondary">⇄ Comparar</button>
        </div>
      </div>

      {error && <div style={{ background: "#ff6b6b22", border: "1px solid #ff6b6b44", borderRadius: 8, padding: 12, color: "#ff6b6b", fontSize: 13 }}>{error}</div>}

      {/* Radio espectral */}
      {rho !== null && (
        <div style={{ display: "flex", gap: 16 }}>
          <div style={{
            flex: 1, background: rho < 1 ? "#4ecdc422" : "#ff6b6b22",
            border: `1px solid ${rho < 1 ? "#4ecdc444" : "#ff6b6b44"}`,
            borderRadius: 10, padding: 16
          }}>
            <div style={{ fontSize: 12, color: "var(--text2)", marginBottom: 4 }}>Radio Espectral ρ(T)</div>
            <div style={{ fontSize: 28, fontWeight: 600, fontFamily: "var(--mono)", color: rho < 1 ? "#4ecdc4" : "#ff6b6b" }}>{rho.toFixed(6)}</div>
            <div style={{ fontSize: 13, marginTop: 8, color: rho < 1 ? "#4ecdc4" : "#ff6b6b" }}>
              {rho < 1 ? "✓ El método CONVERGE (ρ < 1)" : "✗ El método puede NO converger (ρ ≥ 1)"}
            </div>
          </div>
          {result && (
            <div style={{ flex: 1, background: "var(--bg2)", border: "1px solid var(--border)", borderRadius: 10, padding: 16 }}>
              <div style={{ fontSize: 12, color: "var(--text2)", marginBottom: 8 }}>Solución encontrada en {result.rows.length} iteraciones</div>
              {result.x.map((v, i) => (
                <div key={i} style={{ fontFamily: "var(--mono)", fontSize: 14, color: "#4ecdc4", marginBottom: 4 }}>
                  x<sub>{i+1}</sub> = {v.toFixed(8)}
                </div>
              ))}
            </div>
          )}
        </div>
      )}

      {/* Table */}
      {result && (
        <div style={{ background: "var(--bg2)", borderRadius: 12, border: "1px solid var(--border)", overflow: "hidden" }}>
          <div style={{ padding: "16px 20px", borderBottom: "1px solid var(--border)" }}>
            <h3 style={{ fontSize: 14, fontWeight: 500 }}>Tabla de iteraciones — {method}</h3>
          </div>
          <div style={{ overflowX: "auto" }}>
            <table>
              <thead><tr>{cols.map(c => <th key={c}>{c}</th>)}</tr></thead>
              <tbody>{result.rows.map((row, i) => (
                <tr key={i}>{cols.map(c => <td key={c}>{row[c] ?? "---"}</td>)}</tr>
              ))}</tbody>
            </table>
          </div>
        </div>
      )}

      {/* Comparison */}
      {compData.length > 0 && (
        <div style={{ background: "var(--bg2)", borderRadius: 12, padding: 20, border: "1px solid var(--border)" }}>
          <h3 style={{ fontSize: 14, color: "var(--text2)", marginBottom: 16 }}>Comparación de métodos — error final</h3>
          <ResponsiveContainer width="100%" height={280}>
            <BarChart data={compData}>
              <CartesianGrid strokeDasharray="3 3" stroke="#2a2a38" />
              <XAxis dataKey="name" stroke="#5a5a78" tick={{ fontSize: 11 }} />
              <YAxis type="number" scale="log" domain={[dataMin => Math.max(dataMin || 1e-16, 1e-16), 'dataMax']} stroke="#5a5a78" tick={{ fontSize: 11 }} />
              <Tooltip contentStyle={{ background: "#111118", border: "1px solid #2a2a38", borderRadius: 8 }} formatter={value => typeof value === 'number' ? value.toExponential(3) : value} />
              <Bar dataKey="error">
                {compData.map((entry, index) => (
                  <Cell key={entry.name} fill={COLORS[index % COLORS.length]} />
                ))}
                <LabelList dataKey="error" position="top" formatter={value => typeof value === 'number' ? value.toExponential(2) : ''} />
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        </div>
      )}
    </div>
  );
}
