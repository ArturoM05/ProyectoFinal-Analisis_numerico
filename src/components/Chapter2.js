import { useState } from "react";
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from "recharts";

const METHODS = ["Jacobi", "Gauss-Seidel", "SOR"];
const ERROR_TYPES = ["Relativo", "Absoluto", "Condición"];
const COLORS = ["#7c6af7", "#4ecdc4", "#ff6b6b"];

function parseMatrix(text) {
  return text.trim().split("\n").map(row => row.trim().split(/[\s,]+/).map(Number));
}

function parseVector(text) {
  return text.trim().split(/[\s,]+/).map(Number);
}

// Calcula los 3 tipos de error para un par (x, xOld)
function calcAllErrors(x, xOld, A, b) {
  // Absoluto: max diferencia entre iteraciones
  const absoluto = Math.max(...x.map((v, i) => Math.abs(v - xOld[i])));

  // Relativo: absoluto / max(|x|)
  const maxX = Math.max(...x.map(v => Math.abs(v)));
  const relativo = maxX === 0 ? absoluto : absoluto / maxX;

  // Condición: norma del residual ||b - Ax||
  const res = A.map((row, i) => b[i] - row.reduce((s, v, j) => s + v * x[j], 0));
  const condicion = Math.sqrt(res.reduce((s, v) => s + v ** 2, 0));

  return {
    Relativo:  isFinite(relativo)  && relativo  > 0 ? relativo  : null,
    Absoluto:  isFinite(absoluto)  && absoluto  > 0 ? absoluto  : null,
    Condicion: isFinite(condicion) && condicion > 0 ? condicion : null,
  };
}

// Error simple para criterio de parada
function calcError(x, xOld, A, b, errType) {
  if (errType === "Relativo") {
    const maxX = Math.max(...x.map(v => Math.abs(v)));
    const delta = Math.max(...x.map((v, i) => Math.abs(v - xOld[i])));
    return maxX === 0 ? delta : delta / maxX;
  }
  if (errType === "Absoluto") {
    return Math.max(...x.map((v, i) => Math.abs(v - xOld[i])));
  }
  const res = A.map((row, i) => b[i] - row.reduce((s, v, j) => s + v * x[j], 0));
  return Math.sqrt(res.reduce((s, v) => s + v ** 2, 0));
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
    const allErr = calcAllErrors(x, xOld, A, b);
    rows.push({
      iter: k,
      ...Object.fromEntries(x.map((v, i) => [`x${i+1}`, +v.toFixed(6)])),
      error: isFinite(err) ? +err.toFixed(8) : "---",
      ...allErr
    });
    if (isFinite(err) && err < tol) break;
  }
  return { rows, x };
}

function spectralRadius(A, method, w) {
  const n = A.length;
  const D = A.map((r, i) => r.map((v, j) => i === j ? v : 0));
  const L = A.map((r, i) => r.map((v, j) => j < i ? -v : 0));
  const U = A.map((r, i) => r.map((v, j) => j > i ? -v : 0));
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
  return A.map(r => B[0].map((_, j) => r.reduce((s, v, k) => s + v * B[k][j], 0)));
}
function matInv(M) {
  const n = M.length;
  const aug = M.map((r, i) => [...r, ...Array.from({ length: n }, (_, j) => i === j ? 1 : 0)]);
  for (let i = 0; i < n; i++) {
    let pivot = aug[i][i];
    if (Math.abs(pivot) < 1e-12) {
      for (let k = i+1; k < n; k++) {
        if (Math.abs(aug[k][i]) > Math.abs(pivot)) { [aug[i], aug[k]] = [aug[k], aug[i]]; pivot = aug[i][i]; break; }
      }
    }
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
      setRho(spectralRadius(A, method, wVal));
      setCompData([]); // limpiar comparativa
    } catch (e) { setError("Error: " + e.message); }
  };

  // Comparativa: 3 líneas (Relativo, Absoluto, Condición) por iteración del método actual
  const runComparison = () => {
    setError("");
    try {
      const A = parseMatrix(matA);
      const b = parseMatrix(vecB).flat();
      const x0Vec = parseVector(x0);
      const n = A.length;
      if (A.some(r => r.length !== n) || b.length !== n) throw new Error("Dimensiones incorrectas");
      const t = parseFloat(tol), m = parseInt(maxIter), wVal = parseFloat(w);
      const { rows } = runIterative(A, b, method, wVal, x0Vec, t, m, errType);

      // Saltar iter 0 (sin error previo)
      const data = rows
        .filter(row => row.iter > 0)
        .map(row => ({
          iter: row.iter,
          Relativo:  row.Relativo  !== null && row.Relativo  > 0 ? row.Relativo  : null,
          Absoluto:  row.Absoluto  !== null && row.Absoluto  > 0 ? row.Absoluto  : null,
          Condicion: row.Condicion !== null && row.Condicion > 0 ? row.Condicion : null,
        }));

      setCompData(data);
    } catch (e) { setError("Error en comparación: " + e.message); }
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
          <button key={m} onClick={() => { setMethod(m); setResult(null); setCompData([]); }} style={{
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
            <input value={x0} onChange={e => setX0(e.target.value)} placeholder="0 0 0" />
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
            <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>Tipo de error (criterio parada)</span>
            <select value={errType} onChange={e => setErrType(e.target.value)}>
              {ERROR_TYPES.map(e => <option key={e}>{e}</option>)}
            </select>
          </label>
          <button onClick={run} style={{ marginTop: "auto", background: "#4ecdc4", color: "#000" }}>▶ Ejecutar</button>
          <button onClick={runComparison} className="secondary">⇄ Comparar errores</button>
        </div>
      </div>

      {error && <div style={{ background: "#ff6b6b22", border: "1px solid #ff6b6b44", borderRadius: 8, padding: 12, color: "#ff6b6b", fontSize: 13 }}>{error}</div>}

      {/* Radio espectral */}
      {rho !== null && (
        <div style={{ display: "flex", gap: 16 }}>
          <div style={{
            flex: 1,
            background: rho < 1 ? "#4ecdc422" : "#ff6b6b22",
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
              <div style={{ fontSize: 12, color: "var(--text2)", marginBottom: 8 }}>Solución — {result.rows.length} iteraciones</div>
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

      {/* Comparison — 3 líneas de error por iteración */}
      {compData.length > 0 && (
        <div style={{ background: "var(--bg2)", borderRadius: 12, padding: 20, border: "1px solid var(--border)" }}>
          <h3 style={{ fontSize: 14, color: "var(--text2)", marginBottom: 4 }}>
            Comparación de tipos de error — <span style={{ color: "#4ecdc4" }}>{method}</span>
          </h3>
          <p style={{ fontSize: 12, color: "var(--text3)", marginBottom: 16 }}>
            Evolución del error relativo, absoluto y de condición (residual) a través de las iteraciones
          </p>
          <ResponsiveContainer width="100%" height={300}>
            <LineChart data={compData} margin={{ top: 5, right: 20, bottom: 5, left: 10 }}>
              <CartesianGrid strokeDasharray="3 3" stroke="#2a2a38" />
              <XAxis dataKey="iter" stroke="#5a5a78" tick={{ fontSize: 11 }}
                label={{ value: "Iteración", position: "insideBottom", offset: -2, fill: "#5a5a78", fontSize: 11 }} />
              <YAxis
                scale="log"
                domain={['auto', 'auto']}
                stroke="#5a5a78"
                tick={{ fontSize: 10 }}
                tickFormatter={v => v > 0 ? v.toExponential(0) : ""}
                allowDataOverflow
              />
              <Tooltip
                contentStyle={{ background: "#111118", border: "1px solid #2a2a38", borderRadius: 8 }}
                formatter={(value, name) => [value !== null ? value.toExponential(4) : "N/A", name]}
                labelFormatter={l => `Iteración ${l}`}
              />
              <Legend wrapperStyle={{ fontSize: 12, paddingTop: 12 }} />
              <Line type="monotone" dataKey="Relativo"  stroke="#7c6af7" dot={false} strokeWidth={2} connectNulls={false} />
              <Line type="monotone" dataKey="Absoluto"  stroke="#4ecdc4" dot={false} strokeWidth={2} connectNulls={false} />
              <Line type="monotone" dataKey="Condicion" stroke="#ff6b6b" dot={false} strokeWidth={2} connectNulls={false} strokeDasharray="5 3" />
            </LineChart>
          </ResponsiveContainer>

          {/* Resumen valores finales */}
          <div style={{ display: "flex", gap: 12, marginTop: 16 }}>
            {[
              { label: "Error Relativo final",        key: "Relativo",  color: "#7c6af7" },
              { label: "Error Absoluto final",        key: "Absoluto",  color: "#4ecdc4" },
              { label: "Error Condición final (res)", key: "Condicion", color: "#ff6b6b" },
            ].map(({ label, key, color }) => {
              const last = [...compData].reverse().find(d => d[key] !== null);
              return (
                <div key={key} style={{ flex: 1, background: color + "11", border: `1px solid ${color}33`, borderRadius: 8, padding: "10px 14px" }}>
                  <div style={{ fontSize: 11, color: "var(--text2)", marginBottom: 4 }}>{label}</div>
                  <div style={{ fontFamily: "var(--mono)", fontSize: 13, color }}>
                    {last ? last[key].toExponential(4) : "N/A"}
                  </div>
                </div>
              );
            })}
          </div>
        </div>
      )}
    </div>
  );
}