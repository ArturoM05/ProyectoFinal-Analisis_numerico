import { useState } from "react";
import { LineChart, Line, ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from "recharts";

const METHODS = ["Vandermonde", "Newton", "Lagrange", "Spline Lineal", "Spline Cúbico"];
const VAL_PCTS = [10, 20, 30, 40];
const COLORS = ["#7c6af7", "#4ecdc4", "#ff6b6b", "#ffd93d", "#a8ff78"];

// ---- Math helpers ----
function vandermonde(xs, ys) {
  const n = xs.length;
  const V = xs.map(x => Array.from({ length: n }, (_, j) => Math.pow(x, j)));
  return gaussSolve(V, ys);
}

function dividedDiffs(xs, ys) {
  const n = xs.length, table = ys.map(y => [y]);
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

function evalVandermonde(coeffs, x) {
  return coeffs.reduce((s, c, i) => s + c * Math.pow(x, i), 0);
}

function evalLagrange(xs, ys, x) {
  return xs.reduce((sum, xi, i) => {
    const Li = xs.reduce((p, xj, j) => i === j ? p : p * (x - xj) / (xi - xj), 1);
    return sum + ys[i] * Li;
  }, 0);
}

function splineLinear(xs, ys, x) {
  const i = Math.max(0, xs.findIndex((xi, idx) => xi > x) - 1);
  const idx = Math.min(i, xs.length - 2);
  const h = xs[idx+1] - xs[idx];
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
  const i = Math.min(Math.max(xs.findIndex((xi, idx) => xi > x) - 1, 0), xs.length - 2);
  const hi = h[i], t = x - xs[i];
  return ((M[i+1]-M[i])/(6*hi))*t*t*t + (M[i]/2)*t*t + ((ys[i+1]-ys[i])/hi - hi*(2*M[i]+M[i+1])/6)*t + ys[i];
}

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

function generatePlot(xs, ys, method, splineData, vanCoeffs, newtonCoeffs) {
  const xMin = Math.min(...xs), xMax = Math.max(...xs);
  const pts = [];
  for (let i = 0; i <= 300; i++) {
    const x = xMin + (xMax - xMin) * i / 300;
    let y;
    try {
      if (method === "Vandermonde") y = evalVandermonde(vanCoeffs, x);
      else if (method === "Newton") y = evalNewton(xs, newtonCoeffs, x);
      else if (method === "Lagrange") y = evalLagrange(xs, ys, x);
      else if (method === "Spline Lineal") y = splineLinear(xs, ys, x);
      else y = evalCubicSpline(splineData, x);
      if (isFinite(y)) pts.push({ x: +x.toFixed(4), y: +y.toFixed(4) });
    } catch { }
  }
  return pts;
}

export default function Chapter3() {
  const [method, setMethod] = useState("Vandermonde");
  const [pointsText, setPointsText] = useState("1, 3\n2, 5\n3, 4\n4, 7\n5, 9");
  const [evalX, setEvalX] = useState("2.5");
  const [valPct, setValPct] = useState(30);
  const [result, setResult] = useState(null);
  const [plotData, setPlotData] = useState([]);
  const [scatterData, setScatterData] = useState([]);
  const [polynomial, setPolynomial] = useState("");
  const [compErrors, setCompErrors] = useState([]);
  const [error, setError] = useState("");

  const parsePoints = () => {
    const lines = pointsText.trim().split("\n").map(l => l.trim().split(/[\s,]+/).map(Number));
    const xs = lines.map(l => l[0]), ys = lines.map(l => l[1]);
    if (xs.some(isNaN) || ys.some(isNaN)) throw new Error("Formato inválido. Usa: x, y (una por línea)");
    if (xs.length < 2) throw new Error("Se necesitan al menos 2 puntos");
    if (xs.length > 10) throw new Error("Máximo 10 puntos");
    return { xs, ys };
  };

  const run = () => {
    setError("");
    try {
      const { xs, ys } = parsePoints();
      const xEval = parseFloat(evalX);
      let yEval, vanCoeffs, newtonCoeffs, splineData, polyStr;

      if (method === "Vandermonde") {
        vanCoeffs = vandermonde(xs, ys);
        yEval = evalVandermonde(vanCoeffs, xEval);
        polyStr = "P(x) = " + vanCoeffs.map((c, i) => i === 0 ? c.toFixed(4) : `${c >= 0 ? "+" : ""}${c.toFixed(4)}x^${i}`).join(" ");
      } else if (method === "Newton") {
        newtonCoeffs = dividedDiffs(xs, ys);
        yEval = evalNewton(xs, newtonCoeffs, xEval);
        polyStr = "P(x) = " + newtonCoeffs.map((c, i) => {
          if (i === 0) return c.toFixed(4);
          const term = xs.slice(0, i).map(xi => `(x${xi >= 0 ? "-" : "+"}${Math.abs(xi).toFixed(2)})`).join("");
          return `${c >= 0 ? "+" : ""}${c.toFixed(4)}${term}`;
        }).join(" ");
      } else if (method === "Lagrange") {
        yEval = evalLagrange(xs, ys, xEval);
        polyStr = `P(${xEval}) = ${yEval.toFixed(6)} (Lagrange interpolado)`;
      } else if (method === "Spline Lineal") {
        yEval = splineLinear(xs, ys, xEval);
        polyStr = "Spline lineal por tramos";
      } else {
        splineData = buildCubicSpline(xs, ys);
        yEval = evalCubicSpline(splineData, xEval);
        polyStr = `Spline cúbico natural — S(${xEval}) = ${yEval.toFixed(6)}`;
      }

      setPolynomial(polyStr);
      setResult({ yEval, xs, ys });
      const plot = generatePlot(xs, ys, method, splineData, vanCoeffs, newtonCoeffs);
      setPlotData(plot);
      setScatterData(xs.map((x, i) => ({ x, y: ys[i] })));

      // Validation split
      runValidation(xs, ys, valPct);
    } catch (e) { setError("Error: " + e.message); }
  };

  const runValidation = (xs, ys, pct) => {
    try {
      const n = xs.length;
      const errors = VAL_PCTS.map(p => {
        const nVal = Math.max(1, Math.round(n * p / 100));
        const step = Math.floor(n / nVal);
        const idxVal = Array.from({ length: nVal }, (_, i) => Math.min(i * step, n - 1));
        const idxTrain = Array.from({ length: n }, (_, i) => i).filter(i => !idxVal.includes(i));
        if (idxTrain.length < 2) return { pct: p, error: null };
        const xTr = idxTrain.map(i => xs[i]), yTr = idxTrain.map(i => ys[i]);
        const xVal = idxVal.map(i => xs[i]), yVal = idxVal.map(i => ys[i]);
        try {
          const coeffs = vandermonde(xTr, yTr);
          const errs = xVal.map((x, i) => Math.abs(evalVandermonde(coeffs, x) - yVal[i]));
          return { pct: p, error: +(errs.reduce((a, b) => a + b, 0) / errs.length).toFixed(4) };
        } catch { return { pct: p, error: null }; }
      });
      setCompErrors(errors);
    } catch { }
  };

  return (
    <div style={{ display: "flex", flexDirection: "column", gap: 24 }}>
      <div>
        <h1 style={{ fontSize: 24, fontWeight: 600, color: "#ff6b6b" }}>Interpolación</h1>
        <p style={{ color: "var(--text2)", fontSize: 14, marginTop: 4 }}>Ajusta un polinomio a un conjunto de puntos</p>
      </div>

      {/* Method */}
      <div style={{ display: "flex", gap: 8, flexWrap: "wrap" }}>
        {METHODS.map(m => (
          <button key={m} onClick={() => setMethod(m)} style={{
            background: method === m ? "#ff6b6b22" : "var(--bg3)",
            color: method === m ? "#ff6b6b" : "var(--text2)",
            border: method === m ? "1px solid #ff6b6b44" : "1px solid var(--border)",
            borderRadius: 8, padding: "6px 14px", cursor: "pointer",
            fontSize: 13, fontFamily: "var(--sans)"
          }}>{m}</button>
        ))}
      </div>

      {/* Help */}
      <div style={{ background: "#ff6b6b11", border: "1px solid #ff6b6b33", borderRadius: 10, padding: "12px 16px", fontSize: 13, color: "#ff9999" }}>
        💡 Ingresa hasta 10 puntos en formato <code>x, y</code> (uno por línea). Ejemplo: <code>1, 3.5</code>
      </div>

      {/* Inputs */}
      <div style={{ display: "grid", gridTemplateColumns: "2fr 1fr", gap: 16 }}>
        <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
          <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>Puntos (x, y — uno por línea, máx 10)</span>
          <textarea value={pointsText} onChange={e => setPointsText(e.target.value)} rows={8}
            style={{ resize: "vertical", fontFamily: "var(--mono)", fontSize: 13 }} />
        </label>
        <div style={{ background: "var(--bg2)", borderRadius: 12, padding: 20, border: "1px solid var(--border)", display: "flex", flexDirection: "column", gap: 12 }}>
          <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
            <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>Evaluar en x =</span>
            <input value={evalX} onChange={e => setEvalX(e.target.value)} />
          </label>
          <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
            <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>% validación</span>
            <select value={valPct} onChange={e => setValPct(parseInt(e.target.value))}>
              {VAL_PCTS.map(p => <option key={p} value={p}>{p}%</option>)}
            </select>
          </label>
          <button onClick={run} style={{ marginTop: "auto", background: "#ff6b6b", color: "white" }}>▶ Interpolar</button>
        </div>
      </div>

      {error && <div style={{ background: "#ff6b6b22", border: "1px solid #ff6b6b44", borderRadius: 8, padding: 12, color: "#ff6b6b", fontSize: 13 }}>{error}</div>}

      {/* Polynomial */}
      {polynomial && (
        <div style={{ background: "#ff6b6b11", border: "1px solid #ff6b6b33", borderRadius: 10, padding: 16 }}>
          <div style={{ fontSize: 12, color: "var(--text2)", marginBottom: 6 }}>Polinomio — {method}</div>
          <div style={{ fontFamily: "var(--mono)", fontSize: 13, color: "#ff9999", wordBreak: "break-all" }}>{polynomial}</div>
          {result && <div style={{ marginTop: 10, fontSize: 14, color: "#ff6b6b" }}>
            P({evalX}) = <strong>{result.yEval.toFixed(8)}</strong>
          </div>}
        </div>
      )}

      {/* Plot */}
      {plotData.length > 0 && (
        <div style={{ background: "var(--bg2)", borderRadius: 12, padding: 20, border: "1px solid var(--border)" }}>
          <h3 style={{ fontSize: 14, color: "var(--text2)", marginBottom: 16 }}>Gráfica de interpolación</h3>
          <ResponsiveContainer width="100%" height={280}>
            <LineChart data={plotData}>
              <CartesianGrid strokeDasharray="3 3" stroke="#2a2a38" />
              <XAxis dataKey="x" stroke="#5a5a78" tick={{ fontSize: 11 }} />
              <YAxis stroke="#5a5a78" tick={{ fontSize: 11 }} />
              <Tooltip contentStyle={{ background: "#111118", border: "1px solid #2a2a38", borderRadius: 8 }} />
              <Line type="monotone" dataKey="y" stroke="#ff6b6b" dot={false} strokeWidth={2} name="Polinomio" />
            </LineChart>
          </ResponsiveContainer>
          {/* Data points overlay */}
          <div style={{ marginTop: 12, display: "flex", gap: 8, flexWrap: "wrap" }}>
            {scatterData.map((p, i) => (
              <div key={i} style={{ background: "var(--bg3)", borderRadius: 6, padding: "4px 10px", fontSize: 12, fontFamily: "var(--mono)", color: "#ff9999" }}>
                ({p.x}, {p.y})
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Validation comparison */}
      {compErrors.length > 0 && (
        <div style={{ background: "var(--bg2)", borderRadius: 12, padding: 20, border: "1px solid var(--border)" }}>
          <h3 style={{ fontSize: 14, color: "var(--text2)", marginBottom: 16 }}>Error de validación por porcentaje (Vandermonde)</h3>
          <ResponsiveContainer width="100%" height={200}>
            <LineChart data={compErrors}>
              <CartesianGrid strokeDasharray="3 3" stroke="#2a2a38" />
              <XAxis dataKey="pct" stroke="#5a5a78" tick={{ fontSize: 11 }} tickFormatter={v => `${v}%`} />
              <YAxis stroke="#5a5a78" tick={{ fontSize: 11 }} />
              <Tooltip contentStyle={{ background: "#111118", border: "1px solid #2a2a38", borderRadius: 8 }} formatter={(v) => [v, "Error promedio"]} labelFormatter={v => `Validación: ${v}%`} />
              <Line type="monotone" dataKey="error" stroke="#ff6b6b" strokeWidth={2} dot={{ fill: "#ff6b6b", r: 5 }} />
            </LineChart>
          </ResponsiveContainer>
          <div style={{ display: "flex", gap: 12, marginTop: 12 }}>
            {compErrors.map(({ pct, error }) => (
              <div key={pct} style={{ background: "var(--bg3)", borderRadius: 8, padding: "8px 14px", textAlign: "center", flex: 1 }}>
                <div style={{ fontSize: 11, color: "var(--text2)" }}>{pct}% validación</div>
                <div style={{ fontFamily: "var(--mono)", fontSize: 14, color: error === null ? "var(--text3)" : "#ff9999", marginTop: 4 }}>
                  {error === null ? "N/A" : error}
                </div>
              </div>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
