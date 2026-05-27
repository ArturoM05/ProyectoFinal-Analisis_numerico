import { useState } from "react";
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from "recharts";
import * as math from "mathjs";

const METHODS = ["Bisección", "Regla Falsa", "Punto Fijo", "Newton-Raphson", "Secante", "Raíces Múltiples"];
const COLORS = ["#7c6af7", "#4ecdc4", "#ff6b6b", "#ffd93d", "#a8ff78", "#ff9a9e"];

function safeEval(expr, x) {
  try { return math.evaluate(expr, { x }); }
  catch { return NaN; }
}

// ---- Calcula los 3 tipos de error para un par (xNew, xOld) ----
function calcAllErrors(xNew, xOld, f) {
  if (xOld === null || xOld === undefined) return { Relativo: null, Absoluto: null, Condicion: null };

  const absoluto = Math.abs(xNew - xOld);
  const relativo = Math.abs((xNew - xOld) / (xNew !== 0 ? xNew : 1));

  let condicion = null;
  if (f) {
    const fx = safeEval(f, xNew);
    const dfx = (safeEval(f, xNew + 1e-7) - safeEval(f, xNew - 1e-7)) / (2e-7);
    if (Math.abs(fx) > 1e-14 && Math.abs(dfx) > 1e-14) {
      condicion = Math.abs(xNew * dfx / fx);
    }
  }

  return {
    Relativo: isFinite(relativo) ? relativo : null,
    Absoluto: isFinite(absoluto) ? absoluto : null,
    Condicion: condicion !== null && isFinite(condicion) ? condicion : null,
  };
}

// ---- Error simple para criterio de parada ----
function calcError(xNew, xOld, f, errType) {
  if (xOld === null) return Infinity;
  if (errType === "Relativo") return Math.abs((xNew - xOld) / (xNew || 1));
  if (errType === "Absoluto") return Math.abs(xNew - xOld);
  if (errType === "Condición" && f) {
    const fx = safeEval(f, xNew);
    if (Math.abs(fx) < 1e-14) return Infinity;
    const dfx = (safeEval(f, xNew + 1e-7) - safeEval(f, xNew - 1e-7)) / (2e-7);
    return Math.abs(xNew * dfx / fx);
  }
  return Math.abs(xNew - xOld);
}

function bisection(f, a, b, tol, maxIter, errType) {
  const rows = [];
  let xa = a, xb = b, xr, xrPrev = null;
  for (let i = 1; i <= maxIter; i++) {
    xr = (xa + xb) / 2;
    const fa = safeEval(f, xa), fr = safeEval(f, xr);
    const err = calcError(xr, xrPrev, f, errType);
    const allErr = calcAllErrors(xr, xrPrev, f);
    rows.push({ iter: i, xa: +xa.toFixed(6), xb: +xb.toFixed(6), xr: +xr.toFixed(6), fr: +fr.toFixed(6), error: isFinite(err) ? +err.toFixed(6) : "---", ...allErr });
    if (xrPrev !== null && isFinite(err) && err < tol) break;
    if (fa * fr < 0) xb = xr; else xa = xr;
    xrPrev = xr;
    if (Math.abs(fr) < 1e-12) break;
  }
  return { rows, root: xr };
}

function falsePosition(f, a, b, tol, maxIter, errType) {
  const rows = [];
  let xa = a, xb = b, xr, xrPrev = null;
  for (let i = 1; i <= maxIter; i++) {
    const fa = safeEval(f, xa), fb = safeEval(f, xb);
    xr = xb - fb * (xa - xb) / (fa - fb);
    const fr = safeEval(f, xr);
    const err = calcError(xr, xrPrev, f, errType);
    const allErr = calcAllErrors(xr, xrPrev, f);
    rows.push({ iter: i, xa: +xa.toFixed(6), xb: +xb.toFixed(6), xr: +xr.toFixed(6), fr: +fr.toFixed(6), error: isFinite(err) ? +err.toFixed(6) : "---", ...allErr });
    if (xrPrev !== null && isFinite(err) && err < tol) break;
    if (fa * fr < 0) xb = xr; else xa = xr;
    xrPrev = xr;
    if (Math.abs(fr) < 1e-12) break;
  }
  return { rows, root: xr };
}

function fixedPoint(g, x0, tol, maxIter, errType) {
  const rows = [];
  let x = x0, err = tol + 1;
  for (let i = 0; i <= maxIter && err > tol; i++) {
    const fx = safeEval(g, x) - x;
    const xNew = safeEval(g, x);
    err = Math.abs(xNew - x);
    const allErr = calcAllErrors(xNew, x, null);
    rows.push({ iter: i, x: +x.toPrecision(15), fx: +fx.toPrecision(15), error: isFinite(err) ? +err.toPrecision(15) : "---", ...allErr });
    x = xNew;
    if (!isFinite(x)) break;
  }
  return { rows, root: x };
}

function newtonRaphson(f, x0, tol, maxIter, errType) {
  const rows = [];
  let x = x0, err = Infinity;
  for (let i = 1; i <= maxIter && err > tol; i++) {
    const fx = safeEval(f, x);
    const dfx = (safeEval(f, x + 1e-7) - safeEval(f, x - 1e-7)) / (2e-7);
    if (Math.abs(dfx) < 1e-14) break;
    const xNew = x - fx / dfx;
    err = calcError(xNew, x, f, errType);
    const allErr = calcAllErrors(xNew, x, f);
    rows.push({ iter: i, x: +x.toFixed(6), fx: +fx.toFixed(6), dfx: +dfx.toFixed(6), xNew: +xNew.toFixed(6), error: isFinite(err) ? +err.toFixed(6) : "---", ...allErr });
    x = xNew;
    if (Math.abs(fx) < 1e-12) break;
  }
  return { rows, root: x };
}

function secant(f, x0, x1, tol, maxIter, errType) {
  const rows = [];
  let xa = x0, xb = x1, err = Infinity;
  for (let i = 1; i <= maxIter && err > tol; i++) {
    const fa = safeEval(f, xa), fb = safeEval(f, xb);
    if (Math.abs(fb - fa) < 1e-14) break;
    const xNew = xb - fb * (xb - xa) / (fb - fa);
    err = calcError(xNew, xb, f, errType);
    const allErr = calcAllErrors(xNew, xb, f);
    rows.push({ iter: i, xa: +xa.toFixed(6), xb: +xb.toFixed(6), fb: +fb.toFixed(6), xNew: +xNew.toFixed(6), error: isFinite(err) ? +err.toFixed(6) : "---", ...allErr });
    xa = xb; xb = xNew;
    if (Math.abs(fb) < 1e-12) break;
  }
  return { rows, root: xb };
}

function multipleRoots(f, x0, tol, maxIter, errType) {
  const rows = [];
  let x = x0, err = Infinity;
  for (let i = 1; i <= maxIter && err > tol; i++) {
    const fx = safeEval(f, x);
    const dfx = (safeEval(f, x + 1e-7) - safeEval(f, x - 1e-7)) / (2e-7);
    const d2fx = (safeEval(f, x + 1e-7) - 2 * fx + safeEval(f, x - 1e-7)) / (1e-14);
    const denom = dfx * dfx - fx * d2fx;
    if (Math.abs(denom) < 1e-14) break;
    const xNew = x - fx * dfx / denom;
    err = calcError(xNew, x, f, errType);
    const allErr = calcAllErrors(xNew, x, f);
    rows.push({ iter: i, x: +x.toFixed(6), fx: +fx.toFixed(6), xNew: +xNew.toFixed(6), error: isFinite(err) ? +err.toFixed(6) : "---", ...allErr });
    x = xNew;
    if (Math.abs(fx) < 1e-12) break;
  }
  return { rows, root: x };
}

function plotFunction(expr, a, b) {
  const pts = [];
  for (let i = 0; i <= 200; i++) {
    const x = a + (b - a) * i / 200;
    const y = safeEval(expr, x);
    if (isFinite(y)) pts.push({ x: +x.toFixed(4), y: +y.toFixed(4) });
  }
  return pts;
}

const HELP = {
  "Bisección": "Ingresa f(x), intervalo [a,b] donde f(a)·f(b)<0.",
  "Regla Falsa": "Ingresa f(x), intervalo [a,b] donde f(a)·f(b)<0.",
  "Punto Fijo": "Despeja x = g(x) de f(x)=0. Ingresa g(x) y x₀.",
  "Newton-Raphson": "Ingresa f(x) y punto inicial x₀. La derivada se calcula automáticamente.",
  "Secante": "Ingresa f(x) y dos puntos iniciales x₀, x₁.",
  "Raíces Múltiples": "Ingresa f(x) y x₀. Útil cuando la raíz tiene multiplicidad > 1.",
};

export default function Chapter1() {
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
  const [error, setError] = useState("");

  const run = () => {
    setError("");
    try {
      const t = parseFloat(tol), m = parseInt(maxIter);
      let res;
      if (method === "Bisección") res = bisection(fx, parseFloat(a), parseFloat(b), t, m, errType);
      else if (method === "Regla Falsa") res = falsePosition(fx, parseFloat(a), parseFloat(b), t, m, errType);
      else if (method === "Punto Fijo") res = fixedPoint(gx, parseFloat(x0), t, m, errType);
      else if (method === "Newton-Raphson") res = newtonRaphson(fx, parseFloat(x0), t, m, errType);
      else if (method === "Secante") res = secant(fx, parseFloat(x0), parseFloat(x1), t, m, errType);
      else res = multipleRoots(fx, parseFloat(x0), t, m, errType);
      setResult(res);
      setCompData([]); // limpiar comparativa al ejecutar
      const pa = method === "Punto Fijo" ? parseFloat(x0) - 2 : parseFloat(a);
      const pb = method === "Punto Fijo" ? parseFloat(x0) + 2 : parseFloat(b);
      setPlotData(plotFunction(method === "Punto Fijo" ? gx : fx, pa - 1, pb + 1));
    } catch (e) { setError("Error: " + e.message); }
  };

  // Comparativa: 3 líneas (Relativo, Absoluto, Condición) por iteración del método actual
  const runComparison = () => {
    setError("");
    try {
      const t = parseFloat(tol), m = parseInt(maxIter);
      let res;
      if (method === "Bisección") res = bisection(fx, parseFloat(a), parseFloat(b), t, m, errType);
      else if (method === "Regla Falsa") res = falsePosition(fx, parseFloat(a), parseFloat(b), t, m, errType);
      else if (method === "Punto Fijo") res = fixedPoint(gx, parseFloat(x0), t, m, errType);
      else if (method === "Newton-Raphson") res = newtonRaphson(fx, parseFloat(x0), t, m, errType);
      else if (method === "Secante") res = secant(fx, parseFloat(x0), parseFloat(x1), t, m, errType);
      else res = multipleRoots(fx, parseFloat(x0), t, m, errType);

      // Construir datos para la gráfica: una entrada por iteración con los 3 errores
      const data = res.rows
        .filter(row => row.iter > 1) // saltar iter 1 donde no hay error previo
        .map(row => ({
          iter: row.iter,
          Relativo:  row.Relativo  !== null && row.Relativo  > 0 ? row.Relativo  : null,
          Absoluto:  row.Absoluto  !== null && row.Absoluto  > 0 ? row.Absoluto  : null,
          Condicion: row.Condicion !== null && row.Condicion > 0 ? row.Condicion : null,
        }));

      setCompData(data);
    } catch (e) { setError("Error en comparación: " + e.message); }
  };

  const cols = {
    "Bisección": ["iter", "xa", "xb", "xr", "fr", "error"],
    "Regla Falsa": ["iter", "xa", "xb", "xr", "fr", "error"],
    "Punto Fijo": ["iter", "x", "fx", "error"],
    "Newton-Raphson": ["iter", "x", "fx", "dfx", "xNew", "error"],
    "Secante": ["iter", "xa", "xb", "fb", "xNew", "error"],
    "Raíces Múltiples": ["iter", "x", "fx", "xNew", "error"],
  };

  const colLabels = {
    iter: "i", xa: "a", xb: "b", xr: "xr", fr: "f(r)",
    x: "x_i", fx: "F(x_i)", xNew: "x_{i+1}", dfx: "f'(x)", fb: "f(b)", error: "E_i",
  };

  return (
    <div style={{ display: "flex", flexDirection: "column", gap: 24 }}>
      <div>
        <h1 style={{ fontSize: 24, fontWeight: 600, color: "#7c6af7" }}>Raíces de Ecuaciones</h1>
        <p style={{ color: "var(--text2)", fontSize: 14, marginTop: 4 }}>Encuentra raíces de f(x) = 0 con diferentes métodos iterativos</p>
      </div>

      {/* Method selector */}
      <div style={{ display: "flex", gap: 8, flexWrap: "wrap" }}>
        {METHODS.map(m => (
          <button key={m} onClick={() => { setMethod(m); setResult(null); setCompData([]); }} style={{
            background: method === m ? "#7c6af722" : "var(--bg3)",
            color: method === m ? "#7c6af7" : "var(--text2)",
            border: method === m ? "1px solid #7c6af744" : "1px solid var(--border)",
            borderRadius: 8, padding: "6px 14px", cursor: "pointer",
            fontSize: 13, fontFamily: "var(--sans)"
          }}>{m}</button>
        ))}
      </div>

      {/* Help box */}
      <div style={{ background: "#7c6af711", border: "1px solid #7c6af733", borderRadius: 10, padding: "12px 16px", fontSize: 13, color: "#a89cf7" }}>
        💡 <strong>Ayuda:</strong> {HELP[method]} Usa sintaxis: <code>x^2</code>, <code>sin(x)</code>, <code>exp(x)</code>, <code>log(x)</code>, <code>sqrt(x)</code>
      </div>

      {/* Inputs */}
      <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(200px, 1fr))", gap: 16, background: "var(--bg2)", borderRadius: 12, padding: 20, border: "1px solid var(--border)" }}>
        {method !== "Punto Fijo" && (
          <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
            <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>f(x)</span>
            <input value={fx} onChange={e => setFx(e.target.value)} placeholder="x^3 - x - 2" />
          </label>
        )}
        {method === "Punto Fijo" && (
          <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
            <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>g(x)</span>
            <input value={gx} onChange={e => setGx(e.target.value)} placeholder="(x + 2)^(1/3)" />
          </label>
        )}
        {["Bisección", "Regla Falsa"].includes(method) && (<>
          <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
            <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>a (límite inferior)</span>
            <input value={a} onChange={e => setA(e.target.value)} />
          </label>
          <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
            <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>b (límite superior)</span>
            <input value={b} onChange={e => setB(e.target.value)} />
          </label>
        </>)}
        {["Punto Fijo", "Newton-Raphson", "Raíces Múltiples"].includes(method) && (
          <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
            <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>x₀ (punto inicial)</span>
            <input value={x0} onChange={e => setX0(e.target.value)} />
          </label>
        )}
        {method === "Secante" && (<>
          <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
            <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>x₀</span>
            <input value={x0} onChange={e => setX0(e.target.value)} />
          </label>
          <label style={{ display: "flex", flexDirection: "column", gap: 6 }}>
            <span style={{ fontSize: 12, color: "var(--text2)", fontFamily: "var(--mono)" }}>x₁</span>
            <input value={x1} onChange={e => setX1(e.target.value)} />
          </label>
        </>)}
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
            {["Relativo", "Absoluto", "Condición"].map(e => <option key={e}>{e}</option>)}
          </select>
        </label>
        <div style={{ display: "flex", alignItems: "flex-end", gap: 8 }}>
          <button onClick={run} style={{ flex: 1 }}>▶ Ejecutar</button>
          <button onClick={runComparison} className="secondary" style={{ flex: 1 }}>⇄ Comparar errores</button>
        </div>
      </div>

      {error && <div style={{ background: "#ff6b6b22", border: "1px solid #ff6b6b44", borderRadius: 8, padding: 12, color: "#ff6b6b", fontSize: 13 }}>{error}</div>}

      {/* Plot */}
      {plotData.length > 0 && (
        <div style={{ background: "var(--bg2)", borderRadius: 12, padding: 20, border: "1px solid var(--border)" }}>
          <h3 style={{ fontSize: 14, color: "var(--text2)", marginBottom: 16 }}>Gráfica de la función</h3>
          <ResponsiveContainer width="100%" height={250}>
            <LineChart data={plotData}>
              <CartesianGrid strokeDasharray="3 3" stroke="#2a2a38" />
              <XAxis dataKey="x" stroke="#5a5a78" tick={{ fontSize: 11 }} />
              <YAxis stroke="#5a5a78" tick={{ fontSize: 11 }} />
              <Tooltip contentStyle={{ background: "#111118", border: "1px solid #2a2a38", borderRadius: 8 }} />
              <Line type="monotone" dataKey="y" stroke="#7c6af7" dot={false} strokeWidth={2} />
            </LineChart>
          </ResponsiveContainer>
          {result && (
            <div style={{ marginTop: 12, padding: 12, background: "#7c6af711", borderRadius: 8, fontFamily: "var(--mono)", fontSize: 13, color: "#a89cf7" }}>
              ✓ Raíz aproximada: <strong>{result.root?.toFixed(8)}</strong> — encontrada en <strong>{result.rows.length}</strong> iteraciones
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
              <thead><tr>{cols[method].map(c => <th key={c}>{colLabels[c] || c}</th>)}</tr></thead>
              <tbody>{result.rows.map((row, i) => (
                <tr key={i}>{cols[method].map(c => <td key={c}>{row[c] ?? "---"}</td>)}</tr>
              ))}</tbody>
            </table>
          </div>
        </div>
      )}

      {/* Comparison — 3 líneas de error por iteración */}
      {compData.length > 0 && (
        <div style={{ background: "var(--bg2)", borderRadius: 12, padding: 20, border: "1px solid var(--border)" }}>
          <h3 style={{ fontSize: 14, color: "var(--text2)", marginBottom: 4 }}>
            Comparación de tipos de error — <span style={{ color: "#7c6af7" }}>{method}</span>
          </h3>
          <p style={{ fontSize: 12, color: "var(--text3)", marginBottom: 16 }}>
            Evolución del error relativo, absoluto y de condición a través de las iteraciones (escala logarítmica)
          </p>
          <ResponsiveContainer width="100%" height={300}>
            <LineChart data={compData} margin={{ top: 5, right: 20, bottom: 5, left: 10 }}>
              <CartesianGrid strokeDasharray="3 3" stroke="#2a2a38" />
              <XAxis dataKey="iter" stroke="#5a5a78" tick={{ fontSize: 11 }} label={{ value: "Iteración", position: "insideBottom", offset: -2, fill: "#5a5a78", fontSize: 11 }} />
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

          {/* Resumen final */}
          <div style={{ display: "flex", gap: 12, marginTop: 16 }}>
            {[
              { label: "Error Relativo final",  key: "Relativo",  color: "#7c6af7" },
              { label: "Error Absoluto final",  key: "Absoluto",  color: "#4ecdc4" },
              { label: "Error Condición final", key: "Condicion", color: "#ff6b6b" },
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