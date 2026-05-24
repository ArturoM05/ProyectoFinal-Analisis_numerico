import { useState } from "react";
import Chapter1 from "./components/Chapter1";
import Chapter2 from "./components/Chapter2";
import Chapter3 from "./components/Chapter3";

const chapters = [
  { id: 1, label: "Cap. 1", title: "Raíces de Ecuaciones", icon: "√x", color: "#7c6af7" },
  { id: 2, label: "Cap. 2", title: "Sistemas Lineales", icon: "Ax=b", color: "#4ecdc4" },
  { id: 3, label: "Cap. 3", title: "Interpolación", icon: "P(x)", color: "#ff6b6b" },
];

export default function App() {
  const [active, setActive] = useState(1);

  return (
    <div style={{ minHeight: "100vh", display: "flex", flexDirection: "column" }}>
      {/* Header */}
      <header style={{
        background: "var(--bg2)", borderBottom: "1px solid var(--border)",
        padding: "0 2rem", display: "flex", alignItems: "center",
        gap: "2rem", height: 64, position: "sticky", top: 0, zIndex: 100
      }}>
        <div style={{ display: "flex", alignItems: "center", gap: 12 }}>
          <div style={{
            width: 36, height: 36, borderRadius: 10,
            background: "linear-gradient(135deg, #7c6af7, #4ecdc4)",
            display: "flex", alignItems: "center", justifyContent: "center",
            fontSize: 16, fontWeight: 700, fontFamily: "var(--mono)"
          }}>∑</div>
          <span style={{ fontWeight: 600, fontSize: 16, letterSpacing: "-0.3px" }}>
            Métodos Numéricos
          </span>
        </div>
        <nav style={{ display: "flex", gap: 4, flex: 1 }}>
          {chapters.map(ch => (
            <button key={ch.id} onClick={() => setActive(ch.id)}
              style={{
                background: active === ch.id ? ch.color + "22" : "transparent",
                color: active === ch.id ? ch.color : "var(--text2)",
                border: active === ch.id ? `1px solid ${ch.color}44` : "1px solid transparent",
                borderRadius: 8, padding: "6px 16px", cursor: "pointer",
                fontSize: 13, fontWeight: 500, transition: "all 0.2s",
                fontFamily: "var(--sans)"
              }}>
              <span style={{ fontFamily: "var(--mono)", marginRight: 6 }}>{ch.label}</span>
              {ch.title}
            </button>
          ))}
        </nav>
      </header>

      {/* Content */}
      <main style={{ flex: 1, padding: "2rem", maxWidth: 1200, margin: "0 auto", width: "100%" }}>
        {active === 1 && <Chapter1 />}
        {active === 2 && <Chapter2 />}
        {active === 3 && <Chapter3 />}
      </main>
    </div>
  );
}
