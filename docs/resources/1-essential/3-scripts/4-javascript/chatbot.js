// static/script.js

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// LÓGICA DEL CHATBOT CON EMBEDDINGS Y TRANSFORMERS.JS
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//=============================================================================================================
// IMPORTS Y VARIABLES GLOBALES
//=============================================================================================================

import { embeddingsData } from "./embeddings_web_reduced.js";
import { pipeline } from "https://cdn.jsdelivr.net/npm/@xenova/transformers/dist/transformers.min.js";
import predefResponses from "./responses.js";

let embedder = null;





//=============================================================================================================
// CONFIGURACIÓN
//=============================================================================================================

const QA_PRIOR_DISPLAY = 2;
const QA_MAX = 3;
const QMD_MAX = 3;
const PDF_MAX = 2;
const REDUNDANCY_THRESHOLD = 0.88;
const MAX_SUMMARY_SENTENCES = 10;
const BOT_TYPING_DELAY = 600;

const seenFragmentKeys = new Set();
const _embeddingCache = {};





//=============================================================================================================
// FUNCIONES AUXILIARES
//=============================================================================================================

function cosineSim(a, b) {
  let dot = 0,
    normA = 0,
    normB = 0;
  for (let i = 0; i < a.length; i++) {
    dot += a[i] * b[i];
    normA += a[i] ** 2;
    normB += b[i] ** 2;
  }
  return dot / (Math.sqrt(normA) * Math.sqrt(normB));
}

function cleanFragment(text) {
  return text
    .replace(/^(Resumen|Análisis|Revisión del archivo de metadatos)[:\-]*/i, "")
    .replace(/\s{2,}/g, " ")
    .trim();
}

function extractAnswerOnly(text) {
  const m = text.match(/^¿.*?\?\s*(.*)/s);
  return m ? m[1].trim() : text.trim();
}

// Normalización básica (acentos, signos, espacios)
function normalize(str) {
  return str
    .toLowerCase()
    .normalize("NFD")
    .replace(/[\u0300-\u036f]/g, "")
    .replace(/[¿?.,!]/g, "")
    .trim();
}





//=============================================================================================================
// EMBEDDINGS
//=============================================================================================================

let embedderReady = setupEmbedder();

async function getEmbedding(text) {
  await embedderReady;
  if (_embeddingCache[text]) return _embeddingCache[text];
  const out = await embedder(text, { pooling: "mean", normalize: true });
  _embeddingCache[text] = out.data;
  return out.data;
}




//=============================================================================================================
// SEPARAR QA / OTROS
//=============================================================================================================

embeddingsData.qa = embeddingsData.filter((d) =>
  d.file.toLowerCase().startsWith("qa/")
);
embeddingsData.qmd = embeddingsData.filter((d) =>
  d.file.toLowerCase().startsWith("qmd/")
);
embeddingsData.pdf = embeddingsData.filter((d) =>
  d.file.toLowerCase().startsWith("pdf/")
);





//=============================================================================================================
// SETUP TRANSFORMERS
//=============================================================================================================

async function setupEmbedder() {
  embedder = await pipeline("feature-extraction", "Xenova/all-mpnet-base-v2");
  console.log(
    "Transformers.js listo con all-mpnet-base-v2 y embeddings cargados"
  );
}





//=============================================================================================================
// EXTRAER PALABRAS CLAVE
//=============================================================================================================

function extractKeywords(query) {
  return query
    .toLowerCase()
    .replace(/[¿?.,!]/g, "")
    .replace(/^(qué|quien|cómo|cuando|donde|porque|puedes|dime)\s+/i, "")
    .split(/\s+/)
    .filter((w) => w.length > 2);
}





//=============================================================================================================
// BUSCAR FRAGMENTOS Y RE-RANK
//=============================================================================================================

async function searchEmbeddings(query, data, topN = 5) {
  const qNorm = normalize(query);
  const qEmb = await getEmbedding(query);
  const keywords = extractKeywords(query);

  let fragments = data.map((d) => ({
    text: cleanFragment(d.text),
    embedding: d.embedding || null,
    file: d.file,
  }));

  // 1. PRIORIDAD: coincidencia literal exacta
  for (let f of fragments) {
    const fNorm = normalize(f.text);
    if (fNorm.includes(qNorm)) {
      return [f, ...fragments.filter((x) => x !== f).slice(0, topN - 1)];
    }
  }

  // 2. Si no hay match literal, usar similitud coseno + keywords
  let scored = [];
  for (let f of fragments) {
    const fEmb = f.embedding || (await getEmbedding(f.text));
    let simScore = cosineSim(qEmb, fEmb);

    // Coincidencias literales de keywords
    const lowerText = f.text.toLowerCase();
    let keywordHits = 0;
    for (let kw of keywords) if (lowerText.includes(kw)) keywordHits++;

    // Refuerzo adicional
    let keywordBoost = 1 + keywordHits / (keywords.length || 1);
    if (f.file.toLowerCase().startsWith("qa/") && keywordHits > 0) {
      simScore *= 50; // boost fuerte para QA relevantes
    } else {
      simScore *= keywordBoost;
    }

    if (keywordHits === 0) simScore *= 0.85;

    scored.push({ ...f, score: simScore, keywordHits });
  }

  // Orden final
  scored.sort((a, b) => {
    if (b.keywordHits !== a.keywordHits) return b.keywordHits - a.keywordHits;
    return b.score - a.score;
  });

  return scored.slice(0, topN);
}





//=============================================================================================================
// FINALIZAR RESUMEN
//=============================================================================================================

function finalizeSummary(summary) {
  if (!summary) return "";
  summary = summary
    .replace(/\s+\./g, ".")
    .replace(/\s+,/g, ",")
    .replace(/\s+;/g, ";")
    .replace(/\s+:/g, ":")
    .replace(/\s+\?/g, "?")
    .replace(/\s+!/g, "!")
    .replace("▼", "")
    .replace("→", "")
    .replace(":", "")
    .replace(/\bademás,?\b\s*(además,?\s*)+/gi, "Además, ")
    .replace(/\ben definitiva,?\b\s*(en definitiva,?\s*)+/gi, "En definitiva, ")
    .replace(/^\s*Análisis\s*:?/i, "")
    .replace(/^\s*[A-Za-z]+ \d+(\.\d+)*\s*:?\s*/, "")
    .trim();
  return summary[0].toUpperCase() + summary.slice(1);
}





//=============================================================================================================
// GENERAR RESUMEN COMPLETO
//=============================================================================================================

async function generateSummary(query, combined) {
  const selected = [];

  const addUnique = (items, limit, isQA = false) => {
    let count = 0;
    for (let f of items) {
      if (count >= limit) break;
      const key = `${f.file}::${f.text.slice(0, 300)}`;
      if (!seenFragmentKeys.has(key)) {
        seenFragmentKeys.add(key);
        selected.push({ ...f, isQA });
        count++;
      }
    }
  };

  addUnique(
    combined.filter((f) => f.file.toLowerCase().startsWith("qa/")),
    QA_PRIOR_DISPLAY,
    true
  );
  addUnique(
    combined.filter(
      (f) => f.file.toLowerCase().startsWith("qa/") && !selected.includes(f)
    ),
    QA_MAX - QA_PRIOR_DISPLAY,
    true
  );
  addUnique(
    combined.filter((f) => f.file.toLowerCase().startsWith("qmd/")),
    QMD_MAX
  );
  addUnique(
    combined.filter((f) => f.file.toLowerCase().startsWith("pdf/")),
    PDF_MAX
  );

  const summaryText = selected
    .map((f) => (f.isQA ? extractAnswerOnly(f.text) : f.text))
    .join(" ");

  return finalizeSummary(summaryText);
}





//=============================================================================================================
// RESPUESTAS PREDEFINIDAS
//=============================================================================================================

function checkPredefinedResponse(query) {
  // 🔧 Normalización completa
  const q = query
    .toLowerCase()
    .normalize("NFD") // separa letras y tildes
    .replace(/[\u0300-\u036f]/g, "") // elimina tildes/acentos
    .replace(/[¿?¡!.,;:()"'«»]/g, "") // elimina signos de puntuación
    .replace(/\s+/g, " ") // colapsa múltiples espacios
    .trim(); // elimina espacios al inicio y final

  // Buscar coincidencias
  console.log("🔎 Buscando coincidencia con:", q);
  for (let item of predefResponses) {
    const regex = new RegExp(item.pattern, "i");
    if (regex.test(q)) {
      console.log("✅ Coincidencia encontrada con patrón:", item.pattern);
      return item.response;
    }
  }
  console.log("❌ Ninguna coincidencia encontrada");
  return null;
}





//=============================================================================================================
// INTERFAZ CON DOM
//=============================================================================================================

async function handleQuery(query) {
  if (!query.trim()) return;

  const chatInfo = document.getElementById("chat-info");
  if (chatInfo && !chatInfo.classList.contains("hidden")) {
    chatInfo.style.opacity = 0;
    setTimeout(() => {
      chatInfo.classList.add("hidden");
    }, 300);
  }

  displayMessage("user", query);

  const predef = checkPredefinedResponse(query);
  if (predef) {
    displayMessage("bot", predef);
    return;
  }

  const typingBubble = displayMessage(
    "bot",
    `<span class="typing-spinner"></span><span class="typing-text">Generando respuesta...</span>`
  );

  const qaCandidates = await searchEmbeddings(query, embeddingsData.qa, QA_MAX);
  const qmdCandidates = await searchEmbeddings(
    query,
    embeddingsData.qmd,
    QMD_MAX
  );
  const pdfCandidates = await searchEmbeddings(
    query,
    embeddingsData.pdf,
    PDF_MAX
  );

  const combined = [
    ...qaCandidates.map((q) => ({ ...q, isQA: true })),
    ...qmdCandidates,
    ...pdfCandidates,
  ];

  const summary = await generateSummary(query, combined);

  const highlighted = highlightKeywords(summary, query);

  setTimeout(() => {
    typingBubble.innerHTML = highlighted || "No se encontró respuesta.";
    typingBubble.classList.add("bot-final");
    const log = document.getElementById("chat-log");
    log.scrollTop = log.scrollHeight;
  }, BOT_TYPING_DELAY);
}





//=============================================================================================================
// Resaltar palabras clave
//=============================================================================================================

function highlightKeywords(text, query) {
  if (!text) return text;
  const words = extractKeywords(query).map((w) =>
    w.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")
  );
  if (!words.length) return text;
  const regex = new RegExp(`\\b(${words.join("|")})\\b`, "gi");
  return text.replace(regex, "<mark>$1</mark>");
}





//=============================================================================================================
// Mostrar mensajes en el chat
//=============================================================================================================

function displayMessage(sender, text) {
  const log = document.getElementById("chat-log");
  const msg = document.createElement("div");
  msg.className = "message " + sender;
  msg.innerHTML = text;
  log.appendChild(msg);
  log.scrollTop = log.scrollHeight;
  return msg;
}





//=============================================================================================================
// EVENTOS DOM
//=============================================================================================================

document
  .getElementById("chat-input")
  .addEventListener("keypress", async (e) => {
    if (e.key === "Enter") {
      await embedderReady;
      const query = e.target.value.trim();
      if (!query) return;
      e.target.value = "";
      await handleQuery(query);
    }
  });

document.getElementById("chat-toggle").addEventListener("click", () => {
  document.getElementById("chat-box").classList.toggle("active");
});

const chatInput = document.getElementById("chat-input");
const chatSend = document.getElementById("chat-send");

chatSend.addEventListener("click", async () => {
  const query = chatInput.value.trim();
  if (!query) return;
  chatInput.value = "";
  await handleQuery(query);
});

chatInput.addEventListener("keypress", async (e) => {
  if (e.key === "Enter") {
    const query = chatInput.value.trim();
    if (!query) return;
    chatInput.value = "";
    await handleQuery(query);
  }
});
