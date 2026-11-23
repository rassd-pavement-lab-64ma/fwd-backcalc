const $ = (id) => document.getElementById(id);

const DEFAULT_SENSORS = [0, 203, 305, 457, 610, 914, 1219, 1524];

function sensorRow(r = "", d = "") {
  const div = document.createElement("div");
  div.className = "row sensors";
  div.innerHTML = `<label>Sensor offset (mm)<input type="number" step="any" class="s-off" value="${r}"></label>
    <label>Deflection (µm)<input type="number" step="any" class="s-def" value="${d}"></label>`;
  return div;
}

function layerRow(L = {}) {
  const div = document.createElement("div");
  div.className = "row layers";
  div.innerHTML = `
    <label>Thickness (mm)<input type="number" step="any" class="l-th" value="${L.th ?? ""}" placeholder="∞"></label>
    <label>Poisson ν<input type="number" step="any" class="l-nu" value="${L.nu ?? 0.35}"></label>
    <label>E seed (MPa)<input type="number" step="any" class="l-seed" value="${L.seed ?? ""}"></label>
    <label>E min<input type="number" step="any" class="l-min" value="${L.min ?? ""}"></label>
    <label>E max<input type="number" step="any" class="l-max" value="${L.max ?? ""}"></label>
    <button class="btn-icon" title="Remove layer">×</button>`;
  div.querySelector(".btn-icon").onclick = () => {
    if (document.querySelectorAll("#layers .row.layers").length > 2) div.remove();
  };
  return div;
}

function loadExample() {
  $("load").value = 40; $("plate").value = 150;
  $("sensors").innerHTML = "";
  const offs = [0, 203, 305, 610, 914, 1219], defs = [421, 336, 288, 167, 106, 76];
  offs.forEach((o, i) => $("sensors").appendChild(sensorRow(o, defs[i])));
  $("layers").innerHTML = "";
  [
    { th: 152, nu: 0.35, seed: 2700, min: 690, max: 13800 },
    { th: 457, nu: 0.40, seed: 170, min: 35, max: 3500 },
    { th: null, nu: 0.45, seed: 100, min: 35, max: 3500 },
  ].forEach((L) => $("layers").appendChild(layerRow(L)));
}

function collect() {
  const sensors = [...document.querySelectorAll("#sensors .s-off")].map((e) => +e.value);
  const deflections = [...document.querySelectorAll("#sensors .s-def")].map((e) => +e.value);
  const layers = [...document.querySelectorAll("#layers .row.layers")].map((r) => ({
    thickness_mm: r.querySelector(".l-th").value === "" ? null : +r.querySelector(".l-th").value,
    poisson: +r.querySelector(".l-nu").value,
    E_seed_MPa: r.querySelector(".l-seed").value === "" ? null : +r.querySelector(".l-seed").value,
    E_min_MPa: +r.querySelector(".l-min").value,
    E_max_MPa: +r.querySelector(".l-max").value,
  }));
  return { load_kN: +$("load").value, plate_radius_mm: +$("plate").value,
           sensors_mm: sensors, deflections_um: deflections, layers, n_starts: 3 };
}

function drawChart(meas, calc, sensors) {
  const cv = $("chart"), ctx = cv.getContext("2d");
  const W = (cv.width = cv.clientWidth * 2), H = (cv.height = 440);
  ctx.clearRect(0, 0, W, H);
  const pad = 46;
  const xmax = Math.max(...sensors) * 1.05;
  const ymax = Math.max(...meas, ...calc) * 1.15;
  const x = (v) => pad + (v / xmax) * (W - pad - 12);
  const y = (v) => H - pad - (v / ymax) * (H - pad - 16);
  ctx.strokeStyle = "#d9dde3"; ctx.fillStyle = "#6b7280"; ctx.font = "20px sans-serif";
  for (let i = 0; i <= 5; i++) {
    const gy = pad + (i / 5) * (H - pad - 16);
    ctx.beginPath(); ctx.moveTo(pad, gy); ctx.lineTo(W - 12, gy); ctx.stroke();
    ctx.fillText(Math.round(ymax * (1 - i / 5)), 4, gy + 7);
  }
  const plot = (arr, color) => {
    ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = 2.5;
    ctx.beginPath();
    arr.forEach((v, i) => { i ? ctx.lineTo(x(sensors[i]), y(v)) : ctx.moveTo(x(sensors[i]), y(v)); });
    ctx.stroke();
    arr.forEach((v, i) => { ctx.beginPath(); ctx.arc(x(sensors[i]), y(v), 4.5, 0, 7); ctx.fill(); });
  };
  plot(meas, "#1a4f8b"); plot(calc, "#b45309");
  ctx.fillStyle = "#1a4f8b"; ctx.fillText("■ Measured", pad, H - 14);
  ctx.fillStyle = "#b45309"; ctx.fillText("■ Calculated", pad + 140, H - 14);
  ctx.fillStyle = "#6b7280"; ctx.fillText("Deflection (µm) vs sensor offset (mm)", pad, H - 34);
}

async function run() {
  const status = $("run-status");
  status.textContent = "Running…";
  $("btn-run").disabled = true;
  try {
    const res = await fetch("/api/backcalc", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify(collect()),
    });
    const data = await res.json();
    if (data.error) throw new Error(data.error);
    const badge = $("rms-badge");
    badge.classList.remove("hidden", "warn");
    badge.textContent = `RMS ${data.rms_percent}%`;
    if (data.rms_percent > 3) badge.classList.add("warn");
    const names = ["Layer 1", "Layer 2", "Layer 3", "Layer 4", "Layer 5"];
    $("moduli-cards").innerHTML = data.E_MPa.map((e, i) =>
      `<div class="card"><div class="k">${names[i] || "Layer " + (i + 1)}</div>
       <div class="v">${Math.round(e).toLocaleString()} <span class="u">MPa</span></div></div>`).join("");
    const tb = document.querySelector("#fit-table tbody");
    tb.innerHTML = data.meas_defl_um.map((m, i) =>
      `<tr><td>${collect().sensors_mm[i]}</td><td>${m}</td><td>${data.calc_defl_um[i]}</td>
       <td>${data.rel_error_percent[i]}</td></tr>`).join("");
    $("fit-table").classList.remove("hidden");
    drawChart(data.meas_defl_um, data.calc_defl_um, collect().sensors_mm);
    status.textContent = `Done in ${data.n_function_evals} forward evaluations · load ${data.load_kPa} kPa`;
  } catch (e) {
    status.textContent = "Error: " + e.message;
  } finally {
    $("btn-run").disabled = false;
  }
}

$("btn-example").onclick = loadExample;
$("btn-addlayer").onclick = () => $("layers").appendChild(layerRow({}));
$("btn-run").onclick = run;
loadExample();
