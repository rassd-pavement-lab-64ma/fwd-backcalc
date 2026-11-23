"""Vercel serverless handler: POST /api/backcalc

Body: {
  "load_kN": 40, "plate_radius_mm": 150,
  "sensors_mm": [...], "deflections_um": [...],
  "layers": [{"thickness_mm": 135|null, "poisson": 0.35,
              "E_seed_MPa": 2700, "E_min_MPa": 690, "E_max_MPa": 13800}, ...]
}
Returns: moduli, RMS%, per-sensor fit.
"""
import json
from http.server import BaseHTTPRequestHandler

import numpy as np

from api.engine import backcalculate


class handler(BaseHTTPRequestHandler):
    def _send(self, code, obj):
        body = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        try:
            cfg = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
            layers = cfg["layers"]
            h_mm, poisson, E_seed, E_min, E_max = [], [], [], [], []
            for L in layers:
                th = L.get("thickness_mm")
                h_mm.append(float(th) if th not in (None, "") else 0.0)
                poisson.append(float(L["poisson"]))
                E_seed.append(float(L.get("E_seed_MPa") or
                                    0.5 * (float(L["E_min_MPa"]) + float(L["E_max_MPa"]))))
                E_min.append(float(L["E_min_MPa"]))
                E_max.append(float(L["E_max_MPa"]))
            h_mm = h_mm[:-1]
            q_kPa = float(cfg["load_kN"]) / (np.pi * (float(cfg["plate_radius_mm"]) / 1000.0) ** 2)
            res = backcalculate(
                [float(x) for x in cfg["deflections_um"]],
                [float(x) for x in cfg["sensors_mm"]],
                h_mm, poisson, E_seed, E_min, E_max,
                q_kPa, float(cfg["plate_radius_mm"]),
                fix=cfg.get("fix"), n_starts=int(cfg.get("n_starts", 3)))
            res["E_MPa"] = [float(x) for x in res["E_MPa"]]
            res["load_kPa"] = round(float(q_kPa), 1)
            self._send(200, res)
        except Exception as exc:  # noqa: BLE001
            self._send(400, {"error": str(exc)})

    def do_GET(self):
        self._send(200, {"status": "ok"})
