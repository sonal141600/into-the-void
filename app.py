
from flask import Flask, jsonify, request, send_file , send_from_directory
from flask_cors import CORS
import requests
from datetime import date, timedelta
import math
import csv
import io
import os
import datetime
import sqlite3
import json

app = Flask(__name__)
CORS(app)

@app.route('/')
def index():
    return send_from_directory('.', 'index.html')

HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
HYG_URL = "https://raw.githubusercontent.com/astronexus/HYG-Database/master/hyg/v3/hyg_v3.csv"
HYG_CACHE = os.path.join(os.path.dirname(__file__), "hyg_cache.csv")

TODAY     = date.today().isoformat()
TOMORROW  = (date.today() + timedelta(days=1)).isoformat()

BODIES = {
    "sun":     {"id": "10",  "name": "Sun"},
    "moon":    {"id": "301", "name": "Moon"},
    "venus":   {"id": "299", "name": "Venus"},
    "mars":    {"id": "499", "name": "Mars"},
    "jupiter": {"id": "599", "name": "Jupiter"},
    "saturn":  {"id": "699", "name": "Saturn"},
    "ceres":   {"id": "1",   "name": "Ceres"},
}

AU_TO_KM = 149_597_870.7

# ── Caches ──
_positions_cache = None
_stars_cache     = None
asteroid_cache = None
_tex_cache = {}

# ── Database ──
DB_PATH = os.path.join(os.path.dirname(__file__), 'messages.db')

def init_db():
    conn = sqlite3.connect(DB_PATH)
    conn.execute('''
        CREATE TABLE IF NOT EXISTS messages (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            destination TEXT NOT NULL,
            speed TEXT NOT NULL,
            progress REAL NOT NULL,
            message TEXT NOT NULL,
            created_at INTEGER NOT NULL
        )
    ''')
    conn.commit()
    conn.close()
    print("[DB] Messages database ready.")


# ══════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════

def dist_km(x, y, z):
    return math.sqrt(x**2 + y**2 + z**2) * AU_TO_KM

def fmt_dist(km):
    if km >= 9.461e12: return f"{km/9.461e12:.2f} light-years"
    if km >= 1e9:      return f"{km/1e9:.1f} billion km"
    if km >= 1e6:      return f"{round(km/1e6)} million km"
    return f"{round(km):,} km"
def current_jd():
    now = datetime.datetime.utcnow()
    a = (14 - now.month) // 12
    y = now.year + 4800 - a
    m = now.month + 12*a - 3
    jdn = now.day + (153*m+2)//5 + 365*y + y//4 - y//100 + y//400 - 32045
    return jdn + (now.hour-12)/24.0 + now.minute/1440.0 + now.second/86400.0

def solve_kepler(M, e):
    E = M
    for _ in range(50):
        dE = (M - E + e*math.sin(E)) / (1.0 - e*math.cos(E))
        E += dE
        if abs(dE) < 1e-10:
            break
    return E

def helio_xyz(a, e, i_deg, om_deg, w_deg, ma_deg, epoch_jd, jd):
    i  = math.radians(i_deg)
    om = math.radians(om_deg)
    w  = math.radians(w_deg)
    n  = 0.9856076686 / (a**1.5)   # mean motion deg/day (Kepler's 3rd law)
    M  = math.radians((ma_deg + n*(jd - epoch_jd)) % 360)
    E  = solve_kepler(M, e)
    nu = math.atan2(math.sqrt(1-e*e)*math.sin(E), math.cos(E)-e)
    r  = a*(1 - e*math.cos(E))
    xp, yp = r*math.cos(nu), r*math.sin(nu)
    co,so = math.cos(om),math.sin(om)
    cw,sw = math.cos(w), math.sin(w)
    ci,si = math.cos(i), math.sin(i)
    x = (co*cw - so*sw*ci)*xp + (-co*sw - so*cw*ci)*yp
    y = (so*cw + co*sw*ci)*xp + (-so*sw + co*cw*ci)*yp
    z = (si*sw)*xp             + (si*cw)*yp
    return x, y, z


# ══════════════════════════════════════════
# PLANET POSITIONS  (Stage 1)
# ══════════════════════════════════════════

def fetch_vector(body_id: str) -> dict | None:
    params = {
        "format":     "json",
        "COMMAND":    f"'{body_id}'",
        "CENTER":     "'500@10'",
        "MAKE_EPHEM": "'YES'",
        "TABLE_TYPE": "'VECTORS'",
        "START_TIME": f"'{TODAY}'",
        "STOP_TIME":  f"'{TOMORROW}'",
        "STEP_SIZE":  "'1d'",
        "QUANTITIES": "'1'",
        "OUT_UNITS":  "'AU-D'",
        "VEC_TABLE":  "'2'",
        "REF_PLANE":  "'ECLIPTIC'",
    }
    try:
        r = requests.get(HORIZONS_URL, params=params, timeout=12)
        r.raise_for_status()
        data = r.json()
        raw  = data.get("result", "")
        in_data = False
        for line in raw.split("\n"):
            if "$$SOE" in line: in_data = True; continue
            if "$$EOE" in line: break
            if in_data and ("X =" in line or "X=" in line):
                parts = line.replace("X","").replace("Y","").replace("Z","").replace("=","").split()
                nums  = [float(p) for p in parts if _is_float(p)]
                if len(nums) >= 3:
                    return {"x": nums[0], "y": nums[1], "z": nums[2]}
    except Exception as e:
        print(f"  JPL error for {body_id}: {e}")
    return None

def _is_float(s):
    try: float(s); return True
    except: return False


def get_positions():
    global _positions_cache
    if _positions_cache: return _positions_cache

    print(f"\n[Positions] Fetching from NASA JPL for {TODAY} ...")
    earth = fetch_vector("399") or {"x": 1.0, "y": 0.0, "z": 0.0}
    ex, ey, ez = earth["x"], earth["y"], earth["z"]

    result = {}
    for key, body in BODIES.items():
        print(f"  → {body['name']} ...", end=" ")
        vec = fetch_vector(body["id"])
        if vec:
            if key == "sun":
                rx, ry, rz = -ex, -ey, -ez
            else:
                rx = vec["x"] - ex
                ry = vec["y"] - ey
                rz = vec["z"] - ez
            dk = dist_km(rx, ry, rz)
            result[key] = {
                "name":     body["name"],
                "x": round(rx, 6), "y": round(ry, 6), "z": round(rz, 6),
                "dist_km":  round(dk, 0),
                "dist_str": fmt_dist(dk),
                "live":     True,
            }
            print(f"x={rx:.3f} y={ry:.3f} z={rz:.3f} AU  ({fmt_dist(dk)})")
        else:
            print("FAILED")

    # Proxima Centauri — fixed
    # Earth at origin — reference frame center
    result["earth"] = {
        "name":     "Earth",
        "x": 0.0, "y": 0.0, "z": 0.0,
        "dist_km":  149597870,
        "dist_str": "0 km",
        "live":     True,
    }
    pc = 4.24 * 206265  # AU  (1 ly = 63241 AU, 4.24 ly)
    result["proxima"] = {
        "name":     "Proxima Centauri",
        "x": round(pc * math.cos(math.radians(-60)) * math.cos(math.radians(217)), 0),
        "y": round(pc * math.cos(math.radians(-60)) * math.sin(math.radians(217)), 0),
        "z": round(pc * math.sin(math.radians(-60)), 0),
        "dist_km":  40_141_000_000_000,
        "dist_str": "4.24 light-years",
        "live":     False,
    }

    # Asteroid belt — midpoint of Mars & Jupiter
    if "mars" in result and "jupiter" in result:
        mx = (result["mars"]["x"] + result["jupiter"]["x"]) / 2
        my = (result["mars"]["y"] + result["jupiter"]["y"]) / 2
        mz = (result["mars"]["z"] + result["jupiter"]["z"]) / 2
        result["asteroid_belt"] = {
            "name":     "Asteroid Belt",
            "x": round(mx, 4), "y": round(my, 4), "z": round(mz, 4),
            "dist_km":  round(dist_km(mx, my, mz), 0),
            "dist_str": fmt_dist(dist_km(mx, my, mz)),
            "inner_au": 2.2, "outer_au": 3.2,
            "live":     True,
            "note":     "Between Mars and Jupiter — 99.9% empty space",
        }

    result["_meta"] = {
        "date":          TODAY,
        "source":        "NASA JPL Horizons",
        "reference":     "Earth-centered geocentric ecliptic",
        "units":         "AU for position, km for distance",
        "earth_helio_x": round(ex, 6),
        "earth_helio_y": round(ey, 6),
        "earth_helio_z": round(ez, 6),
    }

    _positions_cache = result
    print("[Positions] Done.\n")
    return result


# ══════════════════════════════════════════
# STAR COLORS  (B-V → RGB)
# ══════════════════════════════════════════

def bv_to_rgb(bv: float) -> tuple[float, float, float]:
    """
    Convert B-V color index to approximate RGB (0-1 each).
    Based on stellar temperature mapping.
    O/B stars: blue-white
    A stars:   white
    F stars:   yellow-white
    G stars:   yellow  (Sun is ~0.65)
    K stars:   orange
    M stars:   red-orange
    """
    if bv < -0.40: return (0.60, 0.70, 1.00)   # very hot blue
    if bv < -0.20: return (0.68, 0.78, 1.00)   # hot blue-white
    if bv < 0.00:  return (0.80, 0.88, 1.00)   # blue-white
    if bv < 0.15:  return (0.95, 0.97, 1.00)   # white
    if bv < 0.30:  return (1.00, 1.00, 0.95)   # white-yellow
    if bv < 0.50:  return (1.00, 1.00, 0.85)   # yellow-white (like Sun)
    if bv < 0.70:  return (1.00, 0.97, 0.75)   # yellow
    if bv < 0.90:  return (1.00, 0.90, 0.60)   # yellow-orange
    if bv < 1.20:  return (1.00, 0.78, 0.45)   # orange
    if bv < 1.50:  return (1.00, 0.65, 0.35)   # orange-red
    return             (1.00, 0.55, 0.30)       # red (cool M stars)


def mag_to_size(mag: float) -> float:
    """Map apparent magnitude to point size. Brighter = larger."""
    if mag < -1:  return 5.0
    if mag < 0:   return 4.0
    if mag < 1:   return 3.5
    if mag < 2:   return 3.0
    if mag < 3:   return 2.5
    if mag < 4:   return 2.0
    if mag < 5:   return 1.5
    return               1.0


# ══════════════════════════════════════════
# HYG STAR DATABASE  (Stage 2)
# ══════════════════════════════════════════

def download_hyg() -> str:
    import zipfile
    """Use local ZIP if it exists to bypass 50MB download limits."""
    zip_path = os.path.join(os.path.dirname(__file__), "hyg_cache.csv.zip")
    
    if os.path.exists(zip_path):
        print("[Stars] Using ZIPPED cached HYG database.")
        with zipfile.ZipFile(zip_path, "r") as z:
            # Find the CSV inside the zip file
            for name in z.namelist():
                if name.endswith('.csv'):
                    with z.open(name) as f:
                        return f.read().decode("utf-8")

    # Fallback just in case
    if os.path.exists(HYG_CACHE):
        print("[Stars] Using unzipped cached HYG database.")
        with open(HYG_CACHE, "r", encoding="utf-8") as f:
            return f.read()

    print("[Stars] Downloading HYG database (~50MB) ...")
    r = requests.get(HYG_URL, timeout=60)
    r.raise_for_status()
    return r.text


def get_stars(n: int = 9000):
    global _stars_cache
    if _stars_cache: return _stars_cache

    print(f"\n[Stars] Loading top {n} brightest stars from HYG database ...")
    try:
        content = download_hyg()
    except Exception as e:
        print(f"[Stars] Download failed: {e} — returning empty star list")
        _stars_cache = []
        return _stars_cache

    reader = csv.DictReader(io.StringIO(content))
    raw_stars = []

    for row in reader:
        try:
            mag  = float(row.get("mag", "99"))
            dist = float(row.get("dist", "0"))
            x    = float(row.get("x",   "0"))
            y    = float(row.get("y",   "0"))
            z    = float(row.get("z",   "0"))
            ci   = row.get("ci", "").strip()

            # Skip: unknown distance, mag too dim, or origin (row 0 = Sun)
            if dist <= 0 or dist >= 100000: continue
            if mag > 8.0: continue
            if x == 0 and y == 0 and z == 0: continue

            bv = float(ci) if ci else 0.65  # default Sun-like if unknown

            raw_stars.append({
                "mag": mag, "dist": dist,
                "x": x, "y": y, "z": z,
                "bv": bv,
                "proper": row.get("proper", "").strip(),
            })
        except (ValueError, TypeError):
            continue

    # Sort by magnitude (ascending = brightest first), keep top n
    raw_stars.sort(key=lambda s: s["mag"])
    top_stars = raw_stars[:n]
    print(f"[Stars] Parsed {len(raw_stars)} valid stars, keeping top {len(top_stars)}.")

    # Convert to compact format for JSON response
    stars_out = []
    for s in top_stars:
        r, g, b = bv_to_rgb(s["bv"])
        stars_out.append({
            "x":    round(s["x"], 4),   # parsecs from Sun (≈ from Earth for our scale)
            "y":    round(s["y"], 4),
            "z":    round(s["z"], 4),
            "dist": round(s["dist"], 2),  # parsecs
            "mag":  round(s["mag"], 2),
            "r":    round(r, 3),
            "g":    round(g, 3),
            "b":    round(b, 3),
            "size": mag_to_size(s["mag"]),
            "name": s["proper"],
        })

    _stars_cache = stars_out
    print(f"[Stars] Done. {len(stars_out)} stars ready.\n")
    return stars_out


# ══════════════════════════════════════════
# ROUTES
# ══════════════════════════════════════════

@app.route("/positions", methods=["GET"])
def api_positions():
    return jsonify(get_positions())

@app.route("/stars", methods=["GET"])
def api_stars():
    stars = get_stars(9000)
    return jsonify({
        "count":  len(stars),
        "units":  "parsecs (x,y,z), apparent magnitude (mag), RGB 0-1 (r,g,b)",
        "coords": "Heliocentric ICRS — Sun at origin",
        "note":   "Parallax multiplier applied in frontend for artistic effect",
        "stars":  stars,
    })

@app.route("/distances", methods=["GET"])
def api_distances():
    data = get_positions()
    dests = {}
    for key, obj in data.items():
        if key.startswith("_"): continue
        dests[key] = {"name": obj["name"], "km": obj.get("dist_km", 0), "live": obj.get("live", False)}
    return jsonify({"date": data["_meta"]["date"], "source": data["_meta"]["source"], "destinations": dests})

@app.route('/messages', methods=['GET'])
def get_messages():
    dest     = request.args.get('dest', '')
    speed    = request.args.get('speed', '')
    progress = float(request.args.get('progress', 0))
    radius   = float(request.args.get('radius', 0.02))
    try:
        conn = sqlite3.connect(DB_PATH)
        conn.row_factory = sqlite3.Row
        rows = conn.execute('''
            SELECT id, destination, speed, progress, message, created_at
            FROM messages
            WHERE destination = ?
              AND speed = ?
              AND progress >= ?
              AND progress <= ?
            ORDER BY ABS(progress - ?) ASC
            LIMIT 3
        ''', (dest, speed, progress-radius, progress+radius, progress)).fetchall()
        conn.close()
        now = int(datetime.datetime.utcnow().timestamp())
        out = []
        for r in rows:
            age = now - r['created_at']
            out.append({
                'id':          r['id'],
                'destination': r['destination'],
                'speed':       r['speed'],
                'progress':    round(r['progress'], 4),
                'message':     r['message'],
                'age_seconds': age,
            })
        return jsonify({'messages': out})
    except Exception as e:
        return jsonify({'error': str(e), 'messages': []}), 200


@app.route('/messages', methods=['POST'])
def post_message():
    try:
        data    = request.get_json()
        dest    = str(data.get('destination', ''))[:20]
        speed   = str(data.get('speed', ''))[:20]
        prog    = float(data.get('progress', 0))
        message = str(data.get('message', ''))[:280]
        if not dest or not message.strip():
            return jsonify({'error': 'missing fields'}), 400
        now = int(datetime.datetime.utcnow().timestamp())
        conn = sqlite3.connect(DB_PATH)
        conn.execute(
            'INSERT INTO messages (destination, speed, progress, message, created_at) VALUES (?,?,?,?,?)',
            (dest, speed, prog, message.strip(), now)
        )
        conn.commit()
        conn.close()
        return jsonify({'ok': True})
    except Exception as e:
        return jsonify({'error': str(e)}), 500
    
TEXTURE_URLS = {
    'earth':   'https://www.solarsystemscope.com/textures/download/2k_earth_daymap.jpg',
    'moon':    'https://raw.githubusercontent.com/mrdoob/three.js/master/examples/textures/planets/moon_1024.jpg',
    'venus':   'https://www.solarsystemscope.com/textures/download/2k_venus_atmosphere.jpg',
    'mars':    'https://www.solarsystemscope.com/textures/download/2k_mars.jpg',
    'jupiter': 'https://www.solarsystemscope.com/textures/download/2k_jupiter.jpg',
    'saturn':  'https://www.solarsystemscope.com/textures/download/2k_saturn.jpg',
    'proxima': 'https://www.solarsystemscope.com/textures/download/2k_sun.jpg',
}

@app.route('/texture/<name>')
def serve_texture(name):
    if name not in TEXTURE_URLS:
        return 'Not found', 404
    if name not in _tex_cache:
        print(f"[Texture] Downloading {name}...")
        r = requests.get(TEXTURE_URLS[name], timeout=15)
        _tex_cache[name] = r.content
        print(f"[Texture] {name} cached ({len(r.content)//1024}KB)")
    return send_file(io.BytesIO(_tex_cache[name]), mimetype='image/jpeg')


@app.route("/health", methods=["GET"])
def api_health():
    return jsonify({"status": "ok", "date": TODAY,
                    "positions_cached": _positions_cache is not None,
                    "stars_cached": _stars_cache is not None})

@app.route('/asteroids')
def get_asteroids():
    global asteroid_cache
    if asteroid_cache and asteroid_cache.get('count', 0) > 0:
        return jsonify(asteroid_cache)
    try:
        jd = current_jd()

        # Pull 2000 largest main-belt asteroids from JPL SBDB
        url = (
            'https://ssd-api.jpl.nasa.gov/sbdb_query.api'
            '?fields=full_name,a,e,i,om,w,ma,epoch'
            '&sb-kind=a&limit=5000&sort=H&full-prec=false'
        )
        r = requests.get(url, timeout=30)
        data = r.json()
        fields = data['fields']
        records = data.get('data', [])
        fi = {f: fields.index(f) for f in fields}

        # Convert heliocentric → geocentric
        # geocentric_asteroid = heliocentric_asteroid - heliocentric_Earth
        # heliocentric_Earth  = -geocentric_Sun (from your positions cache)
        ex, ey, ez = 0.0, 0.0, 0.0
        if _positions_cache:
            s = _positions_cache.get('sun', {})
            ex, ey, ez = -s.get('x',0), -s.get('y',0), -s.get('z',0)

        out = []
        for rec in records:
            try:
                a = float(rec[fi['a']])
                e = float(rec[fi['e']])
                if e >= 1.0 or a <= 0:
                    continue
                if a < 2.0 or a > 3.5:
                    continue
                hx,hy,hz = helio_xyz(
                    a, e,
                    float(rec[fi['i']]),
                    float(rec[fi['om']]),
                    float(rec[fi['w']]),
                    float(rec[fi['ma']]),
                    float(rec[fi['epoch']]),
                    jd
                )
                out.append({
                    'name': str(rec[fi['full_name']]).strip(),
                    'x': round(hx - ex, 4),
                    'y': round(hy - ey, 4),
                    'z': round(hz - ez, 4),
                })
            except Exception:
                continue

        asteroid_cache = {'asteroids': out, 'count': len(out)}
        return jsonify(asteroid_cache)

    except Exception as ex:
        # Do NOT cache errors — next request should retry
        return jsonify({'error': str(ex), 'asteroids': []}), 500


# ══════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════

if __name__ == "__main__":
    print("=" * 50)
    print("Space Simulator — Stage 2 Backend")
    print(f"Date : {TODAY}")
    print("Endpoints:")
    print("  /positions — real 3D planet positions (Earth-centered, AU)")
    print("  /stars     — 9,000 real stars with positions & colors")
    print("  /distances — legacy flat distances")
    print("  /health    — server status")
    print("=" * 50 + "\n")

    # Pre-load both datasets at startup
    init_db()
    get_positions()
    get_stars(9000)

    port = int(os.environ.get('PORT', 8080))
    app.run(host='0.0.0.0', port=port, debug=False)
