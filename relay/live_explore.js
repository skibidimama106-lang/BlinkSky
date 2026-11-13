/* ---------- v10: live Explorer in the browser ----------
   Any RA/Dec -> a "little video" of SPHEREx photos, with no helper server.
   NASA/IPAC IRSA sends no cross-site permission header, so requests go through the BlinkSky relay
   (a tiny Cloudflare Worker, proxy_worker/worker.js) that only forwards SPHEREx searches, cutouts and Gaia lookups.
   This is a line-by-line port of data_access/blinksky_explorer.py (same steps, same numbers). */
const RELAY = 'https://blinksky-relay.a-slave-of-allah.workers.dev';
// DR1 (Nov 2026) replaces QR2 and is processed to match QR3. A new release = one more line here.
const LIVE_RELEASES = [['spherex_dr1', '10.26131/IRSA678'], ['spherex_qr3', '10.26131/IRSA662']];
const LIVE = { OUT_SCALE: 3.075, VISIT_GAP_DAYS: 30, SIZE: 0.1,
  BAD_BITS: { 0: 'TRANSIENT', 1: 'OVERFLOW', 2: 'SUR_ERROR', 6: 'NONFUNC', 9: 'MISSING_DATA', 10: 'HOT', 15: 'NONLINEAR',
              17: 'PERSIST', 19: 'OUTLIER', 26: 'BLOOM', 27: 'SNOWBALL', 28: 'HALO' } };

// ---- small number helpers (NaN-aware, like numpy's nan* functions) ----
const liveMedian = a => { const v = a.filter(Number.isFinite).sort((x, y) => x - y), n = v.length; return n ? (n % 2 ? v[(n - 1) / 2] : (v[n / 2 - 1] + v[n / 2]) / 2) : NaN; };
const livePct = (a, p) => { const v = a.filter(Number.isFinite).sort((x, y) => x - y); if (!v.length) return NaN; const k = (v.length - 1) * p / 100, f = Math.floor(k); return f + 1 < v.length ? v[f] + (v[f + 1] - v[f]) * (k - f) : v[f]; };
const D2R = Math.PI / 180;

// ---- FITS reader: just enough for SPHEREx cutouts (images + one binary table) ----
function fitsRead(buf) {
  const dv = new DataView(buf), u8 = new Uint8Array(buf), hdus = {}; let pos = 0;
  while (pos + 2880 <= u8.length) {
    const h = {}; let end = false;
    while (!end && pos + 2880 <= u8.length) {
      for (let i = 0; i < 36; i++) {
        let c = ''; for (let j = 0; j < 80; j++) c += String.fromCharCode(u8[pos + i * 80 + j]);
        const k = c.slice(0, 8).trim();
        if (k === 'END') { end = true; break; }
        if (c[8] !== '=') continue;
        let v = c.slice(10).trim();
        if (v[0] === "'") v = v.slice(1, v.indexOf("'", 1)).trim();
        else { v = v.split('/')[0].trim(); v = v === 'T' ? true : v === 'F' ? false : Number.isFinite(+v) ? +v : v; }
        h[k] = v;
      }
      pos += 2880;
    }
    const naxis = h.NAXIS || 0, dims = []; for (let i = 1; i <= naxis; i++) dims.push(h['NAXIS' + i]);
    const n = naxis ? dims.reduce((a, b) => a * b, 1) : 0, bytes = Math.abs(h.BITPIX) / 8 * (h.GCOUNT || 1) * ((h.PCOUNT || 0) + n);
    hdus[h.EXTNAME || (Object.keys(hdus).length ? 'HDU' + Object.keys(hdus).length : 'PRIMARY')] = { h, start: pos, dims, dv };
    pos += Math.ceil(bytes / 2880) * 2880;
  }
  return hdus;
}
function fitsImage(hdu) {
  const { h, start, dims, dv } = hdu, n = dims[0] * dims[1], out = new Float64Array(n), b = h.BITPIX;
  for (let i = 0; i < n; i++) out[i] = b === -32 ? dv.getFloat32(start + 4 * i) : b === -64 ? dv.getFloat64(start + 8 * i) : b === 32 ? dv.getInt32(start + 4 * i) : dv.getInt16(start + 2 * i);
  return out;   // row-major: index = y * nx + x (y = FITS row - 1)
}
function fitsTable(hdu) {     // first row of a binary table, every column as an array
  const { h, start, dv } = hdu, out = {}; let off = 0;
  for (let c = 1; c <= h.TFIELDS; c++) {
    const m = String(h['TFORM' + c]).match(/^(\d*)([A-Z])/), r = m[1] ? +m[1] : 1, t = m[2];
    const size = { J: 4, E: 4, D: 8, I: 2, K: 8, B: 1 }[t], vals = [];
    for (let i = 0; i < r; i++) { const p = start + off + i * size; vals.push(t === 'J' ? dv.getInt32(p) : t === 'E' ? dv.getFloat32(p) : t === 'D' ? dv.getFloat64(p) : t === 'I' ? dv.getInt16(p) : t === 'K' ? Number(dv.getBigInt64(p)) : dv.getUint8(p)); }
    out[h['TTYPE' + c]] = vals; off += r * size;
  }
  return out;
}

// ---- WCS: gnomonic (TAN) projection with SIP distortion, as in the SPHEREx headers. Pixels are 0-based like astropy. ----
function wcsMake(h) {
  const sip = p => { const o = h[p + '_ORDER'] || 0, t = []; for (let i = 0; i <= o; i++) for (let j = 0; j <= o - i; j++) { const v = h[`${p}_${i}_${j}`]; if (v) t.push([i, j, v]); } return t; };
  const poly = (t, x, y) => t.reduce((a, [i, j, v]) => a + v * x ** i * y ** j, 0);
  const c1 = h.CDELT1 || 1, c2 = h.CDELT2 || 1;
  const m = [(h.PC1_1 ?? 1) * c1, (h.PC1_2 || 0) * c1, (h.PC2_1 || 0) * c2, (h.PC2_2 ?? 1) * c2], det = m[0] * m[3] - m[1] * m[2];
  const A = sip('A'), B = sip('B'), AP = sip('AP'), BP = sip('BP'), a0 = h.CRVAL1 * D2R, d0 = h.CRVAL2 * D2R;
  return {
    pix2world(px, py) {
      const u = px + 1 - h.CRPIX1, v = py + 1 - h.CRPIX2, x = u + poly(A, u, v), y = v + poly(B, u, v);
      const xi = (m[0] * x + m[1] * y) * D2R, eta = (m[2] * x + m[3] * y) * D2R;
      const den = Math.cos(d0) - eta * Math.sin(d0);
      let ra = (a0 + Math.atan2(xi, den)) / D2R; ra = ((ra % 360) + 360) % 360;
      return [ra, Math.atan2(Math.sin(d0) + eta * Math.cos(d0), Math.hypot(xi, den)) / D2R];
    },
    world2pix(ra, dec) {
      const a = ra * D2R, d = dec * D2R, ca = Math.cos(a - a0), D = Math.sin(d) * Math.sin(d0) + Math.cos(d) * Math.cos(d0) * ca;
      const xi = Math.cos(d) * Math.sin(a - a0) / D / D2R, eta = (Math.sin(d) * Math.cos(d0) - Math.cos(d) * Math.sin(d0) * ca) / D / D2R;
      const x = (m[3] * xi - m[1] * eta) / det, y = (-m[2] * xi + m[0] * eta) / det;
      let u = x + poly(AP, x, y), v = y + poly(BP, x, y);
      for (let k = 0; k < 3; k++) { u = x - poly(A, u, v); v = y - poly(B, u, v); }   // polish so forward(u,v) = (x,y) exactly
      return [u + h.CRPIX1 - 1, v + h.CRPIX2 - 1];
    },
  };
}

// ---- 1. find the images (IRSA's image search, through the relay) ----
async function liveFindImages(ra, dec) {
  const out = [], seen = new Set();
  for (const [col, doi] of LIVE_RELEASES) {
    const r = await fetch(`${RELAY}/sia?collection=${col}&ra=${ra}&dec=${dec}`);
    if (!r.ok) throw new Error(`NASA's image search answered ${r.status}`);
    const xml = new DOMParser().parseFromString(await r.text(), 'text/xml');
    const fields = [...xml.getElementsByTagNameNS('*', 'FIELD')].map(f => f.getAttribute('name'));
    for (const tr of xml.getElementsByTagNameNS('*', 'TR')) {
      const tds = [...tr.getElementsByTagNameNS('*', 'TD')], x = {};
      fields.forEach((f, i) => x[f] = tds[i] ? tds[i].textContent : '');
      const id = x.obs_id + '|' + x.energy_bandpassname;     // the same exposure must not come twice from two releases
      if (x.dataproduct_subtype !== 'science' || seen.has(id)) continue;
      seen.add(id);
      out.push({ obs_id: x.obs_id, mjd: +x.t_min, detector: x.energy_bandpassname, url: x.access_url, release: col, doi });
    }
  }
  return out.sort((a, b) => a.mjd - b.mjd);
}
function liveVisits(images) {
  const v = [];
  for (const im of images) { if (v.length && im.mjd - v[v.length - 1].at(-1).mjd <= LIVE.VISIT_GAP_DAYS) v[v.length - 1].push(im); else v.push([im]); }
  return v;
}
function livePick(visits) {     // one photo per detector per visit, from the middle of the visit (quick mode)
  const chosen = [];
  visits.forEach((visit, vi) => {
    const by = {}; for (const im of visit) (by[im.detector] = by[im.detector] || []).push(im);
    for (const det of Object.keys(by).sort()) { const ims = by[det]; chosen.push({ ...ims[Math.floor(ims.length / 2)], visit: vi }); }
  });
  return chosen.sort((a, b) => a.mjd - b.mjd);
}

// ---- 2+3. one photo: re-grid it north-up and measure the brightest source near the centre ----
function liveWavelength(hd, tb, xf, yf) {
  const xs = tb.X, ys = tb.Y, V = tb.VALUES, nx = xs.length;
  const xr = xf - hd.CRPIX1W + hd.CRVAL1W, yr = yf - hd.CRPIX2W + hd.CRVAL2W;
  const ss = (arr, v) => { let i = 0; while (i < arr.length && arr[i] < v) i++; return i; };
  const ix = Math.min(Math.max(ss(xs, xr) - 1, 0), xs.length - 2), iy = Math.min(Math.max(ss(ys, yr) - 1, 0), ys.length - 2);
  const fx = (xr - xs[ix]) / (xs[ix + 1] - xs[ix]), fy = (yr - ys[iy]) / (ys[iy + 1] - ys[iy]);
  const at = (y, x, k) => V[(y * nx + x) * 2 + k];
  return [0, 1].map(k => at(iy, ix, k) * (1 - fx) * (1 - fy) + at(iy, ix + 1, k) * fx * (1 - fy) + at(iy + 1, ix, k) * (1 - fx) * fy + at(iy + 1, ix + 1, k) * fx * fy);
}
function liveMeasure(buf, ra, dec, size) {
  const F = fitsRead(buf), hd = F.IMAGE.h, nx = F.IMAGE.dims[0], ny = F.IMAGE.dims[1];
  const img = fitsImage(F.IMAGE), flags = fitsImage(F.FLAGS), vr = fitsImage(F.VARIANCE), w = wcsMake(hd), N = nx * ny;
  const [x0, y0] = w.world2pix(ra, dec);
  if (!(x0 >= 6 && x0 < nx - 6 && y0 >= 6 && y0 < ny - 6)) return null;
  const bad = new Uint8Array(N), badBits = Object.keys(LIVE.BAD_BITS).map(Number);
  for (let i = 0; i < N; i++) for (const b of badBits) if ((flags[i] >>> b) & 1) { bad[i] = 1; break; }
  const good = []; for (let i = 0; i < N; i++) if (!bad[i]) good.push(img[i]);
  const sky = good.length ? liveMedian(good) : liveMedian([...img]);
  const sig = good.length ? 1.4826 * liveMedian(good.map(v => Math.abs(v - sky))) : 1;
  const clean = img.map(v => Number.isFinite(v) ? v : sky), disp = Float64Array.from(clean);
  let hiOk = -Infinity, anyOver = false; for (let i = 0; i < N; i++) { if ((flags[i] >>> 1) & 1) anyOver = true; else if (clean[i] > hiOk) hiOk = clean[i]; }
  if (anyOver && hiOk > -Infinity) for (let i = 0; i < N; i++) if ((flags[i] >>> 1) & 1) disp[i] = hiOk;
  // the picture: the same north-up grid for every photo, so frames line up
  const n = Math.round(size * 3600 / LIVE.OUT_SCALE), pic = new Float64Array(n * n);
  for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
    const dE = ((n - 1) / 2 - i) * LIVE.OUT_SCALE / 3600, dN = ((n - 1) / 2 - j) * LIVE.OUT_SCALE / 3600, dg = dec + dN;
    const [px, py] = w.world2pix(ra + dE / Math.cos(dg * D2R), dg);
    if (!(px >= 0 && py >= 0 && px <= nx - 1 && py <= ny - 1)) { pic[j * n + i] = NaN; continue; }
    const ix = Math.min(Math.floor(px), nx - 2), iy = Math.min(Math.floor(py), ny - 2), fx = px - ix, fy = py - iy, k = iy * nx + ix;
    pic[j * n + i] = disp[k] * (1 - fx) * (1 - fy) + disp[k + 1] * fx * (1 - fy) + disp[k + nx] * (1 - fx) * fy + disp[k + nx + 1] * fx * fy;
  }
  // the measurement: brightest 3x3 sum within ~18 arcsec of the centre
  const cx = Math.round(x0), cy = Math.round(y0); let best = -Infinity, bx = cx, by = cy;
  for (let yy = cy - 3; yy <= cy + 3; yy++) for (let xx = cx - 3; xx <= cx + 3; xx++) {
    if (Math.hypot(xx - x0, yy - y0) > 3.2) continue;
    let s = 0; for (let a = yy - 1; a <= yy + 1; a++) for (let b = xx - 1; b <= xx + 1; b++) if (a >= 0 && b >= 0 && a < ny && b < nx) s += clean[a * nx + b];
    if (s > best) { best = s; bx = xx; by = yy; }
  }
  const ring = [];
  for (let yy = 0; yy < ny; yy++) for (let xx = 0; xx < nx; xx++) { const r = Math.hypot(xx - bx, yy - by), k = yy * nx + xx; if (r >= 6 && r <= 10 && !bad[k] && Number.isFinite(img[k])) ring.push(img[k]); }
  const bg = ring.length > 20 ? liveMedian(ring) : sky;
  let flux = 0, ev = 0, ps = 0, mx = 0, my = 0;
  for (let yy = Math.max(0, by - 2); yy < Math.min(ny, by + 3); yy++) for (let xx = Math.max(0, bx - 2); xx < Math.min(nx, bx + 3); xx++) {
    const k = yy * nx + xx, c = clean[k] - bg; flux += c; if (Number.isFinite(vr[k])) ev += vr[k];
    if (c > 0) { ps += c; mx += c * xx; my += c * yy; }
  }
  const err = Math.sqrt(ev) || sig * 5;
  if (ps > 0) { mx /= ps; my /= ps; } else { mx = bx; my = by; }
  const [sra, sdec] = w.pix2world(mx, my);
  let dra = sra - ra; if (dra > 180) dra -= 360; if (dra < -180) dra += 360;
  const warn = new Set();
  for (let yy = Math.max(0, by - 1); yy < Math.min(ny, by + 2); yy++) for (let xx = Math.max(0, bx - 1); xx < Math.min(nx, bx + 2); xx++)
    for (const b of badBits) if ((flags[yy * nx + xx] >>> b) & 1) warn.add(LIVE.BAD_BITS[b]);
  const [wave, width] = liveWavelength(hd, fitsTable(F['WCS-WAVE']), mx + 1, my + 1), pixSr = (6.15 / 206265) ** 2, r4 = v => Math.round(v * 1e4) / 1e4;
  return { pic, n, sky, sigma: sig, flux_mjy: r4(flux * pixSr * 1e9), err_mjy: r4(err * pixSr * 1e9), snr: err > 0 ? Math.round(flux / err * 10) / 10 : 0,
    off_e: Math.round(dra * Math.cos(dec * D2R) * 3600 * 100) / 100, off_n: Math.round((sdec - dec) * 3600 * 100) / 100,
    wave_um: r4(wave), width_um: r4(width), flags: [...warn].sort() };
}
function livePicture(pic, n, lo, hi) {    // same stretch as save_png() in the Python tool
  const cv = document.createElement('canvas'); cv.width = cv.height = n;
  const g = cv.getContext('2d'), id = g.createImageData(n, n);
  for (let i = 0; i < n * n; i++) {
    const v = Number.isFinite(pic[i]) ? pic[i] : lo, s = Math.asinh(Math.min(Math.max((v - lo) / Math.max(hi - lo, 1e-9), 0), 1) * 60) / Math.asinh(60), b = Math.floor(s * 255);
    id.data[i * 4] = id.data[i * 4 + 1] = id.data[i * 4 + 2] = b; id.data[i * 4 + 3] = 255;
  }
  g.putImageData(id, 0, 0); return cv;
}

// ---- 4. motion: straight line through the positions, one pass to drop outliers (fit_motion in the Python tool) ----
function liveFitMotion(frames) {
  let good = frames.filter(f => f.snr >= 8);
  if (good.length < 4) return null;
  const ts = good.map(f => f.mjd); if (Math.max(...ts) - Math.min(...ts) < 70) return null;
  const line = (t, y) => { const n = t.length, mt = t.reduce((a, b) => a + b, 0) / n, my = y.reduce((a, b) => a + b, 0) / n; let sxx = 0, sxy = 0; t.forEach((v, i) => { sxx += (v - mt) ** 2; sxy += (v - mt) * (y[i] - my); }); const s = sxy / sxx; return [s, my - s * mt]; };
  const years = g => { const t = g.map(f => f.mjd), m = t.reduce((a, b) => a + b, 0) / t.length; return t.map(v => (v - m) / 365.25); };
  for (let it = 0; it < 2; it++) {
    const t = years(good), dist = new Array(good.length).fill(0);
    for (const key of ['off_e', 'off_n']) { const y = good.map(f => f[key]), [s, c] = line(t, y); y.forEach((v, i) => dist[i] += (v - (s * t[i] + c)) ** 2); }
    const d = dist.map(Math.sqrt), lim = Math.max(3 * liveMedian(d) * 1.4826, 1), keep = d.map(v => v <= lim), nk = keep.filter(Boolean).length;
    if (nk === good.length || nk < 4) break;
    good = good.filter((f, i) => keep[i]);
  }
  const t = years(good), res = {};
  for (const key of ['off_e', 'off_n']) {
    const y = good.map(f => f[key]), [s, c] = line(t, y), rs = y.reduce((a, v, i) => a + (v - (s * t[i] + c)) ** 2, 0);
    const sc = Math.sqrt(rs / Math.max(1, y.length - 2)), mt = t.reduce((a, b) => a + b, 0) / t.length;
    res[key] = [s, sc / Math.sqrt(t.reduce((a, v) => a + (v - mt) ** 2, 0)), sc];
  }
  const [e, n] = [res.off_e, res.off_n], total = Math.hypot(e[0], n[0]), err = total > 0 ? Math.hypot(e[0] * e[1], n[0] * n[1]) / total : Math.hypot(e[1], n[1]);
  const r2 = v => Math.round(v * 100) / 100, tt = good.map(f => f.mjd);
  return { east_arcsec_per_yr: r2(e[0]), north_arcsec_per_yr: r2(n[0]), total_arcsec_per_yr: r2(total), error_arcsec_per_yr: r2(err),
    scatter_arcsec: r2(Math.hypot(e[2], n[2])), frames_used: good.length, frames_with_warning_flags: good.filter(f => f.flags.length).length,
    span_years: r2((Math.max(...tt) - Math.min(...tt)) / 365.25) };
}
const liveIso = mjd => new Date(Date.UTC(1858, 10, 17) + mjd * 864e5).toISOString().slice(0, 16).replace('T', ' ') + ' UTC';

const LIVE_CACHE = new Map();
async function liveExplore(ra, dec, progress, size = LIVE.SIZE) {
  const key = `${ra.toFixed(4)},${dec.toFixed(4)},${size}`;
  if (LIVE_CACHE.has(key)) return LIVE_CACHE.get(key);
  progress('Searching NASA/IPAC IRSA for every SPHEREx photo of this spot...');
  const images = await liveFindImages(ra, dec), visits = liveVisits(images), chosen = livePick(visits);
  if (!chosen.length) throw new Error('SPHEREx has no public photos of this spot yet.');
  progress(`${images.length} photos found in ${visits.length} visit${visits.length === 1 ? '' : 's'}. Cutting out ${chosen.length} small pieces...`);
  let done = 0; const results = [];
  const work = async im => {
    for (let attempt = 0; attempt < 3; attempt++) {      // NASA's cutout service sometimes answers 502 for a moment
      try {
        const r = await fetch(`${RELAY}/cutout?url=${encodeURIComponent(im.url)}&ra=${ra}&dec=${dec}&size=${size}`);
        if (!r.ok) throw new Error(r.status);
        const m = liveMeasure(await r.arrayBuffer(), ra, dec, size);
        if (m) results.push([im, m]);
        break;
      } catch (e) { if (attempt === 2) console.warn('skipped', im.obs_id, e); else await new Promise(ok => setTimeout(ok, 2000 * (attempt + 1))); }
    }
    progress(`Photo ${++done} of ${chosen.length}`);
  };
  const queue = [...chosen];
  await Promise.all(Array.from({ length: 6 }, async () => { while (queue.length) await work(queue.shift()); }));
  if (!results.length) throw new Error('No usable SPHEREx photo of this spot.');
  results.sort((a, b) => a[0].mjd - b[0].mjd);
  progress('Lining up the frames and measuring...');
  const hiBy = {};
  for (const det of new Set(results.map(([im]) => im.detector)))
    hiBy[det] = liveMedian(results.filter(([im]) => im.detector === det).map(([, m]) => livePct(Array.from(m.pic, v => v - m.sky), 99.8)));
  const frames = results.map(([im, m]) => ({
    img: livePicture(Array.from(m.pic, v => v - m.sky), m.n, -1.5 * m.sigma, Math.max(hiBy[im.detector], 8 * m.sigma)),
    mjd: Math.round(im.mjd * 1e5) / 1e5, time_utc: liveIso(im.mjd), visit: im.visit, detector: im.detector.slice(-2), obs_id: im.obs_id, release: im.release, file: im.url,
    flux_mjy: m.flux_mjy, err_mjy: m.err_mjy, snr: m.snr, off_e: m.off_e, off_n: m.off_n, wave_um: m.wave_um, width_um: m.width_um, flags: m.flags }));
  const byDet = {}; frames.forEach(f => (byDet[f.detector] = byDet[f.detector] || []).push(f.snr));
  const rank = d => [byDet[d].length >= 3 ? 1 : 0, liveMedian(byDet[d])];
  const bestDet = Object.keys(byDet).sort((a, b) => { const A = rank(a), B = rank(b); return B[0] - A[0] || B[1] - A[1]; })[0];
  const n = Math.round(size * 3600 / LIVE.OUT_SCALE);
  const data = { name: `RA ${ra.toFixed(4)}, Dec ${dec >= 0 ? '+' : ''}${dec.toFixed(4)}`, ra, dec, size_deg: size, picture_px: n, arcsec_per_px: LIVE.OUT_SCALE,
    source: 'NASA SPHEREx spectral images via NASA/IPAC IRSA', live: true,
    releases: LIVE_RELEASES.map(([c, d]) => ({ collection: c, doi: d, photos: images.filter(i => i.release === c).length })),
    photos_available: images.length,
    visits: visits.map(v => ({ from: liveIso(v[0].mjd).slice(0, 10), to: liveIso(v.at(-1).mjd).slice(0, 10), photos: v.length })),
    default_detector: bestDet, motion: liveFitMotion(frames),
    colours_collected: new Set(frames.filter(f => f.snr >= 3).map(f => f.detector + f.wave_um.toFixed(2))).size, frames };
  LIVE_CACHE.set(key, data);
  return data;
}
// Gaia DR3 stars in the picture, moved to 2026 with their own motion (blinksky_gaia_overlay.py, in the browser)
async function liveGaia(d, year = 2026.0, gmax = 18) {
  const r = await fetch(`${RELAY}/gaia?ra=${d.ra}&dec=${d.dec}&r=${(d.size_deg * 0.75 + 0.03).toFixed(4)}&gmax=${gmax}`);
  if (!r.ok) return null;
  const rows = (await r.text()).trim().split(/\r?\n/), head = rows.shift().split(','), yrs = year - 2016, half = d.picture_px * d.arcsec_per_px / 2, stars = [];
  for (const line of rows) {
    const c = line.split(','), x = Object.fromEntries(head.map((h, i) => [h, c[i]])), num = (v, z = 0) => Number.isFinite(parseFloat(v)) ? parseFloat(v) : z;
    const e = (num(x.ra) - d.ra) * Math.cos(d.dec * D2R) * 3600 + num(x.pmra) / 1000 * yrs, nn = (num(x.dec) - d.dec) * 3600 + num(x.pmdec) / 1000 * yrs;
    if (Math.abs(e) <= half && Math.abs(nn) <= half) stars.push([Math.round(e * 100) / 100, Math.round(nn * 100) / 100, Math.round(num(x.phot_g_mean_mag, 99) * 10) / 10, Math.hypot(num(x.pmra), num(x.pmdec)) / 1000, x.source_id]);
  }
  stars.sort((a, b) => a[2] - b[2]);
  return { source: 'Gaia DR3 (ESA/Gaia/DPAC) via NASA/IPAC IRSA', year, g_mag_limit: gmax, stars };
}

// AllWISE sources in the picture (NASA's WISE telescope, seen around 2010). Not moved: a fast mover sits away from its square today.
async function liveWise(d) {
  const r = await fetch(`${RELAY}/wise?ra=${d.ra}&dec=${d.dec}&r=${(d.size_deg * 0.75).toFixed(4)}`);
  if (!r.ok) return null;
  const rows = (await r.text()).trim().split(/\r?\n/), head = rows.shift().split(','), half = d.picture_px * d.arcsec_per_px / 2, out = [];
  for (const line of rows) {
    const c = line.split(','), x = Object.fromEntries(head.map((h, i) => [h, c[i]]));
    const e = (+x.ra - d.ra) * Math.cos(d.dec * D2R) * 3600, nn = (+x.dec - d.dec) * 3600;
    if (Math.abs(e) <= half && Math.abs(nn) <= half) out.push([Math.round(e * 100) / 100, Math.round(nn * 100) / 100, +x.w1mpro, +x.w2mpro, x.designation]);
  }
  return { source: 'AllWISE (NASA WISE, 2010-2011) via NASA/IPAC IRSA', stars: out };
}
// Aladin Lite (CDS Strasbourg) showing NASA's SPHEREx DR1 all-sky map around the spot, with our patch drawn on it.
const SPHEREX_HIPS = 'https://irsa.ipac.caltech.edu/data/hips/SPHEREx/spherex_allsky_cube';
let aladinMod = null, aladinView = null, aladinBox = null;
async function liveAladin(el, d) {
  if (!aladinMod) { aladinMod = (await import('https://cdn.jsdelivr.net/npm/aladin-lite@3.6.1/dist/aladin.js')).default; await aladinMod.init; }
  const A = aladinMod;
  if (!aladinView) {
    // the survey must be given here: switching it afterwards left only the coarsest tiles loading
    aladinView = A.aladin(el, { fov: 0.6, target: `${d.ra} ${d.dec}`, survey: SPHEREX_HIPS, showReticle: false, showCooGrid: false, showProjectionControl: false, showFrame: false, showLayersControl: false, showFullscreenControl: true });
    aladinBox = A.graphicOverlay({ color: '#6fd3d6', lineWidth: 2 }); aladinView.addOverlay(aladinBox);
  } else { aladinView.gotoRaDec(d.ra, d.dec); aladinView.setFoV(Math.max(0.6, d.size_deg * 5)); }
  aladinBox.removeAll();
  const h = d.picture_px * d.arcsec_per_px / 7200, c = Math.cos(d.dec * D2R);
  aladinBox.add(A.polygon([[d.ra - h / c, d.dec - h], [d.ra + h / c, d.dec - h], [d.ra + h / c, d.dec + h], [d.ra - h / c, d.dec + h]]));
}
