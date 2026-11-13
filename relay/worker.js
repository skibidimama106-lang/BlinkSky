// BlinkSky relay (Cloudflare Worker).
// NASA/IPAC IRSA does not send a cross-site permission header, so a public web page cannot ask it for SPHEREx
// photos directly. This relay does only that, and only for IRSA:
//   GET /sia?collection=spherex_qr2&ra=..&dec=..   -> IRSA's image search (VOTable XML)
//   GET /cutout?url=<IRSA image URL>&ra=..&dec=..&size=0.1 -> a small cutout of that image (FITS),
//        with the big point-spread-function layer (PSF or EPSF) and ZODI removed so ~4 MB becomes ~100 KB.
//   GET /gaia?ra=..&dec=..&r=..                    -> Gaia DR3 stars near the spot (via IRSA's TAP service, CSV)
//   GET /wise?ra=..&dec=..&r=..                    -> AllWISE sources near the spot (NASA WISE, via IRSA's TAP service, CSV)
const CORS = { "Access-Control-Allow-Origin": "*", "Access-Control-Allow-Methods": "GET, OPTIONS" };
const IRSA = "https://irsa.ipac.caltech.edu";
const COLLECTION = /^spherex_[a-z0-9_]{2,20}$/;   // any SPHEREx release, so new ones (DR1, QR4...) work with no change here
const DROP = new Set(["EPSF", "PSF", "ZODI"]);   // big tables the page never uses
const UA = { "User-Agent": "BlinkSky relay (NASA Space Apps 2026)" };

const num = (v, lo, hi) => { const x = Number(v); if (!Number.isFinite(x) || x < lo || x > hi) throw new Error("bad number"); return x; };
const reply = (body, status = 200, type = "text/plain") => new Response(body, { status, headers: { ...CORS, "Content-Type": type } });

// Walk the FITS file HDU by HDU (2880-byte blocks) and keep every HDU except the ones in DROP.
function dropEpsf(buf) {
  const u8 = new Uint8Array(buf), keep = []; let pos = 0;
  while (pos + 2880 <= u8.length) {
    const start = pos; const cards = {}; let end = false;
    while (!end && pos + 2880 <= u8.length) {
      const block = new TextDecoder("latin1").decode(u8.subarray(pos, pos + 2880)); pos += 2880;
      for (let i = 0; i < 36; i++) {
        const c = block.slice(i * 80, i * 80 + 80), k = c.slice(0, 8).trim();
        if (k === "END") { end = true; break; }
        if (c[8] === "=") cards[k] = c.slice(10).split("/")[0].trim().replace(/^'|'\s*$/g, "").trim();
      }
    }
    if (!end) return buf;
    const naxis = +cards.NAXIS || 0; let n = naxis ? 1 : 0;
    for (let i = 1; i <= naxis; i++) n *= +cards["NAXIS" + i];
    const bytes = Math.abs(+cards.BITPIX) / 8 * (+cards.GCOUNT || 1) * ((+cards.PCOUNT || 0) + n);
    pos += Math.ceil(bytes / 2880) * 2880;
    if (!DROP.has(cards.EXTNAME)) keep.push(u8.subarray(start, Math.min(pos, u8.length)));
  }
  const out = new Uint8Array(keep.reduce((a, k) => a + k.length, 0)); let o = 0;
  for (const k of keep) { out.set(k, o); o += k.length; }
  return out;
}

export default {
  async fetch(req) {
    if (req.method === "OPTIONS") return new Response(null, { headers: CORS });
    const u = new URL(req.url), q = u.searchParams;
    try {
      if (u.pathname === "/sia") {
        const col = q.get("collection"); if (!COLLECTION.test(col || "")) return reply("unknown collection", 400);
        const ra = num(q.get("ra"), 0, 360), dec = num(q.get("dec"), -90, 90);
        const r = await fetch(`${IRSA}/SIA?COLLECTION=${col}&POS=CIRCLE+${ra}+${dec}+0.01&MAXREC=3000`, { headers: UA, cf: { cacheTtl: 86400, cacheEverything: true } });
        return reply(await r.text(), r.status, "text/xml");
      }
      if (u.pathname === "/cutout") {
        const url = new URL(q.get("url") || "");
        if (url.origin !== IRSA || !url.pathname.startsWith("/ibe/data/spherex/")) return reply("only SPHEREx images on IRSA", 400);
        const ra = num(q.get("ra"), 0, 360), dec = num(q.get("dec"), -90, 90), size = num(q.get("size") || 0.1, 0.01, 0.5);
        const r = await fetch(`${url.origin}${url.pathname}?center=${ra},${dec}&size=${size}`, { headers: UA, cf: { cacheTtl: 2592000, cacheEverything: true } });
        if (!r.ok) return reply(`IRSA answered ${r.status}`, 502);
        const body = dropEpsf(await r.arrayBuffer());
        return new Response(body, { headers: { ...CORS, "Content-Type": "application/fits", "Cache-Control": "public, max-age=2592000" } });
      }
      if (u.pathname === "/gaia") {
        const ra = num(q.get("ra"), 0, 360), dec = num(q.get("dec"), -90, 90), rad = num(q.get("r") || 0.06, 0.001, 0.2);
        const gmax = num(q.get("gmax") || 18, 5, 21);
        const adql = `SELECT ra,dec,pmra,pmdec,phot_g_mean_mag,source_id FROM gaia_dr3_source WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',${ra},${dec},${rad}))=1 AND phot_g_mean_mag < ${gmax}`;
        const r = await fetch(`${IRSA}/TAP/sync?FORMAT=csv&QUERY=${encodeURIComponent(adql)}`, { headers: UA, cf: { cacheTtl: 2592000, cacheEverything: true } });
        return reply(await r.text(), r.status, "text/csv");
      }
      if (u.pathname === "/wise") {
        const ra = num(q.get("ra"), 0, 360), dec = num(q.get("dec"), -90, 90), rad = num(q.get("r") || 0.06, 0.001, 0.3);
        const adql = `SELECT TOP 500 ra,dec,w1mpro,w2mpro,designation FROM allwise_p3as_psd WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',${ra},${dec},${rad}))=1 AND (w1mpro < 16 OR w2mpro < 15) ORDER BY w2mpro`;
        const r = await fetch(`${IRSA}/TAP/sync?FORMAT=csv&QUERY=${encodeURIComponent(adql)}`, { headers: UA, cf: { cacheTtl: 2592000, cacheEverything: true } });
        return reply(await r.text(), r.status, "text/csv");
      }
      return reply("BlinkSky relay: /sia, /cutout, /gaia, /wise. Data: NASA/IPAC IRSA (SPHEREx, AllWISE), ESA Gaia.");
    } catch (e) {
      return reply(String(e.message || e), 400);
    }
  },
};
