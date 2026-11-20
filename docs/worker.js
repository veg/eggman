import init, { WasmEggman } from './pkg/eggman_wasm.js';

let engine = null;

async function ensureEngine() {
  if (!engine) {
    await init();
    engine = new WasmEggman();
  }
  return engine;
}

self.onmessage = async (ev) => {
  const msg = ev.data || {};
  const { reqId, type, payload } = msg;
  try {
    const wasm = await ensureEngine();
    if (type === 'init') {
      const metas = {};
      for (const id of ['hiv1', 'iav', 'sc2', 'sarbeco']) {
        metas[id] = JSON.parse(wasm.get_pack_meta_json(id));
      }
      self.postMessage({ reqId, ok: true, type: 'init', metas });
      return;
    }

    if (type === 'detect') {
      const detection = JSON.parse(wasm.detect_pathogen(payload.fastaText || ''));
      self.postMessage({ reqId, ok: true, type: 'detect', detection });
      return;
    }

    if (type === 'load_custom_pack') {
      const res = JSON.parse(wasm.load_custom_pack(payload.jsonStr || '', payload.fastaStr || ''));
      self.postMessage({ reqId, ok: true, type: 'load_custom_pack', ...res });
      return;
    }

    if (type === 'export_pack') {
      const bundle = wasm.export_pack_bundle_json(payload.packId || 'hiv1');
      self.postMessage({ reqId, ok: true, type: 'export_pack', bundle });
      return;
    }

    if (type === 'screen') {
      const t0 = performance.now();
      const raw = typeof wasm.screen_with_pack_and_temporal === 'function'
        ? wasm.screen_with_pack_and_temporal(
            payload.packId || 'auto',
            payload.fastaText || '',
            payload.segMode || 'auto',
            Number(payload.gamma || 0),
            Math.trunc(Number(payload.minSpan || 0)),
            Math.trunc(Number(payload.minSnps || 5)),
            Number(payload.pThresh || 1e-4),
            Boolean(payload.maskLoops),
            Boolean(payload.temporalGuard),
            payload.temporalMode || 'sampled',
            Number.isFinite(Number(payload.queryYear)) && Number(payload.queryYear) > 0 ? Number(payload.queryYear) : NaN,
            Number.isFinite(Number(payload.temporalSlackYr)) && Number(payload.temporalSlackYr) >= 0 ? Number(payload.temporalSlackYr) : NaN
          )
        : wasm.screen_with_pack(
            payload.packId || 'auto',
            payload.fastaText || '',
            payload.segMode || 'auto',
            Number(payload.gamma || 0),
            Math.trunc(Number(payload.minSpan || 0)),
            Math.trunc(Number(payload.minSnps || 5)),
            Number(payload.pThresh || 1e-4),
            Boolean(payload.maskLoops)
          );
      const elapsedMs = performance.now() - t0;
      const response = JSON.parse(raw);
      self.postMessage({ reqId, ok: true, type: 'screen', response, elapsedMs });
      return;
    }

    throw new Error(`Unknown worker action: ${type}`);
  } catch (err) {
    self.postMessage({
      reqId,
      ok: false,
      type,
      error: err && err.message ? err.message : String(err),
    });
  }
};
