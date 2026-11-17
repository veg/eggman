/* tslint:disable */
/* eslint-disable */

export class WasmEggman {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Inspect FASTA records and automatically infer whether the dataset is HIV-1, IAV, SARS-CoV-2, or Sarbecovirus.
     */
    detect_pathogen(fasta_text: string): string;
    /**
     * Export any built-in pack (`hiv1`, `iav`, `sc2`, `sarbeco`) or loaded `"custom"` pack as a self-contained JSON bundle.
     */
    export_pack_bundle_json(pack_id: string): string;
    /**
     * Return the embedded HIV-1 200-lineage reference compendium as a JSON array.
     */
    get_compendium_json(): string;
    /**
     * Return metadata and reference compendium JSON for a built-in pack (`hiv1`, `iav`, `sc2`, `sarbeco`).
     */
    get_pack_meta_json(pack_id: string): string;
    /**
     * Load a custom `PathogenProfilePack` JSON (or bare `PathogenSpec` + optional FASTA profiles string) into slot `"custom"`.
     */
    load_custom_pack(json_str: string, fasta_override: string): string;
    constructor();
    /**
     * Legacy HIV-1 FASTA screening entry point preserved for backward compatibility.
     */
    screen_fasta(fasta_text: string, gamma: number, min_span_nt: number, min_informative_snps: number, p_fisher_thresh: number, mask_hypervariable_loops: boolean): string;
    /**
     * Screen a FASTA string using `"auto"` pathogen detection or an explicit pack id,
     * automatically assembling multi-record Influenza A segments (`seg_mode`: `"auto"`, `"single_isolate"`, `"per_record"`).
     */
    screen_with_pack(pack_id: string, fasta_text: string, seg_mode: string, gamma: number, min_span_nt: number, min_informative_snps: number, p_fisher_thresh: number, mask_hypervariable_loops: boolean): string;
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_wasmeggman_free: (a: number, b: number) => void;
    readonly wasmeggman_detect_pathogen: (a: number, b: number, c: number) => [number, number];
    readonly wasmeggman_export_pack_bundle_json: (a: number, b: number, c: number) => [number, number];
    readonly wasmeggman_get_compendium_json: (a: number) => [number, number];
    readonly wasmeggman_get_pack_meta_json: (a: number, b: number, c: number) => [number, number];
    readonly wasmeggman_load_custom_pack: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly wasmeggman_new: () => number;
    readonly wasmeggman_screen_fasta: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => [number, number];
    readonly wasmeggman_screen_with_pack: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number, k: number, l: number) => [number, number];
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
