# Rebuild the original Sam audit drivers

`generate_sam_audit.py` contains our extraction/patch rules, not the author's implementation. It requires the original public `metric_reconstruction_calc_radiative.m` at commit `5c1b42793ff893fd0c65a5b48ceecc02d34c5e35`, SHA-256 `124fc609fb36f46d4f28378ac040655c8ae7d3cff9346ec62134c358ce223ca9`.

From the repository root, after downloading that original source as documented in `docs/LORENZ_PUBLIC_CODE_SEARCH_20260915.md`:

```bash
python docs/environment_reproduction/wolfram/generate_sam_audit.py \
  --variant i6_h5 --output /path/to/local/run/sam_actual_L10_i6_h5
```

The generator verifies each reconstructed source block against the recorded run hash. Place the documented `bhpt_paths.wl` one directory above the output directory. On this machine the main `.wls` is placed on the Windows filesystem and launched with the available Wolfram Kernel 14.0 using `-noprompt -script <output>/run_sam.wls`. The generator itself never launches a numerical calculation.

After i6/h5 writes its completed `matching_state.mx`, generate the refinement:

```bash
python docs/environment_reproduction/wolfram/generate_sam_audit.py \
  --variant i7_h6 \
  --checkpoint /path/to/local/run/sam_actual_L10_i6_h5/matching_state.mx \
  --output /path/to/local/run/sam_actual_L10_i7_h6
```

A fresh DumpSave binary may have a different hash. The generator records the actual checkpoint hash and inserts it into the i7 runtime verification; all background/matching parameters are independently checked after loading. The original and refinement blocks contain the same physical equations. The only original matching-code edit preserves the supplied exact spin instead of rounding it to three decimals; the radial edits add logs. Both drivers preserve the complete L=10 matching but sample only r={2,10,30}.

The generator reads the adjacent `sam_actual_i6_h5_driver_manifest_20260915.json` and `sam_actual_i7_h6_driver_manifest_20260915.json` for exact source line spans, numerical configurations, and logging insertions. Our `run_sam_i6_h5.wls` and `run_sam_i7_h6.wls` are versioned orchestration templates; generated author blocks stay in the output folder and are not redistributed here.
