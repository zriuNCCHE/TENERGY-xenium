# Codex Project: Xenium Metaprogram Analysis

Purpose: continue development of the Xenium-adapted malignant-cell meta-program analysis pipeline.

Main code:
- `src/xenium_metaprogram/malignant_metaprogrammer.py`

Starting point:
- `examples/run_tenergy_mp.py`

Reference material:
- `references/original_MP.py`
- `references/gemini_handoff.txt`
- `references/example_similarity_heatmap.png`

Recommended next development steps:
1. Run the example on the real TENERGY malignant pre-treatment AnnData.
2. Add clinical-response violin/spatial plotting helpers.
3. Add synthetic AnnData tests for centering, NMF generation, clustering, and scoring.
4. Compare MP stability across `k_range`, `top_n_genes`, and consensus thresholds.
