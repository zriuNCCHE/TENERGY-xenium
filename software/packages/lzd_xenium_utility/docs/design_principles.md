# Design Principles

`LZD Xenium Utility` is intentionally a collection package rather than a single
linear analysis package. Each utility should be small, explicit, and easy to use
from a notebook or server script.

Core rules:

- one module per practical Xenium task,
- dataclass configuration for defaults,
- class-like utilities for stateful workflows,
- explicit file paths instead of hidden global variables,
- clear validation before expensive external calls,
- mutation only through named public methods,
- provenance for methods that modify `SpatialData`,
- small tests that mock heavy Xenium dependencies,
- stable public exports from `lzd_xenium_utility.__init__`.

The package should remain expandable by adding modules such as:

```text
cell_metadata.py
anndata_annotation.py
spatialdata_io_helpers.py
transcript_qc.py
image_qc.py
```

Avoid turning the package into one large miscellaneous file. If a utility needs
its own config, tests, and documentation section, it deserves its own module.
