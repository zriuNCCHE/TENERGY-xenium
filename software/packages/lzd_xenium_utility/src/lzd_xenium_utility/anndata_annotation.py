from __future__ import annotations

import warnings
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from typing import Any, Literal


UnmatchedPolicy = Literal["drop", "keep"]


@dataclass
class AnnDataAnnotationConfig:
    """Default parameters for transferring AnnData annotations to SpatialData."""

    table_key: str = "table"
    sdata_cell_id_key: str | None = "cell_id"
    adata_cell_id_key: str | None = "cell_id"
    sample_key: str = "sample_id"
    sample_id: str | None = None
    annotation_columns: tuple[str, ...] | None = None
    unmatched_policy: UnmatchedPolicy = "drop"
    status_column: str | None = "annotation_match_status"
    annotated_status: str = "annotated"
    filtered_out_status: str = "filtered_out"
    overwrite: bool = True
    record_provenance: bool = True


class SpatialDataAnnotationMatcher:
    """
    Transfer annotations from a filtered AnnData object into a SpatialData table.

    The SpatialData table is treated as the destination. By default, unmatched
    raw SpatialData cells are dropped so the result follows the filtered,
    annotated AnnData. Set ``unmatched_policy="keep"`` to keep raw cells and
    mark cells absent from AnnData as filtered out.
    """

    def __init__(
        self,
        sdata: Any,
        adata: Any,
        config: AnnDataAnnotationConfig | None = None,
    ):
        self.sdata = sdata
        self.adata = adata
        self.config = config or AnnDataAnnotationConfig()
        self.provenance: list[dict[str, Any]] = []
        self._validate_inputs()

    def inspect_overlap(self) -> dict[str, int | str | None]:
        """Return match counts without modifying ``sdata``."""

        table = self._get_table()
        table_obs = table.obs
        source_obs = self._get_source_obs()
        table_keys = self._get_keys(table_obs, self.config.sdata_cell_id_key, "SpatialData table")
        source_keys = self._get_keys(source_obs, self.config.adata_cell_id_key, "AnnData")

        table_key_set = set(table_keys)
        source_key_set = set(source_keys)
        matched = table_key_set & source_key_set

        return {
            "table_key": self.config.table_key,
            "sample_key": self.config.sample_key if self.config.sample_id is not None else None,
            "sample_id": self.config.sample_id,
            "sdata_cells": len(table_keys),
            "adata_cells": len(source_keys),
            "matched_cells": len(matched),
            "sdata_only_cells": len(table_key_set - source_key_set),
            "adata_only_cells": len(source_key_set - table_key_set),
        }

    def transfer_annotations(
        self,
        annotation_columns: tuple[str, ...] | list[str] | None = None,
        unmatched_policy: UnmatchedPolicy | None = None,
        overwrite: bool | None = None,
    ) -> dict[str, int | str | None]:
        """Transfer annotation columns into ``sdata.tables[table_key].obs``."""

        resolved_policy = unmatched_policy or self.config.unmatched_policy
        if resolved_policy not in {"drop", "keep"}:
            raise ValueError("unmatched_policy must be either 'drop' or 'keep'.")

        resolved_overwrite = self.config.overwrite if overwrite is None else overwrite
        table = self._get_table()
        source_obs = self._get_source_obs()
        columns = self._resolve_annotation_columns(source_obs, annotation_columns)
        destination_columns = columns
        if self.config.status_column is not None and self.config.status_column not in destination_columns:
            destination_columns = (*destination_columns, self.config.status_column)
        self._validate_destination_columns(table.obs, destination_columns, resolved_overwrite)

        table_keys = self._get_keys(table.obs, self.config.sdata_cell_id_key, "SpatialData table")
        source_keys = self._get_keys(source_obs, self.config.adata_cell_id_key, "AnnData")
        self._validate_unique(source_keys, "AnnData")
        summary = self.inspect_overlap()
        self._warn_if_adata_only_cells(summary)

        source_by_key = source_obs.loc[:, columns].copy()
        source_by_key.index = source_keys
        matched_mask = table_keys.isin(source_by_key.index)

        working_table = table
        working_keys = table_keys
        if resolved_policy == "drop":
            working_table = self._subset_table(table, matched_mask)
            working_keys = table_keys[matched_mask]
            self.sdata.tables[self.config.table_key] = working_table

        for column in columns:
            working_table.obs[column] = list(working_keys.map(source_by_key[column]))

        if self.config.status_column is not None:
            status_values = working_keys.map(
                lambda key: (
                    self.config.annotated_status
                    if key in source_by_key.index
                    else self.config.filtered_out_status
                )
            )
            working_table.obs[self.config.status_column] = list(status_values)

        summary["unmatched_policy"] = resolved_policy
        summary["transferred_columns"] = len(columns)
        self._record_provenance(summary=summary, annotation_columns=columns)
        return summary

    def _validate_inputs(self) -> None:
        if not hasattr(self.sdata, "tables"):
            raise TypeError("sdata must be a SpatialData-like object with a 'tables' mapping.")
        if self.config.table_key not in self.sdata.tables:
            raise KeyError(f"sdata.tables does not contain table key '{self.config.table_key}'.")
        table = self._get_table()
        if not hasattr(table, "obs"):
            raise TypeError("sdata table must be AnnData-like and provide an 'obs' DataFrame.")
        if not hasattr(self.adata, "obs"):
            raise TypeError("adata must be AnnData-like and provide an 'obs' DataFrame.")
        if self.config.sample_id is not None and self.config.sample_key not in self.adata.obs.columns:
            raise KeyError(
                f"sample_id was provided, but AnnData obs has no sample_key column "
                f"'{self.config.sample_key}'."
            )

    def _get_table(self) -> Any:
        return self.sdata.tables[self.config.table_key]

    def _get_source_obs(self) -> Any:
        source_obs = self.adata.obs
        if self.config.sample_id is None:
            return source_obs
        sample_mask = source_obs[self.config.sample_key] == self.config.sample_id
        return source_obs.loc[sample_mask].copy()

    @staticmethod
    def _get_keys(obs: Any, key: str | None, label: str) -> Any:
        if key is None:
            keys = obs.index
        elif key in obs.columns:
            keys = obs[key]
        else:
            raise KeyError(
                f"{label} obs does not contain cell id column '{key}'. "
                "Pass the correct key or use None to match by obs_names."
            )
        if keys.isna().any():
            raise ValueError(f"{label} cell ids contain missing values.")
        return keys.astype(str)

    @staticmethod
    def _validate_unique(keys: Any, label: str) -> None:
        duplicate_mask = keys.duplicated()
        if duplicate_mask.any():
            duplicate_examples = sorted(set(keys[duplicate_mask]))[:5]
            raise ValueError(
                f"{label} cell ids must be unique after sample filtering. "
                f"Duplicate examples: {duplicate_examples}"
            )

    def _resolve_annotation_columns(
        self,
        source_obs: Any,
        annotation_columns: tuple[str, ...] | list[str] | None,
    ) -> tuple[str, ...]:
        if annotation_columns is not None:
            columns = tuple(annotation_columns)
        elif self.config.annotation_columns is not None:
            columns = tuple(self.config.annotation_columns)
        else:
            excluded = {self.config.adata_cell_id_key}
            columns = tuple(column for column in source_obs.columns if column not in excluded)

        missing = [column for column in columns if column not in source_obs.columns]
        if missing:
            raise KeyError(f"AnnData obs is missing annotation columns: {missing}")
        if not columns:
            raise ValueError("No annotation columns were selected for transfer.")
        return columns

    def _validate_destination_columns(
        self,
        destination_obs: Any,
        columns: tuple[str, ...],
        overwrite: bool,
    ) -> None:
        existing = [column for column in columns if column in destination_obs.columns]
        if existing and not overwrite:
            raise KeyError(
                f"Destination obs already contains columns {existing}. "
                "Pass overwrite=True to replace them."
            )

    @staticmethod
    def _warn_if_adata_only_cells(summary: dict[str, int | str | None]) -> None:
        adata_only_cells = summary.get("adata_only_cells", 0)
        if isinstance(adata_only_cells, int) and adata_only_cells > 0:
            warnings.warn(
                "AnnData contains "
                f"{adata_only_cells} cell IDs that are not present in the SpatialData table. "
                "These AnnData-only cells cannot be transferred.",
                UserWarning,
                stacklevel=3,
            )

    @staticmethod
    def _subset_table(table: Any, mask: Any) -> Any:
        mask_values = mask.to_numpy() if hasattr(mask, "to_numpy") else mask
        try:
            return table[mask_values, :].copy()
        except TypeError:
            return table[mask_values].copy()

    def _record_provenance(
        self,
        summary: dict[str, int | str | None],
        annotation_columns: tuple[str, ...],
    ) -> None:
        if not self.config.record_provenance:
            return

        record = {
            "utility": "SpatialDataAnnotationMatcher",
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "config": asdict(self.config),
            "annotation_columns": list(annotation_columns),
            "summary": dict(summary),
        }
        self.provenance.append(record)
        self._sync_latest_provenance_to_sdata()

    def _sync_latest_provenance_to_sdata(self) -> None:
        if not self.config.record_provenance or not self.provenance:
            return
        if not hasattr(self.sdata, "attrs") or self.sdata.attrs is None:
            return
        try:
            package_records = self.sdata.attrs.setdefault("lzd_xenium_utility", {})
            annotation_records = package_records.setdefault("anndata_annotation", [])
            annotation_records.append(dict(self.provenance[-1]))
        except AttributeError:
            return
