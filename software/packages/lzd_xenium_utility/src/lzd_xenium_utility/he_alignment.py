from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
import zipfile


@dataclass
class HEAlignmentConfig:
    """Default parameters for adding an aligned H&E image to SpatialData."""

    image_key: str = "aligned_he"
    alignment_filename: str = "matrix.csv"
    extract_dir: str | Path = "alignment_metadata"
    overwrite: bool = False
    record_provenance: bool = True


class HEImageAligner:
    """
    Align an H&E image to a Xenium SpatialData object.

    The aligner wraps ``spatialdata_io.xenium_aligned_image`` and stores the
    aligned image in ``sdata.images``. It accepts either a direct alignment CSV
    or a ZIP file containing the alignment matrix exported by Xenium Explorer or
    related 10x alignment tools.
    """

    def __init__(self, sdata: Any, config: HEAlignmentConfig | None = None):
        self.sdata = sdata
        self.config = config or HEAlignmentConfig()
        self.provenance: list[dict[str, Any]] = []
        self._validate_sdata()

    def align_from_files(
        self,
        he_image_path: str | Path,
        alignment_file: str | Path,
        image_key: str | None = None,
        overwrite: bool | None = None,
        _provenance_source: str = "files",
        _alignment_zip: str | Path | None = None,
    ) -> Any:
        """Create an aligned H&E image from explicit files and add it to ``sdata.images``."""

        he_image_path = self._require_file(he_image_path, "H&E image")
        alignment_file = self._require_file(alignment_file, "alignment file")
        resolved_image_key = image_key or self.config.image_key
        resolved_overwrite = self.config.overwrite if overwrite is None else overwrite
        self._validate_image_key(resolved_image_key, overwrite=resolved_overwrite)

        aligned_image = self._load_aligned_image(
            he_image_path=he_image_path,
            alignment_file=alignment_file,
        )
        self.sdata.images[resolved_image_key] = aligned_image

        self._record_provenance(
            image_key=resolved_image_key,
            he_image_path=he_image_path,
            alignment_file=alignment_file,
            source=_provenance_source,
            alignment_zip=Path(_alignment_zip) if _alignment_zip is not None else None,
        )
        return aligned_image

    def align_from_zip(
        self,
        he_image_path: str | Path,
        alignment_zip: str | Path,
        image_key: str | None = None,
        extract_dir: str | Path | None = None,
        alignment_filename: str | None = None,
        overwrite: bool | None = None,
    ) -> Any:
        """Extract an alignment matrix from a ZIP, align the H&E image, and add it."""

        he_image_path = self._require_file(he_image_path, "H&E image")
        alignment_zip = self._require_file(alignment_zip, "alignment ZIP")
        resolved_extract_dir = Path(extract_dir or self.config.extract_dir)
        resolved_alignment_filename = alignment_filename or self.config.alignment_filename
        alignment_file = self.extract_alignment_file(
            alignment_zip=alignment_zip,
            extract_dir=resolved_extract_dir,
            alignment_filename=resolved_alignment_filename,
        )

        aligned_image = self.align_from_files(
            he_image_path=he_image_path,
            alignment_file=alignment_file,
            image_key=image_key,
            overwrite=overwrite,
            _provenance_source="zip",
            _alignment_zip=alignment_zip,
        )
        return aligned_image

    def extract_alignment_file(
        self,
        alignment_zip: str | Path,
        extract_dir: str | Path | None = None,
        alignment_filename: str | None = None,
    ) -> Path:
        """Extract the configured alignment CSV from a ZIP file."""

        alignment_zip = self._require_file(alignment_zip, "alignment ZIP")
        resolved_extract_dir = Path(extract_dir or self.config.extract_dir)
        resolved_alignment_filename = alignment_filename or self.config.alignment_filename
        resolved_extract_dir.mkdir(parents=True, exist_ok=True)

        with zipfile.ZipFile(alignment_zip, "r") as zip_ref:
            members = [
                member
                for member in zip_ref.infolist()
                if not member.is_dir() and Path(member.filename).name == resolved_alignment_filename
            ]
            if not members:
                raise FileNotFoundError(
                    f"Could not find '{resolved_alignment_filename}' inside alignment ZIP: {alignment_zip}"
                )
            member = sorted(members, key=lambda item: item.filename)[0]
            target_path = resolved_extract_dir / resolved_alignment_filename
            with zip_ref.open(member) as source, target_path.open("wb") as target:
                target.write(source.read())

        return target_path

    def _load_aligned_image(self, he_image_path: Path, alignment_file: Path) -> Any:
        try:
            from spatialdata_io import xenium_aligned_image
        except ImportError as exc:
            raise ImportError(
                "spatialdata-io is required for HEImageAligner. Install this package "
                "in an environment that provides spatialdata-io."
            ) from exc

        return xenium_aligned_image(
            image_path=he_image_path,
            alignment_file=alignment_file,
        )

    def _validate_sdata(self) -> None:
        if not hasattr(self.sdata, "images"):
            raise TypeError("sdata must be a SpatialData-like object with an 'images' mapping.")
        if self.sdata.images is None:
            raise TypeError("sdata.images must be a mutable mapping, not None.")
        if not hasattr(self.sdata.images, "__setitem__"):
            raise TypeError("sdata.images must support item assignment.")

    def _validate_image_key(self, image_key: str, overwrite: bool) -> None:
        if not image_key:
            raise ValueError("image_key must be a non-empty string.")
        if image_key in self.sdata.images and not overwrite:
            raise KeyError(
                f"sdata.images already contains image key '{image_key}'. "
                "Pass overwrite=True to replace it."
            )

    @staticmethod
    def _require_file(path: str | Path, label: str) -> Path:
        resolved_path = Path(path)
        if not resolved_path.exists():
            raise FileNotFoundError(f"{label} does not exist: {resolved_path}")
        if not resolved_path.is_file():
            raise FileNotFoundError(f"{label} is not a file: {resolved_path}")
        return resolved_path

    def _record_provenance(
        self,
        image_key: str,
        he_image_path: Path,
        alignment_file: Path,
        source: str,
        alignment_zip: Path | None,
    ) -> None:
        if not self.config.record_provenance:
            return

        record = {
            "utility": "HEImageAligner",
            "image_key": image_key,
            "he_image_path": str(he_image_path),
            "alignment_file": str(alignment_file),
            "alignment_zip": str(alignment_zip) if alignment_zip is not None else None,
            "source": source,
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "config": asdict(self.config),
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
            he_records = package_records.setdefault("he_alignment", [])
            he_records.append(dict(self.provenance[-1]))
        except AttributeError:
            return
