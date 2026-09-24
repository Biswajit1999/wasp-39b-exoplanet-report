"""Offline checksum verification for committed research inputs."""
from __future__ import annotations

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPECTRA = ROOT / "data" / "spectra"

# Windows CRLF working-tree hash and canonical-LF Git/archive-entry hash.
# Linux checkouts retain LF; Windows checkouts may convert to CRLF.
EXPECTED = {
    "eureka_transmission_spectrum.txt": (
        "976abf995f100ed34cd1fc0f94f777b0806a068896a2ed0510815a522594c243",
        "18ab790d131d28f4ad97d5a19dc3c5787ecc12269685823bbfd0d38dfe8619a8",
    ),
    "scchimeramodel.txt": (
        "9a376a2e2b966b0b8280906e8d559fb611b93deceb7e6a28549cf65b1eb6cec8",
        "45d014400577a5208f9f699a7811c0cc95aeb2c4f308583db9e90297f35cca1d",
    ),
    "scchimeramodel_no_co2.txt": (
        "176aefe865b59c30998af820fecc7c53683769532156bd68098285bda3ae0ef4",
        "7d41ab15efb4350ecce93f7e646d86cf4bd0fdab559f567d8e3dac51374ff004",
    ),
}


def verify() -> list[dict[str, object]]:
    rows = []
    for filename, (local_expected, archive_expected) in EXPECTED.items():
        raw = (SPECTRA / filename).read_bytes()
        checkout_actual = hashlib.sha256(raw).hexdigest()
        canonical = raw.replace(b"\r\n", b"\n")
        archive_actual = hashlib.sha256(canonical).hexdigest()
        rows.append({"file": filename,
                     "checkout_ok": checkout_actual in {local_expected, archive_expected},
                     "archive_content_ok": archive_actual == archive_expected,
                     "checkout_sha256": checkout_actual,
                     "canonical_lf_sha256": archive_actual})
    return rows


if __name__ == "__main__":
    results = verify()
    for row in results:
        print(f"{row['file']}: checkout={row['checkout_ok']} archive-content={row['archive_content_ok']}")
    if not all(row["checkout_ok"] and row["archive_content_ok"] for row in results):
        raise SystemExit(1)
