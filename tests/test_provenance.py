from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import verify_data_provenance as provenance


def test_committed_spectra_match_local_and_archive_checksums():
    rows = provenance.verify()
    assert len(rows) == 3
    assert all(row["local_ok"] for row in rows)
    assert all(row["archive_content_ok"] for row in rows)
