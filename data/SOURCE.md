# Data sources

## TESS light curve

- File: `tess2022112184951-s0051-0000000181949561-0223-s_lc.fits`
- Archive: Mikulski Archive for Space Telescopes (MAST), TESS SPOC light-curve product
- TESS sector: 51
- TIC target ID: 181949561
- MAST observation ID: 87494793
- MAST data URI: `mast:TESS/product/tess2022112184951-s0051-0000000181949561-0223-s_lc.fits`
- Exact download URL: <https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:TESS%2Fproduct%2Ftess2022112184951-s0051-0000000181949561-0223-s_lc.fits>
- Collection DOI: [10.17909/t9-nmc8-f686](https://doi.org/10.17909/t9-nmc8-f686) (TESS 2-minute light curves, all sectors; sector 51 used here)
- Retrieved: 2026-08-15
- SHA-256: `dc75775fde0187da7f058354635a90529e59128ecbfe683112024ad3e0c04ce8`

The FITS file is stored unmodified. The analysis reads `TIME`, `PDCSAP_FLUX`,
`PDCSAP_FLUX_ERR`, and `QUALITY`. PDCSAP flux is the SPOC light curve with common
instrumental trends removed and aperture/crowding corrections applied; this does
not make it free of residual stellar or instrumental systematics.

## System parameters

- File: `system_parameters.csv`
- Service: NASA Exoplanet Archive TAP, `pscomppars` table
- Exact query: <https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name%2Chostname%2Cra%2Cdec%2Cpl_orbper%2Cpl_tranmid%2Cpl_trandur%2Cpl_rade%2Cpl_bmasse%2Cpl_eqt%2Cpl_orbsmax%2Csy_dist%2Csy_tmag%2Cst_teff%2Cst_rad%2Cst_mass%2Cdisc_year%2Cdiscoverymethod%2Cdisc_refname%2Cdisc_pubdate%2Cdisc_facility+from+pscomppars+where+pl_name%3D%27WASP-39+b%27&format=csv>
- Retrieved: 2026-08-15

The saved row is the input actually used by `scripts/analyze_transit.py`; the
analysis does not query a changing live service at run time.


## Additional TESS sectors for robustness analysis

All are unmodified standard-cadence SPOC light curves from the same [MAST TESS collection](https://doi.org/10.17909/t9-nmc8-f686).

- Sector 51: `tess2022112184951-s0051-0000000181949561-0223-s_lc.fits` (1,800,000 bytes)
  - MAST URI: `mast:TESS/product/tess2022112184951-s0051-0000000181949561-0223-s_lc.fits`
  - SHA-256: `dc75775fde0187da7f058354635a90529e59128ecbfe683112024ad3e0c04ce8`

## Published planetary spectrum

- Archive record: [10.5281/zenodo.6959427](https://zenodo.org/records/6959427)
- Archive object: `JWST_ERS_1st_LOOK_PAPER_DATA.zip` (375,091 bytes)
- Archive MD5 (as served by Zenodo): `578368eb0c86014462f109d1e8699693`
- Data type: transmission; instrument: JWST NIRSpec PRISM
- Retrieved and audited: 2026-09-24

| Committed path | Original ZIP entry | Committed SHA-256 | Original-entry SHA-256 |
|---|---|---|---|
| `data/spectra/eureka_transmission_spectrum.txt` | `ZENODO/TRANSMISSION_SPECTRA_DATA/EUREKA_REDUCTION.txt` | `976abf995f100ed34cd1fc0f94f777b0806a068896a2ed0510815a522594c243` | `18ab790d131d28f4ad97d5a19dc3c5787ecc12269685823bbfd0d38dfe8619a8` |
| `data/spectra/scchimeramodel.txt` | `ZENODO/MODEL_FITS/ScCHIMERA_MODEL.txt` | `9a376a2e2b966b0b8280906e8d559fb611b93deceb7e6a28549cf65b1eb6cec8` | `45d014400577a5208f9f699a7811c0cc95aeb2c4f308583db9e90297f35cca1d` |
| `data/spectra/scchimeramodel_no_co2.txt` | `ZENODO/MODEL_FITS/ScCHIMERA_MODEL_noCO2.txt` | `176aefe865b59c30998af820fecc7c53683769532156bd68098285bda3ae0ef4` | `7d41ab15efb4350ecce93f7e646d86cf4bd0fdab559f567d8e3dac51374ff004` |

The decoded lines match the named archive entries exactly. The committed byte
hashes differ only because the checkout uses CRLF line endings while the ZIP
entries use LF. `python scripts/verify_data_provenance.py` verifies both the
committed bytes and the LF-normalized archive identity without network access.

The archive README describes the ScCHIMERA files as best-fitting spectra and
the `noCO2` file as a “remove one gas at a time” spectrum used for its Figure 3.
That makes the repository comparison a sensitivity diagnostic between supplied
curves, not an independently reproduced retrieval or a nested hypothesis test.
