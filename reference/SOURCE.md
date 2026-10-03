# Public reference provenance

The final model uses the original Replogle et al. (2022) genome-wide K562 raw
pseudobulk file, not the essential-gene differential-expression summary dataset.

- Article: https://plus.figshare.com/articles/dataset/20029387
- DOI: 10.25452/figshare.plus.20029387.v1
- File: `K562_gwps_raw_bulk_01.h5ad`, Figshare file ID 35774443
- Bytes: 374,587,922
- MD5: `4570b53c9d62ff6df281e622f0350060`
- License: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
- Downloader: `../download_reference.py`
- Source metadata: `replogle_figshare_metadata.json`
- Download/checksum provenance: `provenance.json`

`rejected_essential_reference.json` documents zero competition-target coverage in
the previously investigated essential-gene source. The downloaded `README.md`,
`dataset.toml`, `feature_names.json`, parquet snippets and parquet footer in this
local directory belong to that rejected exploration and are not used by the final
pipeline. They are excluded from Git. The large retained H5AD is also excluded
and can be recovered using the pinned downloader. No challenge ground truth was used.

Responses derived by this project are transformed estimates, not original raw
single-cell measurements. The reference supplies pseudobulk means; the final
submission simulates whole single-cell counts on the challenge control profiles.
