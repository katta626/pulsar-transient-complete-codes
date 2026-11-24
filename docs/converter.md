# converter — PSRFITS search mode to SIGPROC filterbank

Converts one or more PSRFITS SEARCH-mode files into a single SIGPROC
filterbank (`.fil`) file. Requires `numpy` and `astropy`.

This was previously the standalone `psr-searchmode-to-filterbank`
repository; its history is preserved here.

## Usage

From the `filterbank/` directory:

```bash
python3 -m converter obs_0001.fits obs_0002.fits -o obs.fil
```

Or from anywhere, using the wrapper script:

```bash
filterbank/converter/fits2fil obs_0001.fits -o obs.fil --overwrite
```

From Python:

```python
from converter import convert_psrfits_to_filterbank
convert_psrfits_to_filterbank(["obs_0001.fits"], "obs.fil")
```
