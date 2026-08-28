# Vendored web runtime

`three.r149.min.js` is the official Three.js r149 browser build.

- Source: https://raw.githubusercontent.com/mrdoob/three.js/r149/build/three.min.js
- Upstream: https://github.com/mrdoob/three.js/tree/r149
- License: MIT (the upstream license notice is retained in the minified file)
- SHA-256: `8A5F7249903B54D30F79F708699D2FED2D6A1D0741A4CD41377D1F01BB5A2271`

`create_lmm_cinematic.py` embeds this file into `docs/lmm_cinematic.html` so the
briefing remains a single offline HTML file with no runtime network dependency.
