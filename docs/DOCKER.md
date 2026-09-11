# Container & Reproducibility — Docker (I90)

*Dockerfile:* `Dockerfile:1` `FROM python:3.11-slim`, `Dockerfile:2` `COPY . /app`, `Dockerfile:3` `WORKDIR /app`, `Dockerfile:4` `CMD python -m http.server 8000`.

## Build & Run

```bash
# build
docker build -t simulation_coarse .

# run — serves on http://127.0.0.1:8000/
docker run --rm -p 8000:8000 simulation_coarse
# open http://127.0.0.1:8000/ in browser
```

The image is `python:3.11-slim` only (no extra deps) — `COPY . /app` + `WORKDIR /app` + `CMD python -m http.server 8000` mirrors the `README.md` quick-start `python3 -m http.server 8123` but containerized for reproducibility. Any `docker run -p` mapping works (e.g. `-p 8123:8000` then `http://127.0.0.1:8123/`).

## Verify

```bash
docker run -d -p 8000:8000 --name coarse_test simulation_coarse
curl -s http://127.0.0.1:8000/ | head -20   # should return index.html
docker stop coarse_test
```

## Reproducibility notes

- No build step — app is static ES modules (`index.html` + `src/*.js?v=10`), so the container just serves files.
- Pin `python:3.11-slim` (Debian slim) for a stable base; `python -m http.server` is stdlib, no pip.
- For exact parity with host run, `python -m http.server 8000` serves from `/app` (`WORKDIR /app`) exactly as `python3 -m http.server 8123` from the checkout root.
```

