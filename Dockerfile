FROM python:3.12-slim-bookworm@sha256:6a0e06f5ba0b373d051aabd3fffb9b873f27f9a2a3a14bcbe37dee6f0e3bbc70

LABEL org.opencontainers.image.title="HelixForge-SMansoni Atlas" \
      org.opencontainers.image.source="https://github.com/GMiguelAlves/HelixForge-SMansoni" \
      org.opencontainers.image.version="2026.10.21" \
      org.opencontainers.image.licenses="Apache-2.0"

ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
WORKDIR /app
COPY scripts/atlas/runtime.py /app/scripts/atlas/runtime.py
COPY results/atlas/rnaseq/ /app/results/atlas/rnaseq/
COPY LICENSE /app/LICENSE
USER 65532:65532
EXPOSE 8080
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
  CMD ["python", "-c", "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8080/healthz', timeout=3).read()"]
ENTRYPOINT ["python", "/app/scripts/atlas/runtime.py", "serve", "--directory", "/app/results/atlas/rnaseq", "--bind", "0.0.0.0", "--port", "8080", "--no-browser"]
