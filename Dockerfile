FROM python:3.11-slim
COPY . /app
WORKDIR /app
CMD python -m http.server 8000
