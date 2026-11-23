#!/usr/bin/env bash
# Render a Markdown document with LaTeX math to a typeset PDF:
#   Markdown -> HTML (pandoc, KaTeX) -> PDF (headless Chromium via Playwright).
# One-time setup in this directory: npm install && npx playwright install chromium
# Usage, from the repository root:
#   tools/render-methods/render.sh docs/methods/PaPeR-methods.md docs/methods/PaPeR-methods.pdf
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
src="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
out="$2"
html="$(dirname "$src")/.render.html"
pandoc "$src" --from markdown+tex_math_dollars-yaml_metadata_block --to html5 --standalone \
  --katex="file://$here/node_modules/katex/dist/" --css="file://$here/paper.css" \
  --metadata pagetitle="PaPeR methods" --resource-path="$(dirname "$src")" -o "$html"
node "$here/print.js" "$html" "$out"
rm -f "$html"
