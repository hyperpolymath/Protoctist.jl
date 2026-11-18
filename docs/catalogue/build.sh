#!/usr/bin/env bash
# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Build the field catalogue: review sheet (CSV), HTML, and PDF if a Chromium
# binary is given.
#
#   docs/catalogue/build.sh                 # review-sheet.csv + HTML in $OUT
#   CHROMIUM=/path/to/chromium docs/catalogue/build.sh   # also the PDF
#
# Needs asciidoctor. The HTML loads MathJax from cdnjs to typeset formulae.

set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
out="${OUT:-$here/out}"
mkdir -p "$out"

# One row per catalogue entry, from the same source the document is built from.
{
    printf 'id,title,decision,need (now/later/never),decision (agree/differ),correct (yes/no/unsure),comment\n'
    awk '
        function flush() {
            if (id != "") { gsub(/"/, "\"\"", title); printf "%s,\"%s\",%s,,,,\n", id, title, dec }
        }
        /^=== PC-[A-N][0-9][0-9] / {
            flush(); id = $2; title = $0; sub(/^=== PC-[A-N][0-9][0-9] · /, "", title); dec = ""
        }
        /^Decision::/ && id != "" && dec == "" {
            line = $0; d = ""
            while (match(line, /#[A-Z]+#/)) {
                d = d (d == "" ? "" : "+") substr(line, RSTART + 1, RLENGTH - 2)
                line = substr(line, RSTART + RLENGTH)
            }
            dec = d
        }
        END { flush() }
    ' "$here"/parts/0[4-6]-*.adoc
} > "$here/review-sheet.csv"

entries=$(($(wc -l < "$here/review-sheet.csv") - 1))

asciidoctor -a entry-count="$entries" -a stylesheet=catalogue.css -a stylesdir="$here" \
    -a linkcss! -a webfonts! --failure-level WARN \
    -o "$out/protoctist-field-catalogue.html" "$here/protoctist-field-catalogue.adoc"
mkdir -p "$out/img" && cp "$here/img/cover-motif.svg" "$out/img/"
echo "entries: $entries"
echo "html:    $out/protoctist-field-catalogue.html"

if [[ -n "${CHROMIUM:-}" ]]; then
    "$CHROMIUM" --headless --no-sandbox --disable-gpu --no-pdf-header-footer \
        --virtual-time-budget=60000 --run-all-compositor-stages-before-draw \
        --print-to-pdf="$out/protoctist-field-catalogue.pdf" \
        "file://$out/protoctist-field-catalogue.html" 2>/dev/null
    echo "pdf:     $out/protoctist-field-catalogue.pdf"
fi
