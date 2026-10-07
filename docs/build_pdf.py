"""Build docs/bouake-transcriptomics.pdf from the manuscript markdown.

The manuscript is a Google Docs export (backslash-escaped punctuation, reference-style data-URI
images), so this script cleans the escapes, converts to HTML with python-markdown, applies a
manuscript style sheet and prints to PDF with headless Chrome.

    python docs/build_pdf.py        # needs the `markdown` package and Google Chrome

Writes docs/bouake-transcriptomics.pdf (the intermediate HTML is removed).
"""

import re
import subprocess
from pathlib import Path

import markdown

DOCS = Path(__file__).resolve().parent
SOURCE = DOCS / "bouake-transcriptomics.md"
HTML = DOCS / "bouake-transcriptomics.html"
PDF = DOCS / "bouake-transcriptomics.pdf"
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

STYLE = """
@page { size: A4; margin: 20mm 18mm 22mm 18mm;
        @bottom-center { content: counter(page); font: 9pt 'Helvetica Neue', Arial, sans-serif; color: #666; } }
html { font-size: 10.5pt; }
body { font-family: 'Charter', 'Georgia', 'Times New Roman', serif; line-height: 1.42; color: #111; }
p { margin: 0 0 0.65em 0; text-align: justify; hyphens: auto; }
body > p:first-of-type { font-size: 1.45em; line-height: 1.25; text-align: left; margin-bottom: 0.5em; }
h1 { font-family: 'Helvetica Neue', Arial, sans-serif; font-size: 1.3em; margin: 1.4em 0 0.5em 0;
     border-bottom: 1px solid #ccc; padding-bottom: 0.2em; page-break-after: avoid; }
h1 strong { font-weight: 600; }
a { color: inherit; text-decoration: none; }
img { max-width: 100%; max-height: 150mm; width: auto; height: auto; display: block; margin: 0 auto; }
table { border-collapse: collapse; width: 100%; margin: 0.6em 0 1.2em 0; font-size: 0.88em; }
th, td { padding: 3px 6px; text-align: left; vertical-align: top; }
/* figure blocks: one-column tables holding an image and its legend */
table:has(img) { page-break-inside: avoid; border: none; margin: 1em 0 1.4em 0; }
table:has(img) th, table:has(img) td { border: none; padding: 2px 0; font-size: 0.9em; line-height: 1.35; }
table:has(img) th { padding-bottom: 8px; }
table:has(img) td { text-align: justify; hyphens: auto; }
table:not(:has(img)) th, table:not(:has(img)) td { border-bottom: 1px solid #ddd; }
table:not(:has(img)) th { border-bottom: 1.5px solid #666; font-family: 'Helvetica Neue', Arial, sans-serif; }
table:not(:has(img)) { page-break-inside: avoid; }
sup { font-size: 0.75em; line-height: 0; }
p:has(+ table) { break-after: avoid; }
"""


def clean_escapes(text: str) -> str:
    """Remove the backslashes Google Docs puts before punctuation that markdown does not need."""
    text = text.replace("\\<", "&lt;").replace("\\>", "&gt;")
    text = text.replace("\n1Department", "\n<sup>1</sup>Department")
    return re.sub(r"\\([=~|:;,?@$%&^\"'/])", r"\1", text)


def build() -> None:
    text = clean_escapes(SOURCE.read_text())
    body = markdown.markdown(text, extensions=["tables", "sane_lists"])
    HTML.write_text(f"<!doctype html><html><head><meta charset='utf-8'>"
                    f"<title>Molecular mechanisms of resistance to pirimiphos-methyl</title>"
                    f"<style>{STYLE}</style></head><body>{body}</body></html>")
    subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--no-sandbox", "--no-pdf-header-footer",
                    f"--print-to-pdf={PDF}", f"file://{HTML}"], check=True, capture_output=True)
    HTML.unlink()
    print(f"wrote {PDF} ({PDF.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    build()
