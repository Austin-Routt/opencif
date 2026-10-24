"""docs_pdf.py: render the reader-facing documents to PDF, for reading away from a Markdown viewer.

The documents are OPTIONS, EVALUATION, ALIGNMENT, LIMITATIONS, DISSERTATION_NOTES, and REPRODUCTIONS. The Markdown
files are the source, and the PDFs are regenerated rather than edited.

Usage: python3 tools/docs_pdf.py   (needs the `markdown` package and Playwright's Chromium)"""
import os
import sys

import markdown

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = ['OPTIONS', 'EVALUATION', 'ALIGNMENT', 'LIMITATIONS', 'DISSERTATION_NOTES', 'REPRODUCTIONS']
CSS = """
body { font-family: 'DejaVu Sans', Arial, sans-serif; font-size: 8.6pt; line-height: 1.35; color: #111; }
h1 { font-size: 15pt; margin: 0 0 6pt 0; } h2 { font-size: 12pt; margin: 14pt 0 4pt 0; border-bottom: 1px solid #bbb; }
h3 { font-size: 10pt; margin: 10pt 0 3pt 0; }
table { border-collapse: collapse; width: 100%; margin: 4pt 0 8pt 0; page-break-inside: auto; font-size: 8pt; }
tr { page-break-inside: avoid; } th, td { border: 1px solid #c8c8c8; padding: 2pt 4pt; vertical-align: top; text-align: left; }
th { background: #eef1f5; } code { font-family: 'DejaVu Sans Mono', monospace; font-size: 7.8pt; background: #f4f4f4; }
pre { background: #f6f6f6; padding: 5pt; font-size: 7.6pt; white-space: pre-wrap; }
"""


def main():
    """Write docs/<name>.pdf for each document on landscape A4 pages, with the document's name and the page number in
    the footer."""
    from playwright.sync_api import sync_playwright
    out = []
    with sync_playwright() as p:
        b = p.chromium.launch()
        pg = b.new_page()
        for name in DOCS:
            src = os.path.join(ROOT, 'docs', name + '.md')
            html = markdown.markdown(open(src).read(), extensions=['tables', 'fenced_code'])
            pg.set_content('<html><head><meta charset="utf-8"><style>%s</style></head><body>%s</body></html>' % (CSS, html))
            pdf = os.path.join(ROOT, 'docs', name + '.pdf')
            pg.pdf(path=pdf, format='A4', landscape=True, print_background=True,
                   margin=dict(top='12mm', bottom='12mm', left='10mm', right='10mm'),
                   display_header_footer=True, header_template='<span></span>',
                   footer_template='<div style="font-size:7pt;width:100%%;text-align:center;color:#666">opencif %s '
                                   '&middot; <span class="pageNumber"></span>/<span class="totalPages"></span></div>' % name)
            out.append(pdf)
        b.close()
    for f in out:
        print('wrote %s' % f)


if __name__ == '__main__':
    sys.exit(main())
