"""docs_check.py: scan documents for the constructions the author's style sheet excludes, and compare the facts of two
versions of a document.

Usage: python -I docs_check.py scan FILE ...
       python -I docs_check.py facts OLD NEW

The scan skips fenced code blocks, inline code spans, and quoted text, and it lists each excluded construction with its
line. The facts comparison collects the numbers, code spans, identifiers of checks and results (M14, P8, C3b, AA6), and
quoted passages of each version and lists those whose count fell in the new version, so that a rewrite that dropped or
changed a fact is found.
"""
import collections, re, sys

EXCLUDED = [
    (r'\u2014|–(?![\d])', 'em or en dash in prose'),
    (r'(?<![-|])\s--\s(?![-|])', 'double-hyphen dash'),
    (r'\byou\b|\byour\b|\byours\b', '"you"'),
    (r'\*\*', 'bold'),
    (r'\b(notably|importantly|crucially|essentially|straightforward\w*|simply|clearly|robust\w*|ensures?|ensuring|'
     r'leverag\w*|comprehensive|overall|interestingly|specifically|namely|in order to|it is worth noting|'
     r'consequently|key takeaways?|overview)\b', 'excluded word'),
    (r"\b\w+n't\b|\b(it's|that's|there's|we're|we've|let's)\b", 'contraction'),
    (r'\bnot\b[^.;|]{1,50}\bbut\b', '"not X but Y"'),
    (r'(?<![\w\])])\?(\s|$)', 'question'),
    (r'[a-z]\)?: [a-z]', 'clause-chaining colon (check)'),
]
BRITISH = r'\b(centre[sd]?|centring|realis\w*|modell\w*|licence\w*|neighbour\w*|grey\w*|metres?|behaviour\w*|' \
          r'colour\w*|favour\w*|normalis\w*|minimis\w*|maximis\w*|organis\w*|recognis\w*|parameteris\w*|analys(e|ed|ing)|' \
          r'catalogue\w*|programme\w*|labelled|labelling|travell\w*|fibre\w*|optimis\w*|summaris\w*|generalis\w*|' \
          r'characteris\w*|equalis\w*|stabilis\w*|utilis\w*|emphasis(e|ed|es)|prioritis\w*|visualis\w*|initialis\w*|' \
          r'finalis\w*|digitis\w*|discretis\w*|factoris\w*|micrometres?|millimetres?|whilst|amongst|judgement)\b'


def prose_lines(text):
    fence = False
    for n, line in enumerate(text.splitlines(), 1):
        if line.lstrip().startswith('```'):
            fence = not fence
            continue
        if fence:
            continue
        t = re.sub(r'`[^`]*`', ' ', line)
        t = re.sub(r'"[^"\n]*"', ' ', t)
        t = re.sub(r'\$[^$]*\$', ' ', t)
        t = re.sub(r'https?://\S+', ' ', t)
        t = re.sub(r'\b\w*_\w*\b', ' ', t)
        yield n, line, t


def scan(path):
    text = open(path, encoding='utf-8').read()
    hits = []
    for n, line, t in prose_lines(text):
        for pat, why in EXCLUDED:
            if re.search(pat, t, 0 if why == 'bold' else re.I):
                hits.append((n, why, line.strip()[:120]))
        m = re.search(BRITISH, t, re.I)
        if m:
            hits.append((n, 'British: ' + m.group(0), line.strip()[:120]))
    return hits


def facts(text):
    c = collections.Counter()
    text = re.sub(r'```.*?```', lambda m: m.group(0).replace('`', "'"), text, flags=re.S)   # fenced blocks as plain text
    for m in re.finditer(r'`[^`\n]+`', text):
        c['code ' + m.group(0)] += 1
    for m in re.finditer(r'"[^"\n]{3,}"', text):
        c['quote ' + m.group(0)] += 1
    t = re.sub(r'(?<=\d)[ ,](?=\d{3}\b)', '', text)             # thousands separators
    for m in re.finditer(r'(?<![\w.])[-−]?\d+(?:\.\d+)?(?:\s?[×x]\s?10[⁻⁰¹²³⁴⁵⁶⁷⁸⁹-]+|e[-+]?\d+)?', t):
        c['number ' + m.group(0).replace('−', '-').replace(' ', '')] += 1
    for m in re.finditer(r'\b[A-Z]{1,3}\d+[a-z]?\b', text):
        c['id ' + m.group(0)] += 1
    return c


def main(argv):
    if argv[1] == 'scan':
        total = 0
        for p in argv[2:]:
            hits = scan(p)
            total += len(hits)
            print('%s: %d flags' % (p, len(hits)))
            for n, why, line in hits:
                print('  %5d  %-30s %s' % (n, why, line))
        return 0
    if argv[1] == 'facts':
        a = facts(open(argv[2], encoding='utf-8').read())
        b = facts(open(argv[3], encoding='utf-8').read())
        lost = [(k, a[k], b[k]) for k in sorted(a) if b[k] < a[k]]
        print('%d fact tokens in the old version, %d of them with a lower count in the new version' % (len(a), len(lost)))
        for k, x, y in lost:
            print('  %-70s %d -> %d' % (k[:70], x, y))
        return 0
    return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv))
