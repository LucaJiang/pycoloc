#!/usr/bin/env python3
"""Check local Markdown/Pages links, HTML anchors and assets offline."""
from __future__ import annotations

import argparse
import json
import re
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.ids = set()
        self.links = []
        self.errors = []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if "id" in a:
            if a['id'] in self.ids:
                self.errors.append(f"duplicate id {a['id']}")
            self.ids.add(a['id'])
        if tag == 'img' and not a.get('alt'):
            self.errors.append('image lacks nonempty alternative text')
        for key in ('href', 'src'):
            if key in a:
                self.links.append(a[key])


def check():
    errors, checked = [], 0
    docs = ROOT / 'docs'
    pages = {p.resolve(): Page(p.read_text()) for p in docs.rglob('*.html')}
    for path, page in pages.items():
        errors.extend(f'{path.relative_to(ROOT)}: {e}' for e in page.errors)
        for target in page.links:
            url = urlsplit(target)
            if url.scheme or url.netloc:
                if target.startswith(('http:', 'https:', '//')) and re.search(r'(?:src|href)="' + re.escape(target) + '"', path.read_text()):
                    # External hyperlinks are fine; runtime assets must stay local.
                    if target in re.findall(r'<(?:script|link|img)[^>]+(?:src|href)="([^"]+)"', path.read_text()):
                        errors.append(f'{path.relative_to(ROOT)}: external runtime asset {target}')
                continue
            checked += 1
            dest = (path.parent / unquote(url.path)).resolve() if url.path else path
            if target.startswith('/') or not dest.is_relative_to(docs.resolve()):
                errors.append(f'{path.relative_to(ROOT)}: link escapes Pages root: {target}')
            elif not dest.is_file():
                errors.append(f'{path.relative_to(ROOT)}: missing {target}')
            elif url.fragment and dest.suffix == '.html' and unquote(url.fragment) not in pages[dest].ids:
                errors.append(f'{path.relative_to(ROOT)}: missing anchor {target}')
    for path in ROOT.rglob('*.md'):
        if any(x in path.parts for x in ('.git', '.venv', 'build', 'dist', 'target', '*.egg-info')):
            continue
        text = path.read_text()
        for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)', text):
            target = target.split('#', 1)[0].split(' ', 1)[0].strip('<>')
            if not target or urlsplit(target).scheme:
                continue
            checked += 1
            if not (path.parent / unquote(target)).exists():
                errors.append(f'{path.relative_to(ROOT)}: missing {target}')
        if '/home/' in text or '/opt/anaconda' in text:
            errors.append(f'{path.relative_to(ROOT)}: private path')
    return dict(passed=not errors, checked_links=checked, html_pages=len(pages), errors=errors,
                external_urls_checked=False, browser_visual_check=False)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, help='Optional JSON check record')
    args = p.parse_args()
    record = check()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))
    raise SystemExit(0 if record['passed'] else 1)
