# Documentation and GitHub Pages

The committed `docs/` directory is a self-contained static site. It uses relative links, local CSS and optional JavaScript, and system fonts. Reading the guides needs no renderer, server-side application, CDN, network connection, or JavaScript.

## Preview locally

From the repository root:

```sh
python -m http.server 8000 --directory docs
```

Open `http://localhost:8000/` in your browser. Opening `docs/index.html` directly also works; the Copy button may be unavailable under file URLs, where code remains selectable. Narrow screens provide horizontally scrollable navigation, tables, and code blocks. Keyboard users can tab to the skip link and scrollable content.

## Configure GitHub Pages

For [pycoloc](https://github.com/LucaJiang/pycoloc), open **Settings → Pages → Build and deployment**.
Select **Source: Deploy from a branch**, **Branch: main**, and **Folder: /docs**,
then click **Save**. The committed `docs/.nojekyll` serves the prebuilt static files.
The expected address after a successful deployment is `https://LucaJiang.github.io/pycoloc/`.

Repository upload and Pages deployment are separate steps. This source delivery
prepares the site but does not claim an enabled Pages setting or a successful live
build. Check the repository's Pages status after saving the setting. SSH Git access
does not grant access to the Pages administration API.

## Regenerate HTML

Markdown files in `docs/` are the maintainable sources. The renderer version is pinned separately from runtime dependencies:

```sh
python -m pip install -r tools/docs-requirements.txt
python tools/build_docs.py
python tools/build_docs.py --check
python tools/check_docs.py
```

The build derives package version from `pyproject.toml`, applies navigation from `tools/site_config.json`, and renders every Markdown guide. It emits no timestamps or remote assets. `--check` fails if a source change has not been rendered. Edit `docs/assets/site.css` for theme changes; optional code-copy behavior is in `docs/assets/site.js`.

The offline link checker checks local page, anchor, and asset targets and package metadata. It does not claim that external websites are reachable. Browser screenshots and visual inspection are recorded separately in the release review; a link check alone is not a visual test.
