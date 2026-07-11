from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

project = "modBHJet"
author = "modBHJet contributors"
copyright = "2026, modBHJet contributors"

extensions = [
    "myst_nb",
    "sphinx_copybutton",
    "sphinx_design",
    "breathe",
]

master_doc = "index"

html_theme = "pydata_sphinx_theme"
html_title = "modBHJet documentation"
html_static_path = ["_static"]
html_theme_options = {
    "github_url": "https://github.com/rubyduncan/modBHJet",
    "show_toc_level": 2,
}

myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "dollarmath",
]
myst_heading_anchors = 3

nb_execution_mode = "off"
nb_execution_timeout = 120
nb_merge_streams = True

breathe_projects = {
    "modbhjet": str(ROOT / "docs" / "_build" / "doxygen" / "xml"),
}
breathe_default_project = "modbhjet"

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
]

suppress_warnings = [
    "mystnb.unknown_mime_type",
]
