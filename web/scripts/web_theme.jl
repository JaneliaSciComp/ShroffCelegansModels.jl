using Bonito: Asset

# Shared theme assets for the Shroff Bonito web apps.
#
# Bonito's `Asset` type, when included anywhere in a returned page DOM,
# is auto-hoisted into the document `<head>` as a `<link rel=stylesheet>`
# (see Bonito's render_asset / jsrender on Asset). The Asset itself
# renders as nothing in the body, so we don't introduce a wrapper div.
#
# Use in an App callback:
#
#     return DOM.main(theme_assets()...,
#         DOM.h2("..."),
#         ...
#     )

# pico.css v2 — classless framework (~10 KB, dark/light auto via prefers-color-scheme).
# Swap the URL to change frameworks; no per-page edit needed.
const PICO_CSS = Asset("https://cdn.jsdelivr.net/npm/@picocss/pico@2/css/pico.min.css")

# Local override file. Bonito hosts it via its own asset server, so the
# served URL is independent of nginx's static routes.
const SHROFF_CSS = Asset(joinpath(@__DIR__, "..", "static", "style.css"))

theme_assets() = (PICO_CSS, SHROFF_CSS)
