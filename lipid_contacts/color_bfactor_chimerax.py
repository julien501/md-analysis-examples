# ChimeraX script: colour B-factor structures from lipid_contacts.py white (no contact) to dark green (always)
# and save an image and a session next to each file.
#   chimerax --offscreen --script "color_bfactor_chimerax.py out/bfactor_all.pdb [more.pdb ...]"
import sys
from pathlib import Path

from chimerax.core.commands import run

for path in sys.argv[1:]:
    stem = Path(path).with_suffix("")
    model = run(session, f"open {path}")[0]          # noqa: F821  (session is provided by ChimeraX)
    m = f"#{model.id_string}"
    run(session, f"hide {m} atoms; show {m} cartoons")   # lipids, water and ions have no cartoon
    run(session, f"color bfactor {m} palette 0,white:100,#006400")
    run(session, "wait 1")                             # draw one frame so the cartoons exist before 'view'
    run(session, f"view {m}; lighting soft")
    run(session, f"save {stem}.png width 1600 supersample 3")
    run(session, f"save {stem}.cxs")
    run(session, f"close {m}")
