"""Regenerate every Paper-II figure as a vector PDF (alongside the existing PNGs).

matplotlib writes a true vector PDF when handed a ``*.pdf`` path, so this driver
monkeypatches ``Figure.savefig`` to emit a ``*.pdf`` next to every ``*.png`` it
writes, then runs the scripts that produce the five Paper-II figures as ``__main__``:

  foundations_fig1_reduction / _fig4_socsignature / _fig5_inferences  (make_foundations_figures)
  criticality_soc      (experiments/run_chr_soc.py)
  deposition_drive     (experiments/run_deposition_drive.py)

build.sh then points \\includegraphics at the .pdf, so xelatex embeds
the vector versions.

Run from anywhere:
    python paper/make_foundations_vector.py
"""
import os
import sys
import runpy

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# src/ (for `import sede`) + src/experiments/ & src/scripts/ (flat-namespace sibling imports;
# run_chr_soc / run_deposition_drive import sede and, via siblings, share the namespace).
for _p in (os.path.join(ROOT, "src", "scripts"), os.path.join(ROOT, "src", "experiments"), os.path.join(ROOT, "src")):
    sys.path.insert(0, _p)
os.chdir(ROOT)  # the run-scripts write to relative results/ paths

import matplotlib
matplotlib.use("Agg")
import matplotlib.figure as mfig

_orig_savefig = mfig.Figure.savefig


def _savefig(self, fname, *args, **kwargs):
    _orig_savefig(self, fname, *args, **kwargs)
    if isinstance(fname, (str, os.PathLike)):
        s = os.fspath(fname)
        if s.lower().endswith(".png"):
            pdf = s[:-4] + ".pdf"
            kw = {k: v for k, v in kwargs.items() if k != "dpi"}  # dpi irrelevant for vector
            _orig_savefig(self, pdf, *args, **kw)
            print(f"  + vector  {pdf}")


mfig.Figure.savefig = _savefig

SCRIPTS = [
    "figures/make_foundations_figures.py", # fig1_reduction, fig4_socsignature, fig5_inferences
    "src/experiments/run_chr_soc.py",      # criticality_soc (Fig 2)
    "src/experiments/run_deposition_drive.py", # deposition_drive (Fig 3)
]

for sc in SCRIPTS:
    print(f"=== {sc} ===")
    try:
        runpy.run_path(os.path.join(ROOT, sc), run_name="__main__")
    except SystemExit:
        pass
    except Exception as e:  # one bad script should not block the rest
        print(f"  {sc} FAILED: {e}")

print("\nDone. Vector PDFs written next to the PNGs in results/.")
