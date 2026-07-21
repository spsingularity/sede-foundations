#!/usr/bin/env bash
# Build the SciPost Physics PDF (SciPost.cls) for Paper I from SEDE_foundations.md.
#   makedoc (title/abstract/keywords -> YAML, mdmath wraps informal notation, fold figures)
#   -> protect $...$ as raw-LaTeX -> pandoc --natbib -> pdflatex + bibtex (SciPost_bibstyle).
set -e
cd "$(dirname "$0")"
mkdir -p tex
BASE=SEDE_foundations
BIB=refs_foundations

python3 tools/makedoc.py $BASE.md .build.md
trap 'rm -f .build.md' EXIT
python3 - <<'PY'
import re
t=open('.build.md',encoding='utf-8').read()
t=re.sub(r'\n##\s+References\s*\n+:::\s*\{#refs\}\s*\n:::\s*\n','\n',t)
t=re.sub(r'(?<!\\)(?<!\$)\$(?!\$)((?:\\.|[^$\\\n]|\\\n)+?)\$(?!\$)',
         lambda m: '`\\('+m.group(1)+'\\)`{=latex}', t)
open('.build.md','w',encoding='utf-8').write(t)
PY

pandoc -f markdown-superscript-subscript .build.md -o tex/$BASE.tex \
  --standalone --shift-heading-level-by=-1 --natbib \
  --template=tools/template_scipost.tex

# SciPost uses the cite package (numeric \cite), not natbib's \citep/\citet
perl -0pi -e 's#\\cite[tp]\{#\\cite{#g' tex/$BASE.tex
# point at this paper's bibliography file
perl -0pi -e "s#\\\\bibliography\\{refs\\}#\\\\bibliography{$BIB}#" tex/$BASE.tex
# number display equations
perl -0pi -e 's/\\\[/\\begin{equation}/g; s/\\\]/\\end{equation}/g' tex/$BASE.tex

( cd tex && \
  pdflatex -interaction=nonstopmode $BASE.tex >$BASE.build.log 2>&1 ; \
  BIBINPUTS="..:$BIBINPUTS" bibtex $BASE      >>$BASE.build.log 2>&1 ; \
  pdflatex -interaction=nonstopmode $BASE.tex >>$BASE.build.log 2>&1 ; \
  pdflatex -interaction=nonstopmode $BASE.tex >>$BASE.build.log 2>&1 ) || true

if [ -f tex/$BASE.pdf ]; then
  cp tex/$BASE.pdf $BASE.pdf
  echo "built paper/$BASE.pdf"
  grep -c "^!" tex/$BASE.build.log | awk '{print $1" LaTeX errors (see tex/'$BASE'.build.log)"}'
  grep -c "Warning--" tex/$BASE.build.log 2>/dev/null | awk '{print $1" bibtex warnings"}'
else
  echo "BUILD FAILED — see tex/$BASE.build.log"; grep -A2 '^!' tex/$BASE.build.log | head -20
fi
