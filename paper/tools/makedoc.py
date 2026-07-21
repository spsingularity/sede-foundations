#!/usr/bin/env python3
"""
Assemble the build-time markdown: (1) convert informal inline math via mdmath, then
(2) lift the H1 title + the Abstract section into YAML metadata and add the author, so the
PDF renders a proper title block (title / author / date / abstract) and --toc lands AFTER
the abstract — the JCAP/PRD convention. SEDE_cosmology.md itself is left untouched.
Run: python makedoc.py SEDE_cosmology.md OUT.md
"""
import re, sys, datetime, mdmath

def main(inp, outp):
    raw = open(inp, encoding="utf-8").read()
    # 1. wrap inline math
    body, spans = mdmath.protect(raw)
    text = mdmath.restore(mdmath.RUN.sub(mdmath.convert_run, body), spans)
    lines = text.split("\n")

    # 2a. title = first '# ' heading
    title = ""
    ti = next(i for i, L in enumerate(lines) if L.startswith("# "))
    title = lines[ti][2:].strip()

    # 2b. abstract = body of the '## Abstract' section, up to the next '## ' or a '---' rule
    ai = next(i for i, L in enumerate(lines) if re.match(r"##\s+Abstract\s*$", L))
    j = ai + 1
    abs = []
    while j < len(lines) and not (lines[j].startswith("## ") or lines[j].strip() == "---"):
        abs.append(lines[j]); j += 1
    # lift the "**Keywords —** …" line out of the abstract into its own field (journals want
    # keywords in a dedicated macro, not in the abstract body).
    keywords, abs2, in_kw = "", [], False
    for L in abs:
        m = re.match(r"\*\*Keywords\s*[—–-]+\*\*\s*(.+)", L.strip())
        if m:
            keywords = m.group(1).strip(); in_kw = True
        elif in_kw and L.strip():          # keyword list wrapped onto continuation line(s)
            keywords += " " + L.strip()
        else:
            in_kw = False
            abs2.append(L)
    keywords = keywords.rstrip(".")
    abstract = "\n".join(abs2).strip()

    # 2c. body = everything from the first numbered section onward
    bi = next(i for i, L in enumerate(lines) if re.match(r"##\s+1\.", L))
    rest = "\n".join(lines[bi:])

    # 2d. fold "![alt](path)\n\n**Figure N.** caption…" into a real pandoc captioned figure
    # ![caption](path){#fig:N} -> \begin{figure}\includegraphics\caption{…}\end{figure}.
    # The source keeps the GitHub-visible image + bold caption; only the build sees the float.
    def fold_figures(t):
        pat = re.compile(
            r"!\[[^\]]*\]\(([^)]+)\)[ \t]*\n[ \t\n]*"     # image (capture path), then any gap
            r"\*\*Figure\s+([0-9A-Za-z]+)\.\*\*\s*(.+?)"  # **Figure N.** caption …
            r"(?=\n[ \t]*\n|\Z)", re.S)
        def esc(cap):
            # escape literal [ ] OUTSIDE $...$ math — a bare ] would close pandoc's ![…] alt
            parts = re.split(r"(\$[^$]*\$)", cap)
            for i in range(0, len(parts), 2):            # even indices = non-math text
                parts[i] = parts[i].replace("[", r"\[").replace("]", r"\]")
            return "".join(parts)
        def repl(m):
            path, num = m.group(1).strip(), m.group(2)
            cap = esc(" ".join(m.group(3).split()))      # join multi-line caption, escape brackets
            return f"![{cap}]({path}){{#fig:{num}}}"
        return pat.sub(repl, t)
    rest = fold_figures(rest)

    # 2e. real section numbering + cross-refs. Strip the baked-in "N." from headings (LaTeX
    # auto-numbers), attach a stable \label id, inject \appendix, and mark back-matter
    # unnumbered — then rewrite prose "§N", "§N.M", "Section N", "App. X" into \ref{}.
    labels = set()
    def number_sections(t):
        out, appendix = [], False
        for L in t.split("\n"):
            if re.match(r"^## References\s*$", L):
                continue   # natbib's \bibliography prints its own "References" heading — drop ours
            m = re.match(r"^## (Code and data availability|Data availability statement|"
                         r"Use of AI tools|Funding and competing interests|Reproducibility|"
                         r"Acknowledgements|Acknowledgments)\s*$", L)
            if m:
                out.append(f"## {m.group(1)} {{.unnumbered}}"); continue
            m = re.match(r"^## Appendix ([A-G]) [—–-]+ (.+)$", L)   # ## Appendix X — Title
            if m:
                if not appendix:
                    out += ["", "```{=latex}", "\\appendix", "```", ""]; appendix = True
                labels.add(f"sec-{m.group(1)}")
                out.append(f"## {m.group(2)} {{#sec-{m.group(1)}}}"); continue
            m = re.match(r"^### ([A-G])\.([0-9A-Za-z]+) (.+)$", L)   # ### X.M  (appendix sub)
            if m:
                labels.add(f"sec-{m.group(1)}-{m.group(2)}")
                out.append(f"### {m.group(3)} {{#sec-{m.group(1)}-{m.group(2)}}}"); continue
            m = re.match(r"^## ([0-9]+)\. (.+)$", L)                 # ## N. Title
            if m:
                labels.add(f"sec-{m.group(1)}")
                out.append(f"## {m.group(2)} {{#sec-{m.group(1)}}}"); continue
            m = re.match(r"^### ([0-9]+)\.([0-9]+) (.+)$", L)        # ### N.M Title
            if m:
                labels.add(f"sec-{m.group(1)}-{m.group(2)}")
                out.append(f"### {m.group(3)} {{#sec-{m.group(1)}-{m.group(2)}}}"); continue
            out.append(L)
        return "\n".join(out)
    # convert a prose ref to \ref ONLY if the target label exists here (cross-paper "Paper I §5.6"
    # refs point at the companion and have no local label -> left as literal text).
    def convert_refs(t):
        # mask "§" directly following a companion-paper cue ("Paper I §6", "companion paper, §6",
        # "the foundations paper's §5", "cosmology paper §3.4") so those cross-paper refs are left
        # as literal text, not \ref'd to a same-numbered LOCAL section.
        CUE = (r"(?:Paper\s+I{1,3}V?|companion(?:\s+paper|\s+note)?(?:['’]s)?|"
               r"(?:foundations|cosmology|count|observational-tests|data|umbrella|framework|"
               r"scale-sector|supplement)\s+paper(?:['’]s)?)")
        t = re.sub(CUE + r"([,;:)\s]{0,3})§", lambda m: m.group(0)[:-1] + "\x01", t)
        # also mask the postfix form: "§5 of the foundations paper"
        t = re.sub(r"§(?=\s*[0-9A-G][0-9A-Za-z.]*\s+of\s+the\s+[\w-]+\s+paper)", "\x01", t)
        def sec(pfx, lab, whole):   # emit \ref if label is real here, else leave untouched
            return f"{pfx}\\ref{{{lab}}}" if lab in labels else whole
        t = re.sub(r"§§\s*([0-9]+)\s*[–-]\s*([0-9]+)",
                   lambda m: f"§§\\ref{{sec-{m[1]}}}--\\ref{{sec-{m[2]}}}"
                   if f"sec-{m[1]}" in labels and f"sec-{m[2]}" in labels else m[0], t)
        t = re.sub(r"§\s*([0-9]+)\.([0-9]+)", lambda m: sec("§", f"sec-{m[1]}-{m[2]}", m[0]), t)
        t = re.sub(r"§\s*([A-G])\.([0-9A-Za-z]+)", lambda m: sec("§", f"sec-{m[1]}-{m[2]}", m[0]), t)
        t = re.sub(r"§\s*([0-9]+)(?!\.[0-9])\b", lambda m: sec("§", f"sec-{m[1]}", m[0]), t)  # §N, not §N.M
        t = re.sub(r"§\s*([A-G])(?!\.[0-9A-Za-z])\b", lambda m: sec("§", f"sec-{m[1]}", m[0]), t)
        # NB use U+00A0 (real nbsp) not "~": pandoc escapes a literal ~ to \textasciitilde
        t = re.sub(r"\bSection\s+([0-9]+)\b", lambda m: sec("Section ", f"sec-{m[1]}", m[0]), t)
        t = re.sub(r"\bApp(?:endix|\.)\s+([A-G])\b", lambda m: sec("Appendix ", f"sec-{m[1]}", m[0]), t)
        return t.replace("\x01", "§")   # restore masked cross-paper §
    rest = convert_refs(number_sections(rest))

    def yml(s):  # escape for a double-quoted YAML scalar
        return s.replace("\\", "\\\\").replace('"', '\\"')

    date = datetime.date.today().strftime("%-d %B %Y")
    abs_block = "\n".join(("  " + L) if L.strip() else "" for L in abstract.split("\n"))
    front = (
        "---\n"
        f'title: "{yml(title)}"\n'
        'author: "Stilian Pandev ([ORCID 0009-0005-8153-071X](https://orcid.org/0009-0005-8153-071X))"\n'
        f'date: "{yml(date)}"\n'
        + (f'keywords: "{yml(keywords)}"\n' if keywords else "")
        + "abstract: |\n"
        f"{abs_block}\n"
        "---\n\n"
    )
    open(outp, "w", encoding="utf-8").write(front + rest + "\n")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
