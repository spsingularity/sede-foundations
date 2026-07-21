#!/usr/bin/env python3
"""
T2.7 — Reference verification against live arXiv  [REVISION_PLAN.md]
====================================================================
CL: AI-assisted lit search hallucinates plausible-but-fake citations; a single fake
reference costs disproportionate credibility. "Verified against INSPIRE-HEP 2026-06"
is an assertion, not a check. This queries the LIVE arXiv API for every arXiv ID cited
in paper/SEDE.md, reports the real title (so it can be eyeballed against the claimed
content), and FLAGS any ID that does not resolve — with post-2024 IDs (the risk set)
called out separately.

Run: python verify_refs.py    (needs network)
"""
import re, sys, time, urllib.request, urllib.parse

PAPER = 'paper/SEDE.md'
API = "https://export.arxiv.org/api/query?id_list=%s&max_results=%d"


def extract_ids(path):
    txt = open(path).read()
    ids = sorted(set(re.findall(r"arXiv:([0-9]{4}\.[0-9]{4,5})", txt)))
    # one-line context for each (first occurrence)
    ctx = {}
    for i in ids:
        m = re.search(r"([^.\n]*arXiv:%s[^.\n]*)" % re.escape(i), txt)
        ctx[i] = (m.group(1).strip()[:130] if m else "")
    return ids, ctx


def query(ids):
    """Return {id: title or None} for a batch of ids."""
    url = API % (",".join(ids), len(ids) * 2 + 5)
    req = urllib.request.Request(url, headers={'User-Agent': 'SEDE-refcheck/1.0'})
    xml = urllib.request.urlopen(req, timeout=30).read().decode()
    # parse <entry> blocks: each has <id>...abs/XXXX</id> and <title>...</title>
    out = {}
    for entry in re.findall(r"<entry>(.*?)</entry>", xml, re.S):
        idm = re.search(r"abs/([0-9]{4}\.[0-9]{4,5})", entry)
        tm = re.search(r"<title>(.*?)</title>", entry, re.S)
        if idm:
            out[idm.group(1)] = re.sub(r"\s+", " ", tm.group(1)).strip() if tm else "(no title)"
    return out


def main():
    ids, ctx = extract_ids(PAPER)
    print(f"=== {len(ids)} arXiv IDs in {PAPER}; querying live arXiv ===\n")
    found = {}
    for k in range(0, len(ids), 15):
        batch = ids[k:k + 15]
        try:
            found.update(query(batch))
        except Exception as e:
            print(f"  [batch {k} query error: {e}]")
        time.sleep(3)   # arXiv politeness
    missing, post2024 = [], []
    for i in ids:
        t = found.get(i)
        post = i[:2] in ('25', '26')
        tag = "POST-2024" if post else ""
        if t is None:
            missing.append(i)
            print(f"  ✗ MISSING  {i}  {tag}\n              cited as: {ctx[i]}")
        else:
            mark = "⚠" if post else "·"
            print(f"  {mark} {i}  {t[:72]}")
            if post:
                post2024.append(i)
    print("\n" + "=" * 70)
    print(f"  resolved: {len(found)}/{len(ids)}   missing/unresolved: {len(missing)}")
    if missing:
        print(f"  ⚠ UNRESOLVED (verify by hand — possible hallucination): {', '.join(missing)}")
    print(f"  post-2024 IDs (eyeball title vs claimed content): {len(post2024)} flagged with ⚠ above")
    print("  NOTE: a resolved ID with a title that does NOT match the cited claim is still a problem —")
    print("        this checks existence + title; the title↔claim match needs a human read.")


if __name__ == '__main__':
    main()
