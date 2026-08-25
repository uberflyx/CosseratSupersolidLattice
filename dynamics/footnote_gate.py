"""Standing gate: every script the monograph footnotes must exist in this repository.

WHY THIS EXISTS.  A footnote pointing at a script that was never pushed is worse than no
footnote, because a reader following it finds nothing and concludes the result was never
checked.  Three such references accumulated in the corpus before anyone looked: one
written by a pass that produced the script only in scratch, and two citing scripts a
parallel pass had described but not committed.  None was caught by rendering, by the
cross-reference validator, or by the literal audit, because none of those tools looks
outside its own repository.

WHAT IT CHECKS.  Both directions, since both are defects:
  - dangling: cited in the corpus, absent here.  A reader hits a dead end.
  - orphaned: present here, cited nowhere.  Either the footnote was lost in an edit or
    the script has no reader, and the two cases need different fixes.

USAGE.  Point it at a checkout of the monograph:

    python footnote_gate.py /path/to/CosseratSupersolid

Exits non-zero if anything is dangling.  Orphans are reported but do not fail, since a
script may legitimately predate the text that will cite it.
"""
import os
import re
import subprocess
import sys


def cited_scripts(corpus_root):
    """Every dynamics/*.py path footnoted anywhere in the corpus, with its file."""
    found = {}
    for dirpath, _, files in os.walk(os.path.join(corpus_root, "quarto")):
        for fn in files:
            if not fn.endswith((".qmd", ".tex")):
                continue
            path = os.path.join(dirpath, fn)
            with open(path, errors="ignore") as fh:
                for m in re.finditer(r"(?:dynamics|neutrinos|gravity)/[\w]+\.py", fh.read()):
                    found.setdefault(m.group(0), set()).add(fn)
    return found


def committed_scripts(repo_root):
    """Every .py tracked by git here, so an uncommitted working file does not pass."""
    out = subprocess.run(["git", "-C", repo_root, "ls-files", "*.py"],
                         capture_output=True, text=True, check=True).stdout
    return {p for p in out.split() if p.endswith(".py")}


def main(corpus_root):
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    cited, have = cited_scripts(corpus_root), committed_scripts(here)
    dangling = {k: v for k, v in cited.items() if k not in have}
    orphaned = sorted(p for p in have
                      if p.startswith("dynamics/") and p not in cited
                      and not os.path.basename(p).startswith("_"))

    print("cited by the corpus: %d    tracked here: %d\n" % (len(cited), len(have)))
    for k in sorted(cited):
        if k not in dangling:
            print("  ok        %-42s cited in %s" % (k, ", ".join(sorted(cited[k]))))
    for k in sorted(dangling):
        print("  DANGLING  %-42s cited in %s" % (k, ", ".join(sorted(dangling[k]))))
    if orphaned:
        print("\n  orphaned (present, cited nowhere; not a failure):")
        for p in orphaned:
            print("      %s" % p)

    if dangling:
        print("\n%d dangling reference(s): the corpus points at scripts that are not here."
              % len(dangling))
        return 1
    print("\nno dangling references")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1]))
