"""Convert report_template.html into ../README.md (GitHub Markdown).

Figures become PNG images from ../plots, Table 1 is built from the results, the report's own
"Files" section is replaced by readme_code_section.md. Run after plot_sim_results.R.
"""
import csv, html, pathlib, re
from html.parser import HTMLParser

here = pathlib.Path(__file__).resolve().parent
sim = here.parent
VOID = {"br", "img", "meta", "link", "hr", "input"}


class Node:
    def __init__(self, tag, attrs, parent=None):
        self.tag, self.attrs, self.parent, self.kids = tag, dict(attrs), parent, []


class Tree(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.root = Node("root", []); self.cur = self.root

    def handle_starttag(self, tag, attrs):
        n = Node(tag, attrs, self.cur); self.cur.kids.append(n)
        if tag not in VOID:
            self.cur = n

    def handle_endtag(self, tag):
        n = self.cur
        while n is not self.root and n.tag != tag:
            n = n.parent
        if n is not self.root:
            self.cur = n.parent

    def handle_data(self, data):
        self.cur.kids.append(data)


def inline(n):
    if isinstance(n, str):
        return re.sub(r"\s+", " ", n)
    t, k = n.tag, "".join(inline(c) for c in n.kids)
    if t == "b" or t == "strong":
        return f"**{k.strip()}**"
    if t == "i" or t == "em":
        return "" if "style" in n.attrs else f"*{k.strip()}*"
    if t == "code":
        return f"`{k}`"
    if t == "sub":
        return f"<sub>{k}</sub>"
    if t == "br":
        return "  \n"
    if t == "a":
        return f"[{k}]({n.attrs.get('href', '')})"
    return k


def cls(n):
    return n.attrs.get("class", "") if not isinstance(n, str) else ""


def table_md(rows):
    out = ["| " + " | ".join(rows[0]) + " |", "|" + "---|" * len(rows[0])]
    out += ["| " + " | ".join(r) + " |" for r in rows[1:]]
    return "\n".join(out)


def calib_table():
    rows = list(csv.DictReader(open(sim / "sim_dubois_calibrated_results.tsv"), delimiter="\t"))
    labels = {"diet_naive": "Diet · RDA (Dubois et al.)", "diet_means": "Diet · individual means",
              "dietGeo_naive": "Diet \\| geography · RDA (Dubois et al.)",
              "dietGeo_block": "Diet \\| geography · RDA (ID permutation)",
              "dietGeo_means": "Diet \\| geography · individual means",
              "nunavut_naive": "Nunavut only · RDA (Dubois et al.)",
              "nunavut_block": "Nunavut only · RDA (ID permutation)"}
    cols = []
    for d in ("euclid", "unifrac"):
        for r in rows:
            key = (d, r["rho_d"], r["R2true"])
            if r["dist"] == d and key not in cols:
                cols.append(key)
    head = ["Analysis"] + [f"{'E' if d == 'euclid' else 'U'} {float(rd):g}" + (f" (R²={float(t):g})" if float(t) > 0 else "")
                           for d, rd, t in cols]
    body = []
    for m, lab in labels.items():
        line = [lab]
        for d, rd, t in cols:
            r = next(x for x in rows if (x["dist"], x["rho_d"], x["R2true"]) == (d, rd, t) and x["method"] == m)
            v = f"{round(float(r['reject']) * 100)}%"
            if float(t) == 0 and not m.startswith("diet_") and float(r["reject"]) > 0.1:
                v = f"**{v}**"
            line.append(v)
        body.append(line)
    note = ("*Table 1.* Share of datasets with p < 0.05 in the study-matched simulation. Columns: distance "
            "(E = Euclidean, U = square-root UniFrac) and cor(diet | ID); R² is the true diet effect where it is not 0. "
            "Diet-alone rows test a real association, since diet and microbiome both differ by region. "
            "Bold: false-positive rate above 10% in a test whose null is true. "
            "Paper values: diet 0.11, diet | geography 0.064, Nunavut only 0.169 (all p ≤ 0.004).")
    return note + "\n\n" + table_md([head] + body)


out = []


def block(n):
    if isinstance(n, str):
        if n.strip():
            out.append(re.sub(r"\s+", " ", n).strip())
        return
    t, c = n.tag, cls(n)
    if t in ("style", "title", "meta", "link", "script"):
        return
    if t == "h1":
        out.append("# " + inline(n).strip())
    elif t == "h2":
        out.append("## " + inline(n).strip())
    elif t == "h3":
        out.append("### " + inline(n).strip())
    elif "eyebrow" in c:
        out.append("*" + inline(n).strip() + "*")
    elif t == "p":
        txt = inline(n).strip()
        if txt == "__CALIB_TABLE__":
            out.append(calib_table())
        elif txt:
            out.append(txt)
    elif t == "ul":
        out.append("\n".join("- " + inline(li).strip() for li in n.kids if not isinstance(li, str) and li.tag == "li"))
    elif "keyfind" in c:
        items = []
        for kf in n.kids:
            if isinstance(kf, str):
                continue
            parts = [k for k in kf.kids if not isinstance(k, str)]
            items.append(f"- **{inline(parts[0]).strip()}** {inline(parts[1]).strip()}")
        out.append("\n".join(items))
    elif t == "dl":
        items, dt = [], None
        for k in n.kids:
            if isinstance(k, str):
                continue
            if k.tag == "dt":
                dt = inline(k).strip()
            elif k.tag == "dd":
                items.append(f"- **{dt}:** {inline(k).strip()}")
        out.append("\n".join(items))
    elif t == "table":
        rows = []
        for tr in iter_tags(n, "tr"):
            rows.append([inline(cell).strip().replace("|", "\\|") for cell in tr.kids
                         if not isinstance(cell, str) and cell.tag in ("th", "td")])
        out.append(table_md(rows))
    elif t == "figure":
        img = next(iter_tags(n, "img"))
        name = re.match(r"__IMG_(\w+)__", img.attrs["src"]).group(1)
        cap = inline(next(iter_tags(n, "figcaption"))).strip()
        out.append(f"![{img.attrs.get('alt', '')}](plots/{name}.png)\n\n{cap}")
    elif t == "blockquote":
        out.append("\n>\n".join("> " + inline(p).strip() for p in n.kids if not isinstance(p, str) and p.tag == "p"))
    elif "model" in c:
        out.append("  \n".join("> " + l.strip() for l in inline(n).split("  \n")))
    elif "tablebox" in c and "__CALIB_TABLE__" in inline(n):
        out.append(calib_table())
    else:
        for k in n.kids:
            block(k)


def iter_tags(n, tag):
    for k in n.kids:
        if isinstance(k, str):
            continue
        if k.tag == tag:
            yield k
        yield from iter_tags(k, tag)


src = (here / "report_template.html").read_text()
src = src[src.index('<div class="wrap">'):]
src = src[:src.index("<h2>Files</h2>")]          # replaced by readme_code_section.md
p = Tree(); p.feed(src)
block(p.root)

md = "\n\n".join(x for x in out if x.strip())
# title first, then the date line
m = re.match(r"(\*Simulation report[^\n]*\*)\n\n(# [^\n]*)", md)
if m:
    md = m.group(2) + "\n\n" + m.group(1) + md[m.end():]
intro = ("The same report as a web page, with figures that follow the light or dark theme: "
         "[report/rda_repeated_report.html](report/rda_repeated_report.html) (download and open in a browser).")
# put the pointer to the HTML version right after the header block
md = md.replace("\n\n## Summary", "\n\n" + intro + "\n\n## Summary", 1)
md += "\n\n## Code\n\n" + (here / "readme_code_section.md").read_text().replace("## Files", "### Files").replace("## Running", "### Running")
(sim / "README.md").write_text(md.rstrip() + "\n")
print("wrote", sim / "README.md", len(md.splitlines()), "lines")
