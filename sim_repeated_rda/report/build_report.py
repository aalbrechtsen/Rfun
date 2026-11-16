"""Fill report_template.html with plots and the study-matched table -> rda_repeated_report.html
Run plot_sim_results.R first (from sim_repeated_rda/)."""
import base64, csv, pathlib, re

here = pathlib.Path(__file__).resolve().parent
sim = here.parent
html = (here / "report_template.html").read_text()

# figures are inlined as SVG markup (the artifact viewer does not load image files)
def inline_svg(name, alt, k):
    svg = (sim / "plots" / f"{name}.svg").read_text()
    svg = svg[svg.index("<svg"):]
    svg = re.sub(r"id='(cp\w+)'", rf"id='f{k}\1'", svg)
    svg = re.sub(r"url\(#(cp\w+)\)", rf"url(#f{k}\1)", svg)
    svg = re.sub(r"<svg ([^>]*?)width='[^']*' height='[^']*'",
                 rf"<svg \1role='img' aria-label='{alt}' style='width:100%;height:auto;display:block'", svg, count=1)
    return svg

for k, (name, alt) in enumerate(re.findall(r'<img class="plot" src="__IMG_(\w+?)__" alt="([^"]*)">', html)):
    html = html.replace(f'<img class="plot" src="__IMG_{name}__" alt="{alt}">', inline_svg(name, alt, k))

rows = list(csv.DictReader(open(sim / "sim_dubois_calibrated_results.tsv"), delimiter="\t"))
labels = {"diet_naive": "Diet · RDA (Dubois et al.)", "diet_means": "Diet · individual means",
          "dietGeo_naive": "Diet | geography · RDA (Dubois et al.)",
          "dietGeo_block": "Diet | geography · RDA (ID permutation)",
          "dietGeo_means": "Diet | geography · individual means",
          "nunavut_naive": "Nunavut only · RDA (Dubois et al.)",
          "nunavut_block": "Nunavut only · RDA (ID permutation)"}
cols = []
for d in ("euclid", "unifrac"):
    for r in rows:
        key = (d, r["rho_d"], r["R2true"])
        if r["dist"] == d and key not in cols:
            cols.append(key)
t = ['<p class="meta" style="margin:0 0 6px">Table 1. Share of datasets with p &lt; 0.05 in the study-matched '
     'simulation. Columns: distance (E = Euclidean, U = UniFrac), diet within-person correlation ρ<sub>d</sub>, '
     'and true diet R² where it is not 0. Diet-alone rows test a real association, since diet and microbiome both '
     'differ by region. Paper values: diet 0.11, diet | geography 0.064, Nunavut only 0.169 (all p ≤ 0.004). Red: false-positive rate above 10% in a test whose null is true.</p>',
     "<table><thead><tr><th>Analysis</th>"]
for d, rd, tr in cols:
    t.append(f"<th>{'E' if d == 'euclid' else 'U'} ρ<sub>d</sub>={float(rd):g}"
             + (f"<br>R²={float(tr):g}" if float(tr) > 0 else "") + "</th>")
t.append("</tr></thead><tbody>")
for m, lab in labels.items():
    t.append(f"<tr><td>{lab}</td>")
    for d, rd, tr in cols:
        r = next(x for x in rows if (x["dist"], x["rho_d"], x["R2true"]) == (d, rd, tr) and x["method"] == m)
        rej = float(r["reject"])
        bad = float(tr) == 0 and not m.startswith("diet_") and rej > 0.1
        t.append(f'<td class="{"bad" if bad else ""}">{round(rej * 100)}%</td>')
    t.append("</tr>")
t.append("</tbody></table>")
html = html.replace("__CALIB_TABLE__", "".join(t))
assert "__IMG_" not in html and "__CALIB" not in html
(here / "rda_repeated_report.html").write_text(html)
print("wrote", here / "rda_repeated_report.html", len(html) // 1024, "KB")
