"""Extract the cooling sample of ns_thermal_catalogue.html into a machine-readable JSON file.

The HTML catalogue is the master list; this script only reads it (requires beautifulsoup4). Usage:
    python3 extract_thermal_sample.py ns_thermal_catalogue.html ns_thermal_sample.json
"""
import json, re, sys
from bs4 import BeautifulSoup

NUM = r"\d[\d,]*(?:\.\d+)?"

def number(text):
    return float(text.replace(",", ""))

def parse_interval(text):
    """'3.7–5.2' -> (3.7, 5.2, False); '< 0.3' or '<0.7' -> (None, 0.3, True); '79.43' -> (79.43, 79.43, False)."""
    t = text.strip().replace("−", "-").replace("–", "-").replace("—", "-")
    m = re.fullmatch(rf"<\s*({NUM})", t)
    if m:
        return None, number(m.group(1)), True
    m = re.fullmatch(rf"({NUM})\s*-\s*({NUM})", t)
    if m:
        return number(m.group(1)), number(m.group(2)), False
    m = re.fullmatch(rf"({NUM})", t)
    if m:
        v = number(m.group(1))
        return v, v, False
    raise ValueError(f"cannot parse interval: {text!r}")

CHOSEN = re.compile(r"^Cooling sample · (?P<age>.+?) kyr \( (?P<age_ref>[A-Z0-9]+) \) · (?P<lum>.+?) × 10³² erg/s \( (?P<lum_ref>[A-Z0-9]+) \) · (?P<note>.*?)\.?$")

def extract(path):
    soup = BeautifulSoup(open(path, encoding="utf-8").read(), "html.parser")
    for t in soup(["style", "script"]):
        t.decompose()
    sources, sample = [], []
    for d in soup.find_all("details", class_="source"):
        name_span = d.find("span", class_="name")
        name = name_span.contents[0].strip() if name_span.contents else ""
        small = name_span.find("small")
        entry = {
            "id": int(d["id"].split("-")[1]),
            "name": name,
            "subtitle": small.get_text(" ", strip=True) if small else "",
            "type": d["data-type"],
            "status": d["data-status"],
            "references": d["data-refs"].split("|"),
            "issues": d["data-issues"].split("|") if d["data-issues"] else [],
            "in_cooling_sample": d["data-sample"] == "yes",
        }
        sources.append(entry)
        ch = d.find("p", class_="chosen")
        if ch:
            raw = re.sub(r"\s+", " ", ch.get_text(" ", strip=True))
            m = CHOSEN.match(raw)
            if not m:
                raise ValueError(f"unexpected 'chosen' line for {name}: {raw!r}")
            a_lo, a_hi, a_up = parse_interval(m["age"])
            l_lo, l_hi, l_up = parse_interval(m["lum"])
            sample.append({**entry,
                "age_yr": [None if a_lo is None else a_lo * 1e3, a_hi * 1e3], "age_is_upper_limit": a_up, "age_reference": m["age_ref"],
                "L_inf_erg_s": [None if l_lo is None else l_lo * 1e32, l_hi * 1e32], "L_is_upper_limit": l_up, "L_reference": m["lum_ref"],
                "emission": m["note"], "chosen_text": raw})
    return sources, sample

if __name__ == "__main__":
    src, out = sys.argv[1], sys.argv[2]
    sources, sample = extract(src)
    json.dump({"note": "Derived from ns_thermal_catalogue.html (the master list) by extract_thermal_sample.py. Ages in yr, L_inf in erg/s, as chosen for the cooling sample of the catalogue; bounds keep the conventions of each reference. Working list: to be improved and checked.",
               "n_catalogue_entries": len(sources), "n_cooling_sample": len(sample),
               "cooling_sample": sample, "catalogue_entries": sources}, open(out, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print(f"{len(sources)} catalogue entries, {len(sample)} in the cooling sample -> {out}")
