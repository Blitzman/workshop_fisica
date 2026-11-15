"""Pairing gap models and superfluid functions.

CATALOGUE OF GAP MODELS (decisions D-033 and D-034)
Every available choice is an entry of GAP_CATALOGUE, numbered as follows:

    1-7     crust neutrons, 1S0   local table
    8-10    core protons, 1S0     local table
    11-16   core neutrons, 3P2    local table
    17-25   core protons, 1S0     Ho, Elshamouty, Heinke and Potekhin 2015 (Ho15), Table II: AO, BCLL, BS, CCDK, CCYms, CCYps, EEHO, EEHOr, T
    26-33   core neutrons, 3P2    Ho15 Table II: AO, BEEHS, EEHO, EEHOr, SYHHP, T, TTav, TToa
    34-36   crust neutrons, 1S0   Ho15 Table II, three of its nine neutron-singlet rows: SCLBL (34), SFB (35), WAP (36). The other six (AWP2,
                                  AWP3, CCDK, CLS, GIPSF, MSH) are not included: the local models and these three cover the range (D-035).

The zero-temperature gap is Delta(k) = Delta0 (k-k0)^2/((k-k0)^2+k1) (k-k2)^2/((k-k2)^2+k3) for k0 < k < k2 and zero outside
(Andersson et al. 2005, Eq. 19; Ho15, Eq. 2); k is the Fermi momentum [fm^-1] of the paired species. The critical temperature is
k_B T_c = r Delta. The ratio r belongs to the entry, because the sources tabulate Delta in two conventions:
- Singlet: r = 0.5669 for every entry.
- Triplet, local table: r = 0.8416, the weak-coupling value for the minimum (equatorial) gap Delta_B of the type-B state, whose gap is
  Delta_B sqrt(1 + 3 cos^2 theta) (Yakovlev et al. 2001, Table 5 and Eqs. 183-188; Andersson et al. 2005, Eq. 20). The deep fit of
  Ho et al. 2012 (model 12) prints 0.8418; the code uses 0.8416 for it (0.024 % in T_c).
- Triplet, Ho15 entries: r = 0.1187 (Ho15, Eq. 1), for their own amplitude Delta_H, which is about 7.09 times Delta_B.
There is no basis in the fits for applying 0.1187 to the local triplet models (literature audit of 15 November 2026). The convention
matters only through T_c: the dimensionless gap v(T/T_c) below is the one of Yakovlev et al. 2001 and does not depend on how a table
normalises Delta.

Default selection of the input files (crust neutrons 7, core neutrons 12, core protons 8): the local entries with the values and ratios
that the code has always used. The input file takes the number of the model; a model can be used only in its own region.
``python3 superfluid.py`` prints the whole catalogue with the source of every entry.

The local numbering is that of the original implementation, not the selectors of Page's NSCool, whose documented numbers differ.
"""

from dataclasses import dataclass

import numpy as np
from constants import BOLTZMANN_MEV_PER_K
from numerical_validation import checked

# Region of a gap model: the key of the input file and the pairing type (1 = singlet, 2 = triplet) of its dimensionless gap.
REGIONS = {
    "crust_n": {"pairing_type": 1, "description": "crust neutrons, 1S0"},
    "core_n": {"pairing_type": 2, "description": "core neutrons, 3P2"},
    "core_p": {"pairing_type": 1, "description": "core protons, 1S0"},
}
RATIO_SINGLET = 0.5669  # k_B T_c / Delta, every singlet entry
RATIO_TRIPLET_LOCAL = 0.8416  # k_B T_c / Delta_B, triplet entries of the local table
RATIO_TRIPLET_HO = 0.1187  # k_B T_c / Delta_H, triplet entries of Ho et al. 2015 (their Eq. 1)


@dataclass(frozen=True)
class GapModel:
    """One gap model: parameters of Delta(k) and the ratio k_B T_c / Delta of its convention."""

    number: int
    region: str  # key of REGIONS
    label: str
    parameters: tuple  # (Delta0 [MeV], k0 [fm^-1], k1 [fm^-2], k2 [fm^-1], k3 [fm^-2])
    tc_ratio: float
    source: str
    remark: str = ""  # one line, for the README table
    note: str = ""  # what is known about the entry

    @property
    def pairing_type(self):
        return REGIONS[self.region]["pairing_type"]


_ACG = "Andersson, Comer and Glampedakis 2005, Nucl. Phys. A 763, 212 (ACG05), Table I, model {}"
_HGA = "Ho, Glampedakis and Andersson 2012, MNRAS 422, 2632 (HGA12), Table 1, {}"
_HO = "Ho, Elshamouty, Heinke and Potekhin 2015, Phys. Rev. C 91, 015806 (arXiv:1412.7759), Table II, entry '{}'; k_B T_c = r Delta (Eq. 1)"
_ELG = "Elgaroy et al. 1996, Nucl. Phys. A 607, 425"

# (number, region, label, parameters, ratio, source, remark, note) of the local table, as always used by the code.
# Sources and conventions: literature audit of 15 November 2026 (fit tuples and conversions checked against ACG05 and HGA12).
_LOCAL_MODELS = [
    (1, "crust_n", "local 1", (9.3, 0.02, 0.6, 1.55, 0.32), RATIO_SINGLET, _ACG.format("A") + " (bare interaction)", "ACG05 A", ""),
    (2, "crust_n", "local 2", (68.0, 0.1, 4.0, 1.7, 4.0), RATIO_SINGLET, _ACG.format("a") + "; also " + _HGA.format("crust fit"),
     "ACG05 a = HGA12 crust; = 7", "alias: the parameters are identical to model 7"),
    (3, "crust_n", "local 3", (4.0, 0.4, 1.5, 1.65, 0.05), RATIO_SINGLET, _ACG.format("b"), "ACG05 b", "Ho15 SCLBL (model 34) is a different fit"),
    (4, "crust_n", "local 4", (22.0, 0.3, 0.09, 1.05, 4.0), RATIO_SINGLET, _ACG.format("c"), "ACG05 c", ""),
    (5, "crust_n", "local 5", (2.9, 0.3, 0.017, 1.3, 0.07), RATIO_SINGLET, _ACG.format("d"), "ACG05 d", ""),
    (6, "crust_n", "local 6", (82.11, 0.05, 2.31, 1.63, 5.81), RATIO_SINGLET, "unidentified: the fit source was not established by the audit",
     "source unresolved", "no identical Ho15 Table II tuple (Ho15 SFB, model 35, is (45, 0.10, 4.5, 1.55, 2.5)); an association needs further evidence"),
    (7, "crust_n", "local 7", (68.0, 0.1, 4.0, 1.7, 4.0), RATIO_SINGLET, _ACG.format("a") + "; also " + _HGA.format("crust fit"),
     "default; ACG05 a; = 2", "default of the input files; alias: the parameters are identical to model 2"),
    (8, "core_p", "local 8", (120.0, 0.0, 9.0, 1.3, 1.8), RATIO_SINGLET, _HGA.format("proton fit (Chen-like)"),
     "default; HGA12; ~ 20 (CCDK)",
     "default of the input files; Ho15 CCDK (model 20) differs in the amplitude (120 against 102) and the last offset (1.8 against 1.5); peak T_c 6.795e9 K against 6.706e9 K"),
    (9, "core_p", "local 9", (55.0, 0.15, 4.0, 1.27, 4.0), RATIO_SINGLET,
     _ACG.format("f") + " (Amundsen and Ostgaard 1985, Nucl. Phys. A 437, 487; OPEG, m*/m = 0.6); also Aguilera, Pons and Miralles 2008, A&A 486, 255, Table 2",
     "ACG05 f", "a distinct fit from Ho15 AO (model 17), which shares the attribution to Amundsen and Ostgaard: peaks 1.912e9 K and 2.433e9 K"),
    (10, "core_p", "local 10", (2.27, 0.1, 0.07, 1.05, 0.25), RATIO_SINGLET, _ACG.format("g"), "ACG05 g; ~ 18 (BCLL)",
     "Ho15 BCLL (model 18) differs in the amplitude (2.27 against 1.69), the lower end (0.10 against 0.05) and the last offset (0.25 against 0.16); peaks 5.691e9 K and 5.407e9 K"),
    (11, "core_n", "local 11", (0.068, 1.28, 0.1, 2.37, 0.02), RATIO_TRIPLET_LOCAL, _ACG.format("l") + " (" + _ELG + ", nonrelativistic); also " + _HGA.format("shallow fit"),
     "ACG05 l = HGA12 shallow; = 15",
     "alias: the parameters are identical to model 15. Ho15 EEHO (model 28) has the same shape with the amplitude 0.48; converted by 0.1187/0.8416 it gives 0.0677 MeV, 0.45 % below the amplitude here"),
    (12, "core_n", "local 12", (0.15, 2.0, 0.1, 3.1, 0.02), RATIO_TRIPLET_LOCAL, _HGA.format("deep fit") + " (empirical, motivated by cooling studies); its conversion is printed as 0.8418",
     "default; HGA12 deep",
     "default of the input files; Ho15 (note d of its Table II) says that SYHHP (model 30) replaces this model; the code uses 0.8416, 0.024 % below the printed conversion; the peak T_c 1.082e9 K (at k_F = 2.69 fm^-1, outside the 1.30-2.41 range of the example star) is above the 0.88e9 K of the Ho15 triplets"),
    (13, "core_n", "local 13", (2.2, 1.05, 1.0, 2.82, 0.6), RATIO_TRIPLET_LOCAL, _ACG.format("j") + " (Baldo et al. 1998, Phys. Rev. C 58, 1921; Paris potential, free spectrum)",
     "ACG05 j (Baldo et al.)", "a large-gap published example (peak T_c 5.39e9 K with 0.8416); Ho15 BEEHS cites the same paper with a BHF spectrum and a peak of 0.388e9 K"),
    (14, "core_n", "local 14", (0.425, 1.1, 0.5, 2.7, 0.5), RATIO_TRIPLET_LOCAL, _ACG.format("k") + " (" + _ELG + ", nonrelativistic)", "ACG05 k (Elgaroy et al.)",
     "not to be identified with Ho15 TTav on resemblance"),
    (15, "core_n", "local 15", (0.068, 1.28, 0.1, 2.37, 0.02), RATIO_TRIPLET_LOCAL, _ACG.format("l") + "; also " + _HGA.format("shallow fit"), "alias of 11",
     "alias: the parameters are identical to model 11"),
    (16, "core_n", "local 16", (2.9, 1.21, 0.5, 1.62, 0.5), RATIO_TRIPLET_LOCAL, _ACG.format("m") + " (" + _ELG + ", relativistic)", "ACG05 m (Elgaroy et al.)",
     "Ho15 EEHOr (model 29) is a refit, (0.23, 1.2, 0.026, 1.6, 0.008), with other shape and amplitude; peaks 0.170e9 K and 0.166e9 K"),
]

# Ho et al. 2015, Table II: Delta0 [MeV], k0 [fm^-1], k1 [fm^-2], k2 [fm^-1], k3 [fm^-2]. The proton and triplet rows were transcribed in the
# old GitHub code and all 17 reproduce Table II (literature audit); the three neutron-singlet rows included are taken from the audit.
_HO_PROTONS = {
    "AO": (14, 0.15, 0.22, 1.05, 3.8), "BCLL": (1.69, 0.05, 0.07, 1.05, 0.16), "BS": (17, 0.0, 2.9, 0.8, 0.08),
    "CCDK": (102, 0.0, 9.0, 1.3, 1.5), "CCYms": (35, 0.0, 5.0, 1.1, 0.5), "CCYps": (34, 0.0, 5.0, 0.95, 0.3),
    "EEHO": (4.5, 0.0, 0.57, 1.2, 0.35), "EEHOr": (61, 0.0, 6.0, 1.1, 0.6), "T": (48, 0.15, 2.1, 1.2, 2.8),
}
_HO_TRIPLETS = {
    "AO": (4.0, 1.2, 0.45, 3.3, 5.0), "BEEHS": (0.45, 1.0, 0.40, 3.2, 0.25), "EEHO": (0.48, 1.28, 0.1, 2.37, 0.02),
    "EEHOr": (0.23, 1.2, 0.026, 1.6, 0.0080), "SYHHP": (1.0, 2.08, 0.04, 2.7, 0.013), "T": (1.2, 1.55, 0.05, 2.35, 0.07),
    "TTav": (3.0, 1.1, 0.60, 2.92, 3.0), "TToa": (2.1, 1.1, 0.60, 3.2, 2.4),
}
# Three of the nine neutron-singlet rows of Ho15 Table II (the others are not included, decision D-035).
_HO_NEUTRON_SINGLET = {"SCLBL": (4.1, 0.35, 1.7, 1.67, 0.06), "SFB": (45, 0.10, 4.5, 1.55, 2.5), "WAP": (69, 0.15, 3, 1.4, 3)}

GAP_CATALOGUE = {}
for _n, _region, _label, _par, _ratio, _source, _remark, _note in _LOCAL_MODELS:
    GAP_CATALOGUE[_n] = GapModel(_n, _region, _label, _par, _ratio, _source, _remark, _note)
for _i, (_name, _par) in enumerate(_HO_PROTONS.items()):
    GAP_CATALOGUE[17 + _i] = GapModel(17 + _i, "core_p", f"Ho2015 {_name}", _par, RATIO_SINGLET, _HO.format(_name))
for _i, (_name, _par) in enumerate(_HO_TRIPLETS.items()):
    GAP_CATALOGUE[26 + _i] = GapModel(26 + _i, "core_n", f"Ho2015 {_name}", _par, RATIO_TRIPLET_HO, _HO.format(_name))
for _i, (_name, _par) in enumerate(_HO_NEUTRON_SINGLET.items()):
    GAP_CATALOGUE[34 + _i] = GapModel(34 + _i, "crust_n", f"Ho2015 {_name}", _par, RATIO_SINGLET, _HO.format(_name),
                                      "Ho15 neutron singlet", "tuple taken from the literature audit, which checked it against Table II")

# Compatibility: the parameters by model number (the tuple order is Delta0, k0, k1, k2, k3) and the numbers of each region.
GAP_MODELS = {number: model.parameters for number, model in GAP_CATALOGUE.items()}
GAP_MODEL_NUMBERS = {region: tuple(n for n, m in GAP_CATALOGUE.items() if m.region == region) for region in REGIONS}
# k_B T_c / Delta of the local table by pairing type (type 3 is not used). The ratio used for a model is GapModel.tc_ratio.
# 0.5669 and 0.8416 follow Yakovlev et al. (2001), Table 5, and Andersson et al. (2005), Eq. (20).
TC_GAP_RATIO = {1: RATIO_SINGLET, 2: RATIO_TRIPLET_LOCAL, 3: 0.4926}


def fermi_quantities(rho, particle_fraction):
    """Return the particle number density [cm^-3] and the Fermi momentum kF [fm^-1]."""
    rho, particle_fraction = np.broadcast_arrays(
        np.asarray(rho, float), np.asarray(particle_fraction, float)
    )
    if not (np.isfinite(rho).all() and np.isfinite(particle_fraction).all()):
        raise ValueError("rho and particle_fraction must be finite")
    if np.any(rho < 0) or np.any((particle_fraction < 0) | (particle_fraction > 1)):
        raise ValueError("rho must be non-negative and particle_fraction in [0, 1]")
    with np.errstate(over="raise", invalid="raise"):
        number_density = rho / 1.66e-24 * particle_fraction
        kF = (3.0 * np.pi**2 * number_density) ** (1.0 / 3.0) * 1e-13
    return number_density, kF


def pairing_gap_from_kF(gap_model, kF):
    """Zero-temperature gap Delta(T=0) in MeV, with each model's own normalization."""
    if gap_model not in GAP_CATALOGUE:
        raise ValueError(f"Unknown gap model {gap_model}; the models are {sorted(GAP_CATALOGUE)}")
    scale, k1, k2, k3, k4 = GAP_CATALOGUE[gap_model].parameters
    kF = np.asarray(kF, float)
    if not np.isfinite(kF).all() or np.any(kF < 0):
        raise ValueError("kF must be finite and non-negative")
    if scale < 0 or k2 <= 0 or k4 <= 0 or k1 >= k3:
        raise ValueError("Invalid gap-model parameters")
    gap_MeV = np.zeros_like(kF)
    mask = (kF > k1) & (kF <= k3)
    left, right = (kF[mask] - k1) ** 2, (kF[mask] - k3) ** 2
    gap_MeV[mask] = scale * (left / (left + k2)) * (right / (right + k4))
    return gap_MeV


def pairing_gap(gap_model, pairing_type, rho, particle_fraction):
    """Return (Delta0 [MeV], Tc [K]) for a gap model and pairing channel."""
    if pairing_type not in TC_GAP_RATIO:
        raise ValueError("Unknown pairing type")
    _, kF = fermi_quantities(rho, particle_fraction)
    gap_MeV = pairing_gap_from_kF(gap_model, kF)
    model = GAP_CATALOGUE[gap_model]
    if model.pairing_type != pairing_type:
        raise ValueError(f"Gap model {gap_model} ({model.label}) is a {REGIONS[model.region]['description']} gap, not pairing type {pairing_type}")
    return gap_MeV, model.tc_ratio * gap_MeV / BOLTZMANN_MEV_PER_K


def pairing_critical_temperatures(cfg, star):
    """Critical temperatures [K] at the faces of crust neutrons, core neutrons and core protons.

    ``cfg`` holds the on/off switches and gap models of the input file; ``star``
    is a dictionary with the arrays ``rho``, ``yn``, ``yp`` and ``xh`` at the faces. Returns ``tc_crust``, ``tc_n`` and
    ``tc_p`` (zero outside the region of each species or when it is switched off).
    """
    result = {}
    for key, name, pairing_type, fraction, region in (
        ("crust_n", "tc_crust", 1, star["yn"], star["xh"] > 0.0),
        ("core_n", "tc_n", 2, star["yn"], star["xh"] == 0.0),
        ("core_p", "tc_p", 1, star["yp"], star["xh"] == 0.0),
    ):
        model = cfg[key]["model"]
        if cfg[key]["enabled"] and model in GAP_CATALOGUE and GAP_CATALOGUE[model].region != key:
            raise ValueError(f"Gap model {model} ({GAP_CATALOGUE[model].label}) cannot be used for {REGIONS[key]['description']}; "
                             f"the models of that region are {list(GAP_MODEL_NUMBERS[key])}")
        _, Tc_K = pairing_gap(model, pairing_type, star["rho"], fraction)
        result[name] = np.where(region & cfg[key]["enabled"], Tc_K, 0.0)
    return result


def dimensionless_pairing_gap(tau, pairing_type):
    """Dimensionless gap v for pairing types 1 and 2, valid for 0 <= tau <= 1.

    For tau <= 0.1 a constant 1e-10 is added to tau, so that tau = 0 is finite. Outside
    this domain the normal phase must be used instead of this function.
    """
    tau = checked(tau, "tau")
    if np.any(tau > 1) or pairing_type not in (1, 2):
        raise ValueError("tau <= 1 and pairing_type 1 or 2 required")
    alpha = 0.5669 if pairing_type == 1 else 0.8416
    v = np.empty_like(tau)
    low = tau <= 0.1
    v[low] = 1.0 / alpha / (tau[low] + 1e-10)
    x = tau[~low]
    if pairing_type == 1:
        v[~low] = np.sqrt(1 - x) * (1.456 - 0.157 / np.sqrt(x) + 1.764 / x)
    else:
        v[~low] = np.sqrt(1 - x) * (0.7893 + 1.188 / x)
    return v


def peak_critical_temperature(gap_model):
    """Largest critical temperature [1e9 K] of a model over the Fermi momentum, and the k_F [fm^-1] where it is reached."""
    model = GAP_CATALOGUE[gap_model]
    k = np.linspace(0.0, 3.6, 3601)
    tc = model.tc_ratio * pairing_gap_from_kF(gap_model, k) / BOLTZMANN_MEV_PER_K / 1e9
    return float(tc.max()), float(k[tc.argmax()])


def describe_gap_models(long=False):
    """The catalogue as text: number, region, label, T_c ratio, peak T_c, k_F range and remark (``long``: sources and notes)."""
    lines = [f"{'No.':>3} {'Region':20s} {'Label':13s} {'r':>6} {'Peak Tc [1e9 K]':>15} {'kF range [fm^-1]':>17}  Remark"]
    for number, model in GAP_CATALOGUE.items():
        peak, _ = peak_critical_temperature(number)
        lines.append(f"{number:3d} {REGIONS[model.region]['description']:20s} {model.label:13s} {model.tc_ratio:6.4f} {peak:15.3f} "
                     f"{model.parameters[1]:7.2f} - {model.parameters[3]:5.2f}  {model.remark}")
        if long:
            lines.append(f"      source: {model.source}")
            if model.note:
                lines.append(f"      note: {model.note}")
    return "\n".join(lines)


if __name__ == "__main__":
    print(describe_gap_models(long=True))
