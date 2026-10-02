"""
Acetone -> Splitter -> Heater (vaporisation totale) - version 2

Ameliorations par rapport a la v1 :
  1. Cp liquide dependant de T (Rowlinson-Bondi), integre entre T_in et Tb
  2. Zone de flash masquee (NaN) et ombree dans les graphiques
  3. Surchauffe optionnelle de la vapeur, perte de charge du Heater,
     travail de pompe optionnel (indicatif, hors puissance du Heater)
  4. Validation np.isfinite de toutes les entrees (nan / inf refuses)

Mode interactif (par defaut) :
    python acetone_heater_param_v2.py
Mode non interactif :
    python acetone_heater_param_v2.py --flow 50 --ratio 0.7 --P 2 --T 25 --superheat 10 --dP 0.1
    python acetone_heater_param_v2.py --flow 38 --ratio 0.6 --P 3 --Psource 1 --no-plot
"""
import argparse
import math
import os
import sys
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

# ---------------- Donnees acetone ----------------
M, TC, OMEGA = 58.08, 508.1, 0.307       # g/mol, K, facteur acentrique
DH_NB        = 29.1e3                    # J/mol a T_NB
A, B, C      = 4.42448, 1312.253, -32.445  # Antoine NIST (bar, K)
ATM          = 1.01325                   # bar
T_MIN_ANT, T_MAX_ANT = 259.0, 507.0      # domaine Antoine (K)
RHO_LIQ      = 790.0                     # kg/m3 (pompe)
R_GAS        = 8.314462618               # J/(mol.K)
CP_LIQ_298   = 2.17e3                    # J/(kg.K) valeur d'ancrage a 298.15 K
T_FREEZE     = 178.5                     # K, point de fusion
TR_MAX       = 0.95                      # Rowlinson-Bondi diverge pres de Tc
# Shomate gaz parfait acetone (NIST, 298-1200 K), t = T/1000, J/(mol.K)
SH = (6.301, 261.1, -150.2, 32.23, 0.1925)


# ---------------- Capacites calorifiques ----------------
def cp_gas_mol(T):
    """Cp gaz parfait [J/(mol.K)] (Shomate)."""
    t = T / 1000.0
    a, b, c, d, e = SH
    return a + b * t + c * t**2 + d * t**3 + e / t**2

def _cp_liq_raw(T):
    """Cp liquide [J/(mol.K)] : Rowlinson-Bondi."""
    Tr = T / TC
    x = 1.0 - Tr
    return cp_gas_mol(T) + R_GAS * (
        1.586 + 0.49 / x
        + OMEGA * (4.2775 + 6.3 * x ** (1 / 3) / Tr + 0.4355 / x))

# Calibrage : on impose Cp(298.15 K) = valeur experimentale (2.17 kJ/kg.K)
_SCALE = CP_LIQ_298 * M * 1e-3 / _cp_liq_raw(298.15)

def cp_liq(T):
    """Cp liquide [J/(kg.K)] dependant de T."""
    return _cp_liq_raw(T) * _SCALE / (M * 1e-3)

def cp_gas(T):
    """Cp vapeur [J/(kg.K)] (gaz parfait)."""
    return cp_gas_mol(T) / (M * 1e-3)

def check_T_liq(T, label):
    if not (T_FREEZE <= T <= TR_MAX * TC):
        raise ValueError(f"{label} = {T - 273.15:.1f} C hors du domaine de validite "
                         f"du Cp liquide ({T_FREEZE - 273.15:.0f} a {TR_MAX * TC - 273.15:.0f} C)")


# ---------------- Thermodynamique ----------------
def p_sat(T):
    return 10 ** (A - B / (T + C))       # bar

def t_boil(P_atm, T_init=330.0, tol=1e-10, itmax=100):
    """Newton-Raphson sur ln Psat(T) - ln P, repli sur Brent. Retourne (Tb [K], n_iter)."""
    lnP = np.log(P_atm * ATM)
    T = T_init
    for i in range(1, itmax + 1):
        if not (T_MIN_ANT - 50 < T < T_MAX_ANT + 50):
            break
        f  = np.log(p_sat(T)) - lnP
        df = np.log(10) * B / (T + C) ** 2
        dT = f / df
        T -= dT
        if abs(dT) < tol:
            if not (T_MIN_ANT <= T <= T_MAX_ANT):
                raise ValueError(f"P = {P_atm} atm hors du domaine d'Antoine (Tb = {T:.1f} K)")
            return T, i
    try:
        T = brentq(lambda t: np.log(p_sat(t)) - lnP, T_MIN_ANT, T_MAX_ANT, xtol=1e-12)
    except ValueError:
        raise ValueError(f"P = {P_atm} atm hors du domaine d'Antoine")
    return T, -1

def t_boil_check(P_atm):
    return brentq(lambda T: p_sat(T) - P_atm * ATM, T_MIN_ANT, T_MAX_ANT, xtol=1e-12)

T_NB = brentq(lambda T: p_sat(T) - ATM, T_MIN_ANT, T_MAX_ANT, xtol=1e-12)

def dh_vap(T):
    """Chaleur latente [J/kg], relation de Watson."""
    if T >= TC:
        raise ValueError("T >= Tc : pas de chaleur latente")
    return DH_NB / (M * 1e-3) * ((1 - T / TC) / (1 - T_NB / TC)) ** 0.38


# ---------------- Modele du procede ----------------
def simulate(flow_kg_h=38.0, ratio_heater=0.60, P_atm=1.0, T_in_C=10.0,
             superheat_K=0.0, dP_atm=0.0, P_source_atm=None, eta_pump=0.6):
    """
    flow_kg_h    : debit massique total avant le splitter [kg/h]
    ratio_heater : fraction massique vers le Heater (0..1)
    P_atm        : pression d'entree du Heater [atm]
    T_in_C       : temperature d'entree [C]
    superheat_K  : surchauffe de la vapeur au-dessus de la saturation [K]
    dP_atm       : perte de charge du Heater [atm] ; la vaporisation est evaluee
                   a la pression moyenne P - dP/2, la sortie est a P - dP
    P_source_atm : pression amont de la pompe [atm] (None = pas de pompe).
                   Travail de pompe indicatif, NON inclus dans Q.
    """
    vals = (flow_kg_h, ratio_heater, P_atm, T_in_C, superheat_K, dP_atm, eta_pump)
    if not all(math.isfinite(v) for v in vals) or \
       (P_source_atm is not None and not math.isfinite(P_source_atm)):
        raise ValueError("toutes les entrees doivent etre des nombres finis")
    if not (0.0 <= ratio_heater <= 1.0):
        raise ValueError("ratio_heater doit etre entre 0 et 1")
    if flow_kg_h < 0 or P_atm <= 0:
        raise ValueError("debit >= 0 et pression > 0 requis")
    if superheat_K < 0 or dP_atm < 0:
        raise ValueError("surchauffe >= 0 et perte de charge >= 0 requis")
    if dP_atm >= P_atm:
        raise ValueError("la perte de charge doit etre < P")
    if not (0 < eta_pump <= 1):
        raise ValueError("rendement de pompe dans ]0, 1]")

    m_heater = flow_kg_h * ratio_heater
    m_bypass = flow_kg_h - m_heater
    P_vap = P_atm - dP_atm / 2.0
    P_out = P_atm - dP_atm
    Tb, n_it = t_boil(P_vap)
    T_in = T_in_C + 273.15
    check_T_liq(T_in, "T entree")
    check_T_liq(Tb, "T ebullition")
    T_out = Tb + superheat_K

    m = m_heater / 3600.0                              # kg/s
    q_sens  = m * quad(cp_liq, T_in, Tb)[0]            # W (Cp(T) integre)
    q_lat   = m * dh_vap(Tb)                           # W
    q_super = m * quad(cp_gas, Tb, T_out)[0] if superheat_K > 0 else 0.0
    w_pump = 0.0
    if P_source_atm is not None and P_atm > P_source_atm:
        w_pump = m / RHO_LIQ * (P_atm - P_source_atm) * ATM * 1e5 / eta_pump
    return dict(flow=flow_kg_h, ratio=ratio_heater, P=P_atm, P_out=P_out, T_in_C=T_in_C,
                m_heater=m_heater, m_bypass=m_bypass,
                Tb_C=Tb - 273.15, T_out_C=T_out - 273.15, n_iter=n_it,
                Q_sens_W=q_sens, Q_lat_W=q_lat, Q_super_W=q_super,
                Q_W=q_sens + q_lat + q_super, W_pump_W=w_pump,
                cp_avg=(q_sens / (m * (Tb - T_in)) if m > 0 and Tb != T_in else cp_liq(T_in)),
                flash_warning=bool(Tb < T_in))


def report(r):
    it = f"{r['n_iter']} it." if r["n_iter"] > 0 else "Brent"
    print("=" * 56)
    print(f" Entree : {r['flow']:.2f} kg/h | {r['T_in_C']:.1f} C | {r['P']:.3f} atm")
    print(f" Splitter : {r['ratio']*100:.1f} % Heater / {(1-r['ratio'])*100:.1f} % bypass")
    print(f"   -> Heater : {r['m_heater']:.3f} kg/h")
    print(f"   -> Bypass : {r['m_bypass']:.3f} kg/h")
    print("-" * 56)
    print(f" T ebullition ({it}) : {r['Tb_C']:.3f} C")
    print(f" Cp liquide moyen    : {r['cp_avg']:.0f} J/(kg.K)")
    if r["flash_warning"]:
        print(" /!\\ T_entree > T_ebullition : le liquide serait en flash a cette")
        print("     pression ; Q non valide (le modele suppose un liquide sous-refroidi).")
    else:
        print(f" Q sensible : {r['Q_sens_W']:9.1f} W")
        print(f" Q latente  : {r['Q_lat_W']:9.1f} W")
        if r["Q_super_W"] > 0:
            print(f" Q surchauffe: {r['Q_super_W']:8.1f} W")
        print(f" Q TOTAL    : {r['Q_W']:9.1f} W = {r['Q_W']/1000:.3f} kW = {r['Q_W']*3.6:.0f} kJ/h")
    print(f" Sortie Heater : {r['T_out_C']:.1f} C, {r['P_out']:.3f} atm")
    if r["W_pump_W"] > 0:
        print(f" Pompe (indicatif, hors Q) : {r['W_pump_W']:.2f} W")
    print("=" * 56)


# ---------------- Sensibilites ----------------
def sweep(base, name, values):
    """Fait varier un parametre. Hors domaine -> NaN ; zone de flash -> Q = NaN."""
    keymap = {"flow": "flow_kg_h", "ratio": "ratio_heater", "P": "P_atm", "T": "T_in_C"}
    out = []
    for v in values:
        kw = dict(base); kw[keymap[name]] = float(v)
        try:
            r = simulate(**kw)
            if r["flash_warning"]:
                r["Q_W"] = np.nan
            out.append(r)
        except ValueError:
            out.append(dict(Tb_C=np.nan, Q_W=np.nan, flash_warning=False))
    return out

def _shade(a, x, res):
    mask = np.array([bool(r["flash_warning"]) for r in res])
    if mask.any():
        a.fill_between(x, 0, 1, where=mask, transform=a.get_xaxis_transform(),
                       color="red", alpha=0.12)
    return mask.any()

def plot_all(base, Pmin=0.1, Pmax=3.0, npts=60, show=True):
    import matplotlib
    if not show or (sys.platform.startswith("linux")
                    and not os.environ.get("DISPLAY")
                    and not os.environ.get("WAYLAND_DISPLAY")):
        matplotlib.use("Agg")
        show = False
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    flash_patch = Patch(color="red", alpha=0.12, label="zone de flash (Q non valide)")

    P = np.linspace(Pmin, Pmax, npts)
    res = sweep(base, "P", P)
    a = ax[0, 0]
    l1, = a.plot(P, [r["Tb_C"] for r in res], "b-", label="T ebullition")
    l2 = a.axhline(base["T_in_C"], color="r", ls=":", label="T entree")
    a.set(xlabel="Pression (atm)", ylabel="T ebullition (C)",
          title=f"Tb et Q = f(P)   [{Pmin} -> {Pmax} atm]")
    a.grid(alpha=.3)
    flashed = _shade(a, P, res)
    a2 = a.twinx()
    l3, = a2.plot(P, [r["Q_W"] / 1000 for r in res], "g--", label="Q (kW)")
    a2.set_ylabel("Q (kW)", color="g")
    a.legend(handles=[l1, l2, l3] + ([flash_patch] if flashed else []), loc="best")

    x = np.linspace(0, 1, 51)
    res = sweep(base, "ratio", x)
    a = ax[0, 1]; a.plot(x * 100, [r["Q_W"] / 1000 for r in res], "g-")
    a.set(xlabel="% du flux vers le Heater", ylabel="Q (kW)",
          title="Q = f(ratio splitter)"); a.grid(alpha=.3)

    F = np.linspace(1, 2 * max(base["flow_kg_h"], 1), 50)
    res = sweep(base, "flow", F)
    a = ax[1, 0]; a.plot(F, [r["Q_W"] / 1000 for r in res], "g-")
    a.set(xlabel="Debit total (kg/h)", ylabel="Q (kW)", title="Q = f(debit)"); a.grid(alpha=.3)

    T = np.linspace(-20, 50, 50)
    res = sweep(base, "T", T)
    a = ax[1, 1]; a.plot(T, [r["Q_W"] / 1000 for r in res], "g-")
    a.set(xlabel="T entree (C)", ylabel="Q (kW)", title="Q = f(T entree)"); a.grid(alpha=.3)
    if _shade(a, T, res):
        a.legend(handles=[flash_patch], loc="best")

    fig.suptitle(f"Sensibilites autour de : {base['flow_kg_h']} kg/h, ratio {base['ratio_heater']}, "
                 f"{base['P_atm']} atm, {base['T_in_C']} C"
                 + (f", surchauffe {base['superheat_K']} K" if base.get("superheat_K") else ""))
    plt.tight_layout()
    plt.savefig("sensibilites_acetone.png", dpi=150)
    print(" Figure enregistree : sensibilites_acetone.png")
    if show:
        plt.show()
    plt.close(fig)


# ---------------- Saisie terminal ----------------
def ask(prompt, default, cast=float, check=None, err="valeur invalide"):
    """Entree = defaut ; redemande tant que invalide (nan/inf refuses)."""
    while True:
        try:
            s = input(f"  {prompt} [{default}] : ").strip().replace(",", ".")
        except EOFError:
            raise SystemExit("\nEntree interrompue.")
        if s == "":
            return default
        try:
            v = cast(s)
        except ValueError:
            print("   -> Veuillez entrer un nombre.")
            continue
        if not math.isfinite(v):
            print("   -> Valeur non finie refusee.")
            continue
        if check is not None and not check(v):
            print(f"   -> {err}")
            continue
        return v

def ask_yes_no(prompt, default=True):
    d = "O/n" if default else "o/N"
    try:
        s = input(f"  {prompt} [{d}] : ").strip().lower()
    except EOFError:
        raise SystemExit("\nEntree interrompue.")
    if s == "":
        return default
    return s in ("o", "oui", "y", "yes")

def interactive():
    print("=" * 56)
    print("  ACETONE : Splitter + Heater (vaporisation totale)")
    print("  (Entree = valeur par defaut entre crochets)")
    print("=" * 56)
    while True:
        print("\n-- Flux entrant --")
        flow  = ask("Debit massique total [kg/h]", 38.0, check=lambda v: v >= 0,
                    err="le debit doit etre >= 0")
        P     = ask("Pression [atm]", 1.0, check=lambda v: v > 0, err="la pression doit etre > 0")
        T     = ask("Temperature d'entree [C]", 10.0, check=lambda v: -90 < v < 200,
                    err="temperature hors plage (-90 a 200 C)")
        print("\n-- Splitter --")
        ratio = ask("Fraction envoyee au Heater (0 a 1)", 0.60,
                    check=lambda v: 0 <= v <= 1, err="le ratio doit etre entre 0 et 1")
        print("\n-- Options Heater (0 = ignore) --")
        sh = ask("Surchauffe de la vapeur [K]", 0.0, check=lambda v: v >= 0, err="doit etre >= 0")
        dP = ask("Perte de charge du Heater [atm]", 0.0, check=lambda v: 0 <= v < P,
                 err=f"doit etre entre 0 et {P}")
        Ps = ask("Pression amont de la pompe [atm] (= P : pas de pompe)", P,
                 check=lambda v: v > 0, err="doit etre > 0")

        try:
            r = simulate(flow, ratio, P, T, sh, dP, Ps)
        except ValueError as e:
            print(f"\n Erreur : {e}\n Veuillez recommencer.")
            continue
        print()
        report(r)
        print(f" Verification Brent : Tb = {t_boil_check(P - dP / 2) - 273.15:.3f} C")

        if ask_yes_no("\n  Tracer le diagramme de sensibilite ?", True):
            print("\n-- Domaine de variation de la pression --")
            Pmin = ask("Pression min [atm]", 0.1, check=lambda v: v > 0, err="Pmin doit etre > 0")
            Pmax = ask("Pression max [atm]", 3.0, check=lambda v: v > Pmin,
                       err=f"Pmax doit etre > Pmin ({Pmin})")
            npts = ask("Nombre de points", 60, cast=int, check=lambda v: v >= 2,
                       err="au moins 2 points")
            plot_all(dict(flow_kg_h=flow, ratio_heater=ratio, P_atm=P, T_in_C=T,
                          superheat_K=sh, dP_atm=dP, P_source_atm=Ps),
                     Pmin=Pmin, Pmax=Pmax, npts=npts)

        if not ask_yes_no("\n  Nouvelle simulation ?", False):
            print("\n Au revoir.")
            break


# ---------------- Main ----------------
def _finite(s):
    v = float(s)
    if not math.isfinite(v):
        raise argparse.ArgumentTypeError("valeur non finie")
    return v

def cli(argv):
    ap = argparse.ArgumentParser(description="Acetone : splitter + Heater")
    ap.add_argument("--flow",  type=_finite, default=38.0, help="debit total [kg/h]")
    ap.add_argument("--ratio", type=_finite, default=0.60, help="fraction vers le Heater (0-1)")
    ap.add_argument("--P",     type=_finite, default=1.0,  help="pression d'entree Heater [atm]")
    ap.add_argument("--T",     type=_finite, default=10.0, help="temperature d'entree [C]")
    ap.add_argument("--superheat", type=_finite, default=0.0, help="surchauffe vapeur [K]")
    ap.add_argument("--dP",    type=_finite, default=0.0,  help="perte de charge Heater [atm]")
    ap.add_argument("--Psource", type=_finite, default=None, help="pression amont pompe [atm]")
    ap.add_argument("--Pmin",  type=_finite, default=0.1,  help="borne min balayage P [atm]")
    ap.add_argument("--Pmax",  type=_finite, default=3.0,  help="borne max balayage P [atm]")
    ap.add_argument("--npts",  type=int,   default=60,   help="nombre de points Tb=f(P)")
    ap.add_argument("--no-plot", action="store_true", help="ne pas tracer les courbes")
    args = ap.parse_args(argv)

    if args.flow < 0:
        ap.error("--flow doit etre >= 0")
    if not (0.0 <= args.ratio <= 1.0):
        ap.error("--ratio doit etre entre 0 et 1")
    if args.P <= 0:
        ap.error("--P doit etre > 0")
    if args.Pmin <= 0 or args.Pmax <= 0:
        ap.error("Pmin et Pmax doivent etre > 0")
    if args.Pmin >= args.Pmax:
        ap.error("il faut Pmin < Pmax")
    if args.npts < 2:
        ap.error("npts doit etre >= 2")

    try:
        r = simulate(args.flow, args.ratio, args.P, args.T,
                     args.superheat, args.dP, args.Psource)
    except ValueError as e:
        raise SystemExit(f"Erreur : {e}")
    report(r)
    print(f" Verification Brent : Tb = {t_boil_check(args.P - args.dP / 2) - 273.15:.3f} C")
    print(f" Domaine P balaye   : {args.Pmin} -> {args.Pmax} atm ({args.npts} points)")
    if not args.no_plot:
        plot_all(dict(flow_kg_h=args.flow, ratio_heater=args.ratio, P_atm=args.P,
                      T_in_C=args.T, superheat_K=args.superheat, dP_atm=args.dP,
                      P_source_atm=args.Psource),
                 Pmin=args.Pmin, Pmax=args.Pmax, npts=args.npts)


def main():
    if len(sys.argv) > 1:
        cli(sys.argv[1:])
    else:
        try:
            interactive()
        except KeyboardInterrupt:
            print("\n\n Interrompu. Au revoir.")


if __name__ == "__main__":
    main()
