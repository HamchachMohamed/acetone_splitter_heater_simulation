"""
Acetone -> Splitter -> Heater (vaporisation totale)

Mode interactif (par defaut) : le programme pose les questions dans le terminal.
    python acetone_heater_param.py
(Entree seule = valeur par defaut entre crochets)

Mode non interactif (arguments) :
    python acetone_heater_param.py --flow 50 --ratio 0.7 --P 2 --T 25 --Pmin 0.2 --Pmax 5
    python acetone_heater_param.py --flow 38 --ratio 0.6 --P 1 --T 10 --no-plot
"""
import argparse
import os
import sys
import numpy as np
from scipy.optimize import brentq

# ---------------- Donnees acetone ----------------
M, TC       = 58.08, 508.1               # g/mol, K
DH_NB       = 29.1e3                     # J/mol a T_NB
CP_LIQ      = 2.17e3                     # J/(kg.K)
A, B, C     = 4.42448, 1312.253, -32.445 # Antoine NIST (bar, K)
ATM         = 1.01325                    # bar
T_MIN_ANT, T_MAX_ANT = 259.0, 507.0      # domaine utilisable (K)


# ---------------- Thermodynamique ----------------
def p_sat(T):
    return 10 ** (A - B / (T + C))       # bar

def t_boil(P_atm, T_init=330.0, tol=1e-10, itmax=100):
    """Newton-Raphson sur ln Psat(T) - ln P, avec repli sur Brent.
    Retourne (Tb [K], n_iter)."""
    lnP = np.log(P_atm * ATM)
    T = T_init
    for i in range(1, itmax + 1):
        if not (T_MIN_ANT - 50 < T < T_MAX_ANT + 50):
            break                         # Newton diverge -> repli Brent
        f  = np.log(p_sat(T)) - lnP
        df = np.log(10) * B / (T + C) ** 2
        dT = f / df
        T -= dT
        if abs(dT) < tol:
            if not (T_MIN_ANT <= T <= T_MAX_ANT):
                raise ValueError(
                    f"P = {P_atm} atm hors du domaine de validite d'Antoine "
                    f"(Tb = {T:.1f} K)")
            return T, i
    try:
        T = brentq(lambda t: np.log(p_sat(t)) - lnP, T_MIN_ANT, T_MAX_ANT, xtol=1e-12)
    except ValueError:
        raise ValueError(f"P = {P_atm} atm hors du domaine de validite d'Antoine")
    return T, -1                          # -1 : solution obtenue par Brent

def t_boil_check(P_atm):
    """Verification independante (Brent)."""
    return brentq(lambda T: p_sat(T) - P_atm * ATM, T_MIN_ANT, T_MAX_ANT, xtol=1e-12)

# T d'ebullition normale calculee avec Antoine (coherence interne)
T_NB = brentq(lambda T: p_sat(T) - ATM, T_MIN_ANT, T_MAX_ANT, xtol=1e-12)

def dh_vap(T):
    """Chaleur latente [J/kg], relation de Watson."""
    if T >= TC:
        raise ValueError("T >= Tc : pas de chaleur latente")
    return DH_NB / (M * 1e-3) * ((1 - T / TC) / (1 - T_NB / TC)) ** 0.38


# ---------------- Modele du procede ----------------
def simulate(flow_kg_h=38.0, ratio_heater=0.60, P_atm=1.0, T_in_C=10.0):
    """
    flow_kg_h    : debit massique total d'acetone avant le splitter [kg/h]
    ratio_heater : fraction massique envoyee au Heater (0..1)
    P_atm        : pression du flux et du Heater [atm] (pas de perte de charge)
    T_in_C       : temperature du flux avant le splitter [C]
    """
    if not (0.0 <= ratio_heater <= 1.0):
        raise ValueError("ratio_heater doit etre entre 0 et 1")
    if flow_kg_h < 0 or P_atm <= 0:
        raise ValueError("debit >= 0 et pression > 0 requis")

    m_heater = flow_kg_h * ratio_heater
    m_bypass = flow_kg_h - m_heater
    Tb, n_it = t_boil(P_atm)
    T_in = T_in_C + 273.15

    m = m_heater / 3600.0                      # kg/s
    q_sens = m * CP_LIQ * (Tb - T_in)          # W
    q_lat  = m * dh_vap(Tb)                    # W
    return dict(flow=flow_kg_h, ratio=ratio_heater, P=P_atm, T_in_C=T_in_C,
                m_heater=m_heater, m_bypass=m_bypass,
                Tb_C=Tb - 273.15, n_iter=n_it,
                Q_sens_W=q_sens, Q_lat_W=q_lat, Q_W=q_sens + q_lat,
                flash_warning=bool(Tb < T_in))


def report(r):
    it = f"{r['n_iter']} it." if r["n_iter"] > 0 else "Brent"
    print("=" * 52)
    print(f" Entree : {r['flow']:.2f} kg/h | {r['T_in_C']:.1f} C | {r['P']:.3f} atm")
    print(f" Splitter : {r['ratio']*100:.1f} % Heater / {(1-r['ratio'])*100:.1f} % bypass")
    print(f"   -> Heater : {r['m_heater']:.3f} kg/h")
    print(f"   -> Bypass : {r['m_bypass']:.3f} kg/h")
    print("-" * 52)
    print(f" T ebullition ({it}) : {r['Tb_C']:.3f} C")
    print(f" Q sensible : {r['Q_sens_W']:9.1f} W")
    print(f" Q latente  : {r['Q_lat_W']:9.1f} W")
    print(f" Q TOTAL    : {r['Q_W']:9.1f} W = {r['Q_W']/1000:.3f} kW = {r['Q_W']*3.6:.0f} kJ/h")
    if r["flash_warning"]:
        print(" /!\\ T_entree > T_ebullition : le liquide serait en flash a cette pression;")
        print("     Q sensible negatif = difference d'enthalpie (fonction d'etat).")
    print("=" * 52)


# ---------------- Sensibilites ----------------
def sweep(base, name, values):
    """Fait varier un seul parametre, garde les autres a leur valeur de base.
    Les points hors domaine de validite sont renvoyes sous forme de NaN."""
    keymap = {"flow": "flow_kg_h", "ratio": "ratio_heater", "P": "P_atm", "T": "T_in_C"}
    out = []
    for v in values:
        kw = dict(base); kw[keymap[name]] = float(v)
        try:
            out.append(simulate(**kw))
        except ValueError:
            out.append(dict(Tb_C=np.nan, Q_W=np.nan))
    return out


def plot_all(base, Pmin=0.1, Pmax=3.0, npts=60, show=True):
    import matplotlib
    # Pas d'ecran (serveur, SSH, etc.) -> backend sans affichage
    if not show or (sys.platform.startswith("linux")
                    and not os.environ.get("DISPLAY")
                    and not os.environ.get("WAYLAND_DISPLAY")):
        matplotlib.use("Agg")
        show = False
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    # 1) Pression : Tb et Q
    P = np.linspace(Pmin, Pmax, npts)
    res = sweep(base, "P", P)
    a = ax[0, 0]
    l1, = a.plot(P, [r["Tb_C"] for r in res], "b-", label="T ebullition")
    l2 = a.axhline(base["T_in_C"], color="r", ls=":", label="T entree")
    a.set(xlabel="Pression (atm)", ylabel="T ebullition (C)",
          title=f"Tb et Q = f(P)   [{Pmin} -> {Pmax} atm]")
    a.grid(alpha=.3)
    a2 = a.twinx()
    l3, = a2.plot(P, [r["Q_W"] / 1000 for r in res], "g--", label="Q (kW)")
    a2.set_ylabel("Q (kW)", color="g")
    a.legend(handles=[l1, l2, l3], loc="best")

    # 2) Ratio vers le Heater
    x = np.linspace(0, 1, 51)
    res = sweep(base, "ratio", x)
    a = ax[0, 1]; a.plot(x * 100, [r["Q_W"] / 1000 for r in res], "g-")
    a.set(xlabel="% du flux vers le Heater", ylabel="Q (kW)",
          title="Q = f(ratio splitter)"); a.grid(alpha=.3)

    # 3) Debit total
    F = np.linspace(1, 2 * max(base["flow_kg_h"], 1), 50)
    res = sweep(base, "flow", F)
    a = ax[1, 0]; a.plot(F, [r["Q_W"] / 1000 for r in res], "g-")
    a.set(xlabel="Debit total (kg/h)", ylabel="Q (kW)", title="Q = f(debit)"); a.grid(alpha=.3)

    # 4) Temperature d'entree
    T = np.linspace(-20, 50, 50)
    res = sweep(base, "T", T)
    a = ax[1, 1]; a.plot(T, [r["Q_W"] / 1000 for r in res], "g-")
    a.set(xlabel="T entree (C)", ylabel="Q (kW)", title="Q = f(T entree)"); a.grid(alpha=.3)

    fig.suptitle(f"Sensibilites autour de : {base['flow_kg_h']} kg/h, ratio {base['ratio_heater']}, "
                 f"{base['P_atm']} atm, {base['T_in_C']} C")
    plt.tight_layout()
    plt.savefig("sensibilites_acetone.png", dpi=150)
    print(" Figure enregistree : sensibilites_acetone.png")
    if show:
        plt.show()
    plt.close(fig)


# ---------------- Saisie terminal ----------------
def ask(prompt, default, cast=float, check=None, err="valeur invalide"):
    """Demande une valeur ; Entree = defaut ; redemande tant que invalide."""
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
    print("=" * 52)
    print("  ACETONE : Splitter + Heater (vaporisation totale)")
    print("  (Entree = valeur par defaut entre crochets)")
    print("=" * 52)
    while True:
        print("\n-- Flux entrant --")
        flow  = ask("Debit massique total [kg/h]", 38.0, check=lambda v: v >= 0,
                    err="le debit doit etre >= 0")
        P     = ask("Pression [atm]", 1.0, check=lambda v: v > 0,
                    err="la pression doit etre > 0")
        T     = ask("Temperature d'entree [C]", 10.0, check=lambda v: -90 < v < 235,
                    err="temperature hors plage raisonnable (-90 a 235 C)")
        print("\n-- Splitter --")
        ratio = ask("Fraction envoyee au Heater (0 a 1)", 0.60,
                    check=lambda v: 0 <= v <= 1, err="le ratio doit etre entre 0 et 1")

        try:
            r = simulate(flow, ratio, P, T)
        except ValueError as e:
            print(f"\n Erreur : {e}\n Veuillez recommencer.")
            continue
        print()
        report(r)
        print(f" Verification Brent : Tb = {t_boil_check(P) - 273.15:.3f} C")

        if ask_yes_no("\n  Tracer le diagramme de sensibilite ?", True):
            print("\n-- Domaine de variation de la pression --")
            Pmin = ask("Pression min [atm]", 0.1, check=lambda v: v > 0,
                       err="Pmin doit etre > 0")
            Pmax = ask("Pression max [atm]", 3.0, check=lambda v: v > Pmin,
                       err=f"Pmax doit etre > Pmin ({Pmin})")
            npts = ask("Nombre de points", 60, cast=int, check=lambda v: v >= 2,
                       err="au moins 2 points")
            plot_all(dict(flow_kg_h=flow, ratio_heater=ratio, P_atm=P, T_in_C=T),
                     Pmin=Pmin, Pmax=Pmax, npts=npts)

        if not ask_yes_no("\n  Nouvelle simulation ?", False):
            print("\n Au revoir.")
            break


# ---------------- Main ----------------
def cli(argv):
    ap = argparse.ArgumentParser(description="Acetone : splitter + Heater")
    ap.add_argument("--flow",  type=float, default=38.0, help="debit total [kg/h]")
    ap.add_argument("--ratio", type=float, default=0.60, help="fraction vers le Heater (0-1)")
    ap.add_argument("--P",     type=float, default=1.0,  help="pression de base [atm]")
    ap.add_argument("--T",     type=float, default=10.0, help="temperature d'entree [C]")
    ap.add_argument("--Pmin",  type=float, default=0.1,  help="borne min du balayage en P [atm]")
    ap.add_argument("--Pmax",  type=float, default=3.0,  help="borne max du balayage en P [atm]")
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
        r = simulate(args.flow, args.ratio, args.P, args.T)
    except ValueError as e:
        raise SystemExit(f"Erreur : {e}")
    report(r)
    print(f" Verification Brent : Tb = {t_boil_check(args.P) - 273.15:.3f} C")
    print(f" Domaine P balaye   : {args.Pmin} -> {args.Pmax} atm ({args.npts} points)")
    if not args.no_plot:
        plot_all(dict(flow_kg_h=args.flow, ratio_heater=args.ratio,
                      P_atm=args.P, T_in_C=args.T),
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
