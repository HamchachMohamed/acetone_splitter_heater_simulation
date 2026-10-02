Acetone Splitter + Heater

## Description

Simulation d'un procédé de **vaporisation totale d'acétone** : un flux liquide est divisé par un splitter, une fraction est envoyée dans un **Heater** où elle est entièrement vaporisée, le reste part en bypass.

Le programme calcule :
- la **température d'ébullition** Tb(P) via l'équation d'Antoine (résolution Newton-Raphson, repli sur Brent),
- la **chaleur latente** ΔH_vap(Tb) par la relation de Watson,
- la **puissance thermique** du Heater : Q = ṁ·[cp_liq·(Tb − T₀) + ΔH_vap(Tb)],
- et trace un **diagramme de sensibilité** sur 4 paramètres.

---

## Fichier

```
acetone_heater_param.py
```

Dépendances : `numpy`, `scipy`, `matplotlib`.

```bash
pip install numpy scipy matplotlib
```

---

## Utilisation

Le script fonctionne en **deux modes**.

### 1. Mode interactif (par défaut)

```bash
python acetone_heater_param.py
```

Le programme pose les questions dans le terminal. **Appuyer sur Entrée** accepte la valeur par défaut affichée entre crochets.

Exemple de session :

```
-- Flux entrant --
  Debit massique total [kg/h] [38.0] :
  Pression [atm] [1.0] :
  Temperature d'entree [C] [10.0] :

-- Splitter --
  Fraction envoyee au Heater (0 a 1) [0.6] :

  Tracer le diagramme de sensibilite ? [O/n] : o

-- Domaine de variation de la pression --
  Pression min [atm] [0.1] :
  Pression max [atm] [3.0] :
  Nombre de points [60] :
```

Puis possibilité d'enchaîner sur une **nouvelle simulation** sans quitter le programme.

### 2. Mode non interactif (arguments CLI)

```bash
python acetone_heater_param.py --flow 50 --ratio 0.7 --P 2 --T 25 --Pmin 0.2 --Pmax 5
python acetone_heater_param.py --flow 38 --ratio 0.6 --P 1 --T 10 --no-plot
```

---

## Arguments CLI

| Argument | Type | Défaut | Description |
|---|---|---|---|
| `--flow` | float | `38.0` | Débit massique total avant splitter [kg/h] |
| `--ratio` | float | `0.60` | Fraction massique envoyée au Heater (0 à 1) |
| `--P` | float | `1.0` | Pression de base du procédé [atm] |
| `--T` | float | `10.0` | Température d'entrée [°C] |
| `--Pmin` | float | `0.1` | Borne **inférieure** du balayage en pression [atm] |
| `--Pmax` | float | `3.0` | Borne **supérieure** du balayage en pression [atm] |
| `--npts` | int | `60` | Nombre de points pour la courbe Tb = f(P) |
| `--no-plot` | flag | — | Ne pas tracer les courbes (mode headless) |

### Validations

- `--flow >= 0`
- `0 <= --ratio <= 1`
- `--P > 0`, `--Pmin > 0`, `--Pmax > 0`
- `--Pmin < --Pmax`
- `--npts >= 2`

---

## Sortie console

```
====================================================
 Entree : 38.00 kg/h | 10.0 C | 1.000 atm
 Splitter : 60.0 % Heater / 40.0 % bypass
   -> Heater : 22.800 kg/h
   -> Bypass : 15.200 kg/h
----------------------------------------------------
 T ebullition (5 it.) : 56.095 C
 Q sensible :    6244.6 W
 Q latente  :   15970.3 W
 Q TOTAL    :   22214.9 W = 22.215 kW = 79974 kJ/h
====================================================
 Verification Brent : Tb = 56.095 C
 Domaine P balaye   : 0.1 -> 3.0 atm (60 points)
 Figure enregistree : sensibilites_acetone.png
```

Si `T_entrée > T_ébullition`, un **avertissement flash** est affiché : le liquide serait déjà partiellement vaporisé à cette pression, et Q sensible devient négatif (différence d'enthalpie).

---

## Diagramme de sensibilité

Une figure **2×2** est générée et sauvegardée dans `sensibilites_acetone.png` (répertoire courant).

| Sous-plot | Courbe | Domaine |
|---|---|---|
| Haut-gauche | **Tb = f(P)** (bleu) + **Q = f(P)** (vert, axe droit) | `Pmin → Pmax`, `npts` points |
| Haut-droit | **Q = f(ratio splitter)** | ratio 0 → 1 |
| Bas-gauche | **Q = f(débit total)** | 1 → 2×débit base |
| Bas-droit | **Q = f(T entrée)** | −20 → 50 °C |

Sur le sous-plot haut-gauche, la température d'entrée est tracée en **pointillé rouge** comme référence.

---

## Modèle physique

### Équation d'Antoine (NIST)

$$\log_{10} P_{sat}[\text{bar}] = A - \frac{B}{T + C}$$

avec `A = 4.42448`, `B = 1312.253`, `C = −32.445`, valable **259 K à 507 K**.

### Résolution de Tb

1. **Newton-Raphson** sur `f(T) = ln Psat(T) − ln P`
2. **Repli sur Brent** si Newton sort du domaine ou diverge
3. **Vérification indépendante** par `brentq` à la fin

### Chaleur latente (Watson)

$$\Delta H_{vap}(T) = \frac{\Delta H_{nb}}{M} \left[\frac{1 - T/T_c}{1 - T_{nb}/T_c}\right]^{0.38}$$

avec `ΔH_nb = 29.1 kJ/mol`, `Tc = 508.1 K`, `M = 58.08 g/mol`.

### Puissance

$$Q = \dot{m} \left[ c_{p,liq} (T_b - T_0) + \Delta H_{vap}(T_b) \right]$$

avec `cp_liq = 2.17 kJ/(kg·K)`, `ṁ` en kg/s.

---

## Données physiques (acétone)

| Constante | Valeur | Unité |
|---|---|---|
| Masse molaire M | 58.08 | g/mol |
| Température critique Tc | 508.1 | K |
| ΔH_vap à T_nb | 29.1 | kJ/mol |
| cp liquide | 2.17 | kJ/(kg·K) |
| Antoine A | 4.42448 | — |
| Antoine B | 1312.253 | K |
| Antoine C | −32.445 | K |
| Domaine Antoine | 259 – 507 | K |

---

## Comportement headless

Sur un serveur sans écran (SSH, CI, Docker), le programme bascule automatiquement en **backend `Agg`** de matplotlib : la figure est sauvegardée mais **pas affichée**. Détection via `DISPLAY` / `WAYLAND_DISPLAY`.

---

## Exemples

```bash
# Cas de base, mode interactif
python acetone_heater_param.py

# Cas rapide sans graphique
python acetone_heater_param.py --flow 38 --ratio 0.6 --P 1 --T 10 --no-plot

# Balayage étroit en pression, 100 points
python acetone_heater_param.py --Pmin 0.5 --Pmax 2.0 --npts 100

# Domaine large
python acetone_heater_param.py --Pmin 0.05 --Pmax 10

# Combinaison complète
python acetone_heater_param.py --flow 50 --ratio 0.7 --P 1.5 --T 20 \
                               --Pmin 0.3 --Pmax 4.0 --npts 80
```

---

## Notes

- **Pas de perte de charge** : la pression est supposée identique en amont et dans le Heater.
- **Vaporisation totale** : le flux envoyé au Heater ressort entièrement vapeur saturée à Tb(P).
- **Bypass** : le flux non envoyé au Heater reste liquide à T₀ ; aucun mélange aval n'est modélisé ici.
- Les points hors domaine de validité d'Antoine sont renvoyés en `NaN` dans le balayage (pas de crash du graphe).
