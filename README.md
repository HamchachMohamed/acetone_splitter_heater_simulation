## Description

Simulation d’un procédé **acétone → splitter → heater** avec vaporisation totale
et options de surchauffe, perte de charge et pompe.

Le script calcule la puissance thermique nécessaire au heater, la température
d’ébullition, les contributions sensible / latente / surchauffe, et trace des
courbes de sensibilité.

---

## Fonctionnalités

- Cp liquide dépendant de la température (Rowlinson-Bondi), intégré entre
  `T_in` et `Tb`.
- Zone de flash masquée (`NaN`) et ombrée dans les graphiques.
- Surchauffe optionnelle de la vapeur.
- Perte de charge du heater optionnelle.
- Travail de pompe optionnel (indicatif, **non inclus dans Q**).
- Validation `np.isfinite` de toutes les entrées : `nan` / `inf` refusés.
- Mode interactif ou ligne de commande.
- Génération d’une figure de sensibilité : `sensibilites_acetone.png`.

---

## Installation

Dépendances Python :

```bash
pip install numpy scipy matplotlib
Python 3.8+ recommandé.

Utilisation
Mode interactif (par défaut)
bash
python acetone_heater_param_v2.py
Le programme pose une série de questions. Une valeur par défaut est proposée
entre crochets ; appuyez sur Entrée pour l’accepter.

Mode non interactif
bash
python acetone_heater_param_v2.py --flow 50 --ratio 0.7 --P 2 --T 25 --superheat 10 --dP 0.1
bash
python acetone_heater_param_v2.py --flow 38 --ratio 0.6 --P 3 --Psource 1 --no-plot
Arguments en ligne de commande
Argument	Défaut	Description
--flow	38.0	Débit massique total avant splitter [kg/h]
--ratio	0.60	Fraction massique envoyée au heater (0 à 1)
--P	1.0	Pression d’entrée du heater [atm]
--T	10.0	Température d’entrée [°C]
--superheat	0.0	Surchauffe de la vapeur au-dessus de la saturation [K]
--dP	0.0	Perte de charge du heater [atm]
--Psource	None	Pression amont de la pompe [atm]. None = pas de pompe
--Pmin	0.1	Borne min du balayage en pression [atm]
--Pmax	3.0	Borne max du balayage en pression [atm]
--npts	60	Nombre de points pour Tb = f(P)
--no-plot	—	Ne pas tracer / enregistrer les courbes
Paramètres interactifs
Le mode interactif demande :

Débit massique total [kg/h]

Pression [atm]

Température d’entrée [°C]

Fraction envoyée au heater (0 à 1)

Surchauffe de la vapeur [K]

Perte de charge du heater [atm]

Pression amont de la pompe [atm]

Puis il affiche le bilan et propose de tracer le diagramme de sensibilité.

Modèle physique
Propriétés de l’acétone
Masse molaire : M = 58.08 g/mol

Température critique : TC = 508.1 K

Facteur acentrique : OMEGA = 0.307

Chaleur latente à T_NB : DH_NB = 29.1e3 J/mol

Coefficient Antoine : A = 4.42448, B = 1312.253, C = -32.445

Domaine Antoine : 259.0 K à 507.0 K

Masse volumique liquide pour la pompe : RHO_LIQ = 790 kg/m³

Point de fusion : T_FREEZE = 178.5 K

Cp liquide de référence à 298.15 K : 2.17 kJ/(kg.K)

Corrélations utilisées
Pression de saturation : Antoine
log10(Psat[bar]) = A - B / (T + C)

Température d’ébullition : Newton-Raphson sur ln(Psat) - ln(P),
avec repli sur Brent.

Chaleur latente : relation de Watson
dh_vap(T) = DH_NB / M * ((1 - T/TC) / (1 - T_NB/TC))^0.38

Cp gaz parfait : Shomate, valable de 298 K à 1200 K.

Cp liquide : Rowlinson-Bondi, recalé pour imposer
Cp(298.15 K) = 2.17 kJ/(kg.K).

Bilan thermique
Pour le débit envoyé au heater :

text
Q_sensible = m * ∫ Cp_liq(T) dT   de T_in à Tb
Q_latente  = m * dh_vap(Tb)
Q_surchauffe = m * ∫ Cp_gas(T) dT de Tb à T_out   si surchauffe > 0
Q_total = Q_sensible + Q_latente + Q_surchauffe
Le travail de pompe est calculé séparément :

text
W_pompe = m / RHO_LIQ * (P - P_source) * ATM * 1e5 / eta_pompe
Il est indicatif et n’est pas ajouté à Q_total.

Perte de charge
La vaporisation est évaluée à la pression moyenne : P_vap = P - dP/2.

La sortie heater est à : P_out = P - dP.

Zone de flash
Si T_in > Tb, le liquide serait en flash à cette pression.
Le modèle suppose un liquide sous-refroidi : Q n’est alors pas valide.
Dans les graphiques, cette zone est masquée (NaN) et ombrée en rouge.

Sorties
Rapport console
Le script affiche :

Débit total, température et pression d’entrée

Répartition splitter : heater / bypass

Température d’ébullition

Cp liquide moyen

Q_sensible, Q_latente, Q_surchauffe, Q_total

Température et pression de sortie heater

Travail de pompe indicatif si applicable

Vérification Brent de Tb

Figure
Si le tracé est activé, une figure est enregistrée dans le répertoire courant :

text
sensibilites_acetone.png
Elle contient 4 sous-graphiques :

Tb et Q = f(P)

Q = f(ratio splitter)

Q = f(débit total)

Q = f(T entrée)

Limites et précautions
Modèle stationnaire.

Vaporisation totale supposée dans le heater.

Le bypass n’est pas mélangé à la vapeur en aval : le modèle calcule
uniquement le heater.

Cp gaz parfait : Shomate valable 298–1200 K.

Cp liquide : Rowlinson-Bondi recalé, domaine de validité limité.

La pompe est indicative et hors bilan thermique du heater.

Les entrées nan / inf sont refusées.

En zone de flash, Q ne doit pas être utilisé.

Licence
Aucune licence n’est spécifiée dans le code source.
