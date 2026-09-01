# Étape 1 — extraire les paramètres normatifs sur le cluster

## Pourquoi ne pas rapatrier les `.rds` directement

Un objet `gamlss` embarque les données d'entraînement (`$y`, model frame, matrices de
design). 2292 objets × ~4800 sujets, c'est à la fois plusieurs Go à transférer et des
données de participants qui finiraient sur shinyapps.io. En plus, `predict.gamlss` sur
des termes `pb()` a besoin des données d'origine pour reconstruire la base de splines :
un `.rds` "allégé" ne prédirait pas correctement de toute façon.

## Ce qu'on exporte à la place

Pour chaque modèle, les quatre paramètres SHASH (μ, σ, ν, τ) évalués sur une grille
âge × sexe, au site de référence. C'est agrégé, ça ne contient aucune donnée
individuelle, et ça suffit pour tout ce que fait BrainChart :

- centile d'un nouveau sujet = `pSHASHo2(valeur, mu, sigma, nu, tau)` à son âge
- z-score = `qNO(centile)`
- courbes de trajectoire = quantiles de la même distribution le long de la grille

## Exécution

```bash
module load r/4.4     # ou l'équivalent sur ta grappe
Rscript export_normative_params.R /chemin/vers/les/modeles /chemin/de/sortie 300
```

Vérifier d'abord les noms de colonnes des covariables en haut du script
(`COL_AGE`, `COL_SEX`, `COL_SITE`, `COL_QC`) : ils doivent correspondre à ceux
utilisés à l'ajustement.

## Ce que je récupère ensuite

Deux fichiers, quelques dizaines de Mo :

- `normative_params.parquet` (ou `.csv.gz`) — colonnes : `marker, age, sex, mu, sigma, nu, tau`
- `model_manifest.csv` — famille de distribution, n, plage d'âge, site de référence,
  offsets de site par modèle

Transfert :

```bash
rsync -avz user@cluster:/chemin/de/sortie/ ~/D-EEG/deeg-app/data/
```

## Point de contrôle avant de construire l'app

Sur le cluster, prendre ~50 sujets d'entraînement, recalculer leur z-score via la table
exportée, et comparer aux résidus du modèle. Si la corrélation n'est pas ~1.0, la grille
d'âge est trop grossière ou une covariable manque. Ce test reste sur le cluster : il
utilise des données individuelles.
