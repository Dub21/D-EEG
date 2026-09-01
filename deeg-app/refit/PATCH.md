# Réajustement : site en effet aléatoire

Quatre changements à `get_vector_region_optimized`. Les points 2 et 4 échouent
silencieusement si on les oublie : aucune erreur, des résultats faux.

---

## 1. Les formules

```python
m.gamlss_normative_model(
    mu=f'{roi_var} ~ ps(age) + as.factor(female_bin) + random(as.factor(Unique_Site_ID)) + ratio_ch_good',
    sigma='~ ps(age) + as.factor(female_bin) + random(as.factor(Unique_Site_ID)) + ratio_ch_good',
    method='mixed(50,50)',          # voir point 2
    train_col=[...],                 # inchangé
)
```

PyNM documente cette syntaxe : « If using `random()` to model a random effect in any
of the formulas, it must be passed a column of the dataframe with categorical values
as a factor ». C'est exactement la forme ci-dessus.

## 2. n.cyc n'est pas atteignable par l'API de PyNM

`GAMLSS.fit()` construit l'appel R en dur :

```python
r(f'gamlss({mu_f}, sigma.formula=..., family=..., data=train_data, method={method})')
```

Il n'y a pas de paramètre `control`, donc `gamlss.control(n.cyc=100)` ne peut pas
être passé. Le seul levier exposé est `method`, qui accepte `'mixed(n,m)'` :
RS pendant n itérations puis CG pendant m. **`method='mixed(50,50)'` porte le budget
à 100 itérations** et règle les 54 % de non-convergence sans toucher à PyNM.

Si tu préfères le contrôle explicite, l'autre voie est d'appeler gamlss directement
via rpy2 plutôt que de passer par PyNM.

## 3. Sauvegarder les constantes de standardisation

Elles n'ont jamais été écrites (cellule 38 : `del ... all_norm_params ...`).
Dans la boucle par marqueur, après l'appel :

```python
pd.DataFrame([{'marker': marker,
               'mean': all_norm_params[marker]['mean'],
               'std':  all_norm_params[marker]['std']}]) \
  .to_csv(f'{STATS_DIR}/norm_params.csv', mode='a', index=False,
          header=not os.path.exists(f'{STATS_DIR}/norm_params.csv'))
```

Sans ce fichier, personne hors du labo ne peut convertir des valeurs brutes.

## 4. `correct_site_effect` cesse de fonctionner, sans rien dire

Le code actuel :

```python
coefs = list(robjects.r('coef(model, what="mu")'))
site_coef = {n.split('Unique_Site_ID)')[-1]: v for n, v in coef_dict.items()
             if 'Unique_Site_ID' in n}
site_vals = df['Unique_Site_ID'].astype(str).map(site_coef).fillna(0.0).values
```

Avec `random()`, les effets de site **ne sont plus dans `coef()`** : ils vivent dans
l'objet du lisseur. Le dictionnaire sort vide, `.map()` produit des NaN, `.fillna(0.0)`
les convertit en zéros, et `total_effect` vaut zéro partout. La correction de site
est donc désactivée en silence, sans exception ni avertissement, et les trajectoires
sortent non corrigées alors que le code prétend le contraire.

Remplacement :

```python
robjects.globalenv['model'] = m.gamlss_model.model
ro_coef = robjects.r('''
    sm <- getSmo(model, what="mu")
    setNames(as.numeric(coef(sm)), names(coef(sm)))
''')
site_coef = dict(zip(list(ro_coef.names), list(ro_coef)))
```

À vérifier après le premier modèle : `len(site_coef)` doit valoir 14, pas 0.

---

## Deux contrôles à faire sur le premier modèle réajusté

1. **Les z-scores des contrôles restent ~N(0,1).** `predict.gamlss` avec `newdata`
   sur des termes `random()` est fragile. Si un site apparaît chez les patients mais
   pas chez les contrôles, la prédiction échouera pour ces sujets : le modèle est
   entraîné sur les contrôles seuls, alors qu'avec `as.factor` les niveaux venaient
   du tableau complet.
2. **Comparer les z-scores fixe contre aléatoire sur un marqueur.** Le rétrécissement
   modifie les valeurs. Avec 14 sites bien peuplés, l'écart devrait être faible, mais
   il n'est pas nul : les z-scores et les tailles d'effet de l'article vont bouger.
   Mieux vaut le mesurer maintenant que le découvrir en révision.

## Pourquoi ce changement, en une phrase pour les Methods

Un effet fixe de site n'est pas extrapolable à un site jamais observé, ce qui rend
le modèle inutilisable hors des cohortes d'origine ; l'effet aléatoire rétrécit vers
la moyenne de population et donne une prédiction définie pour un nouveau site.
Mesure de l'enjeu : l'étendue médiane des offsets de site est de 1,02 écart-type,
et dépasse 1 écart-type pour 53 % des modèles.
