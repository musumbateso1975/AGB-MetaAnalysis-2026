# CHANGES.md — AGB Meta-Analysis R Project

## v_corrected — 2026-06-06 (Corrections post-audit)

### Données corrigées (`data/`)

**Doublons supprimés dans les trois CSV :**
| Study_ID supprimé | Raison | Doublon de |
|---|---|---|
| S6 | Observations identiques à S4 (Liu et al. 2024, Gabon) | S4 |
| S89 | Observations identiques à S88 (Wang et al. 2024) | S88 |
| S90 | Observations identiques à S88 (Wang et al. 2024) | S88 |
| S111 | 5 modèles identiques à S110 (Zhang et al. 2019, Conghua) | S110 |
| S101 | 27 observations identiques à S99 (Tian et al. 2024, Lin'an) | S99 |
| S308 | Identique à S304 (Migolet et al. 2022) — FisherZ/RMS only | S304 |
| S309 | Identique à S305 (Rodda et al. 2024) — FisherZ/RMS only | S305 |
| S310 | Identique à S306 (Wan et al. 2025) — FisherZ/RMS only | S306 |
| S311 | Identique à S307 (Lamulamu et al. 2022) — FisherZ/RMS only | S307 |

**Références disambiguées :**
| Study_ID | Ancienne référence | Nouvelle référence | Journal |
|---|---|---|---|
| S5 | Yang et al. 2025 | Yang et al. 2025a | IEEE JSTARS |
| S102 | Yang et al. 2025 | Yang et al. 2025b | J. Applied Remote Sensing |
| S123 | Yang et al. 2025 | Yang et al. 2025c | Ecological Indicators |

**Dimensions après correction :**
| Fichier | Avant | Après |
|---|---|---|
| DATA_FisherZ_latest.csv | 387 obs, 79 études | 314 obs, 70 études |
| DATA_RMS_latest.csv | 381 obs, 86 études | 308 obs, 77 études |
| AGB_master_database.csv | 384 obs, 82 études | 342 obs, 77 études |

---

### Scripts corrigés (`R/`)

#### `R/01_data_import.R`
- Commentaires d'en-tête mis à jour : counts n et études corrigés
- Commentaire doublons supprimés ajouté dans la documentation inline

#### `R/05_meta_regression.R`  ⚠️ CORRECTION CRITIQUE
- **Bug corrigé** : `sigma2.reml_total` n'est pas un slot valide de `rma.mv`
  (aurait provoqué une erreur R à l'exécution)
- **Nouvelle implémentation** : pseudo-R²_meta calculé correctement via
  comparaison avec un modèle null :
  ```r
  R2_meta = (tau2_null - tau2_model) / tau2_null
  ```
  Conforme à Viechtbauer (2010) et Nakagawa et al. (2022)
- Modèle null ajouté avant le modèle complet (section 3)

#### `R/09_sensitivity.R`
- Label S1 dynamique : `paste0("S1_Primary (n=", nrow(dat_fz), ")")`
  (était hardcodé `"S1_Primary (n=325)"`)

---

### Notes pour le manuscrit

Les chiffres suivants dans le texte doivent être mis à jour :
- Abstract : "75 studies" → 77 (RMSE) / 70 (FisherZ)
- Abstract : "387 Fisher's Z observations" → 314
- Abstract : "381 RMSE observations" → 308
- Section 2.1 : idem
- Table 1 : recalculer depuis les scripts après correction

