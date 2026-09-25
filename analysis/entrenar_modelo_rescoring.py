# -*- coding: utf-8 -*-
"""Entrenamiento reproducible de predictor MM-GBSA.

Uso: python entrenar_modelo_rescoring.py
Entrena ExtraTrees y RandomForest con Vina + descriptores RDKit.
La validacion es GroupKFold por ligand para evitar fuga entre targets.
"""
from pathlib import Path
import json, warnings
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.model_selection import GroupKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.metrics import make_scorer, mean_absolute_error, r2_score
import joblib

warnings.filterwarnings('ignore')
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'rescoring_datos'; OUT=ROOT/'analysis'/'modelo_rescoring'
OUT.mkdir(exist_ok=True)
FILES=['rescoring_mmgbsa_v5.csv','rescoring_nuevo.csv','rescoring_tdp43_top100.csv','rescoring_tdp43_top100_v2.csv','rescoring_tdp43_mmgbsa.csv','rescoring_mmgbsa.csv']
SMILES=['candidatos_limpios.csv','candidatos_mmgbsa.csv']

def descriptors(s):
 m=Chem.MolFromSmiles(str(s))
 if m is None: return None
 return {'MolWt':Descriptors.MolWt(m),'LogP':Crippen.MolLogP(m),'TPSA':rdMolDescriptors.CalcTPSA(m),'HBD':Lipinski.NumHDonors(m),'HBA':Lipinski.NumHAcceptors(m),'RotB':Lipinski.NumRotatableBonds(m),'Rings':Lipinski.RingCount(m),'AromRings':rdMolDescriptors.CalcNumAromaticRings(m),'HeavyAtom':Lipinski.HeavyAtomCount(m),'FractionCSP3':rdMolDescriptors.CalcFractionCSP3(m)}

parts=[]
for f in FILES:
 p=DATA/f
 if p.exists():
  d=pd.read_csv(p); d.columns=[c.strip() for c in d.columns]
  if 'vina_affinity' not in d.columns and 'affinity' in d.columns:
   d=d.rename(columns={'affinity':'vina_affinity'})
  need={'ligand','target','vina_affinity','mmgbsa_dG'}
  if need.issubset(d.columns): parts.append(d[list(need)].assign(source=f))
if not parts: raise SystemExit('No hay archivos de rescoring validos')
raw=pd.concat(parts,ignore_index=True)
raw['vina_affinity']=pd.to_numeric(raw['vina_affinity'],errors='coerce'); raw['mmgbsa_dG']=pd.to_numeric(raw['mmgbsa_dG'],errors='coerce')
raw=raw.dropna(subset=['ligand','target','vina_affinity','mmgbsa_dG'])
raw=raw[raw.mmgbsa_dG.abs()>1e-9].copy()
n_poses=len(raw)
raw=raw.sort_values('mmgbsa_dG').drop_duplicates(['ligand','target'],keep='first')  # mejor pose por par
raw=raw[raw.mmgbsa_dG.abs()>1e-9].copy()
# smiles desde ambas tablas, normalizando columnas
maps={}
for f in SMILES:
 p=DATA/f
 if not p.exists(): p=ROOT/'analysis'/f
 if p.exists():
  d=pd.read_csv(p); d.columns=[str(c).strip().lower() for c in d.columns]
  if 'smiles' in d.columns:
   idc=next((c for c in ['ligand','molecule_chembl_id','chembl_id','compound_id','id'] if c in d.columns),None)
   if idc:
    for _,r in d[[idc,'smiles']].dropna().iterrows(): maps.setdefault(str(r[idc]).strip(),r['smiles'])
# SMILES extra de ChEMBL para ligandos sin cobertura local
p=DATA/'chembl_smiles_extra.csv'
if p.exists():
 d=pd.read_csv(p); d.columns=[str(c).strip().lower() for c in d.columns]
 if {'ligand','smiles'}.issubset(d.columns):
  for _,r in d[['ligand','smiles']].dropna().iterrows(): maps.setdefault(str(r['ligand']).strip(),r['smiles'])
raw['smiles']=raw.ligand.astype(str).str.strip().map(maps)
# fallback: remove blanks; no imputacion estructural
raw=raw.dropna(subset=['smiles']).copy()
raw['desc']=raw.smiles.map(descriptors); raw=raw[raw.desc.notna()].copy()
D=pd.DataFrame(raw.pop('desc').tolist(),index=raw.index)
X=pd.concat([raw[['vina_affinity','target']].reset_index(drop=True),D.reset_index(drop=True)],axis=1)
y=raw.mmgbsa_dG.reset_index(drop=True); groups=raw.ligand.reset_index(drop=True)
num=['vina_affinity']+list(D.columns); cat=['target']
pre=ColumnTransformer([('num',SimpleImputer(strategy='median'),num),('cat',OneHotEncoder(handle_unknown='ignore'),cat)])
models={'ExtraTrees':ExtraTreesRegressor(n_estimators=400,min_samples_leaf=2,max_features=0.8,random_state=42,n_jobs=-1),'RandomForest':RandomForestRegressor(n_estimators=400,min_samples_leaf=2,max_features=0.8,random_state=42,n_jobs=-1)}
n_splits=min(5,groups.nunique()); cv=GroupKFold(n_splits=n_splits)
sc={'MAE':make_scorer(mean_absolute_error,greater_is_better=False),'R2':make_scorer(r2_score)}
report={'dataset':{'raw_rows':int(sum(len(x) for x in parts)),'poses_validas':n_poses,'poses_duplicadas_eliminadas':n_poses-len(raw),'usable_rows':len(X),'unique_ligands':int(groups.nunique()),'targets':raw.target.value_counts().to_dict(),'excluded_zero_dG':int(sum(pd.to_numeric(pd.concat(parts).mmgbsa_dG,errors='coerce').fillna(0).abs()<1e-9))},'features':num+cat,'models':{}}
best=None
for name,model in models.items():
 pipe=Pipeline([('pre',pre),('model',model)])
 z=cross_validate(pipe,X,y,cv=cv,groups=groups,scoring=sc,n_jobs=1)
 vals={'MAE_mean':float(-z['test_MAE'].mean()),'MAE_std':float(z['test_MAE'].std()),'R2_mean':float(z['test_R2'].mean()),'R2_std':float(z['test_R2'].std()),'fold_MAE':(-z['test_MAE']).tolist(),'fold_R2':z['test_R2'].tolist()}
 report['models'][name]=vals
 if best is None or vals['MAE_mean']<best[1]['MAE_mean']: best=(name,vals)
name,_=best
final=Pipeline([('pre',pre),('model',models[name])]).fit(X,y)
joblib.dump(final,OUT/'modelo_mmgbsa.joblib')
raw.to_csv(OUT/'dataset_entrenamiento.csv',index=False)
report['selected_model']=name; report['selected_metrics']=report['models'][name]; report['caveat']='Prediccion exploratoria: dataset pequeno y no sustituye MM-GBSA ni validacion experimental.'
(OUT/'reporte_entrenamiento.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
(OUT/'REPORTE_ENTRENAMIENTO.md').write_text('# Modelo de rescoring MM-GBSA\n\n'+json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
print(json.dumps(report,indent=2,ensure_ascii=False))
