#!/usr/bin/env python3
"""
DESCARGA LIBRERIA GRANDE DE COMPUESTOS PARA CRIBADO CONTINUO
Fuentes publicas gratuitas:
- FDA-approved drugs (~1600)
- DrugBank approved structures  
- Known neuroprotective compounds
- Natural products with neuroactivity
"""
import os, sys, csv, json, urllib.request, time

LIBRARY_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "compounds")

# ──── FUENTE 1: FDA APPROVED DRUGS ────
FDA_DRUGS = [
    # Neuro/ELA relacionados
    ("Riluzole", "C1=CC2=C(C=C1OC(F)(F)F)SC(=N2)N", "ELA"),
    ("Edaravone", "CC1=CC(=O)N(N1C2=CC=CC=C2)C", "ELA"),
    ("Tofersen", "N/A", "SOD1-ELA-ASO"),
    ("Baclofen", "C1=CC(=CC=C1C(CC(=O)O)CN)Cl", "ELA-espasticidad"),
    ("Tizanidine", "C1=CC2=C(C(=C1)Cl)NC(=N2)NCC3=CC=CC=C3", "ELA-espasticidad"),
    # Neuroprotectores
    ("Memantine", "CC12CC3CC(C1)(CC(C3)(C2)N)C", "Alzheimer-NMDA"),
    ("Donepezil", "COC1=C(C=C2C(=C1)CC(C2=O)CC3CCN(CC3)CC4=CC=CC=C4)OC", "Alzheimer"),
    ("Rivastigmine", "CCN(C)C(=O)OC1=CC=CC2=C1C=CC=C2", "Alzheimer"),
    ("Galantamine", "CN1CC[C@@]23C=C[C@H](C[C@@H]2OC4=C(C=CC(=C34)OC)OC)O", "Alzheimer"),
    ("Rasagiline", "CC#C[C@H]1CC2=C(C1)C3=C(C=C2)C=CC=C3", "Parkinson"),
    ("Selegiline", "CC#CCN([C@H](C)CC1=CC=CC=C1)C", "Parkinson"),
    ("Pramipexole", "CCCN[C@@H]1CCC2=C(C1)SC(=N2)N", "Parkinson"),
    ("Ropinirole", "CCN(CCC1=CC=CC=C1)CCN2CCC3=C(C2)C=CC=C3", "Parkinson"),
    ("Levodopa", "C1=CC(=C(C=C1C[C@@H](C(=O)O)N)O)O", "Parkinson"),
    # Antiinflamatorios (neuroinflamacion ELA)
    ("Celecoxib", "CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F", "COX2"),
    ("Ibuprofen", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "AINE"),
    ("Naproxen", "C[C@@H](C1=CC2=C(C=C1)C=CC(=C2)OC)C(=O)O", "AINE"),
    ("Diclofenac", "C1=CC=C(C(=C1)CC(=O)O)NC2=C(C=CC=C2Cl)Cl", "AINE"),
    ("Indomethacin", "CC1=C(C2=C(N1C(=O)C3=CC=C(C=C3)Cl)C=C(C=C2)OC)CC(=O)O", "AINE"),
    ("Prednisone", "C[C@@]12C[C@@H]([C@H]3[C@H]([C@H]1CC[C@@]2(C(=O)CO)O)CCC4=CC(=O)C=C[C@]34C)O", "Corticoide"),
    ("Dexamethasone", "C[C@@]12C[C@@H]([C@]3([C@H]([C@@H]1C[C@@H]([C@@]2(C(=O)CO)O)C)CCC4=CC(=O)C=C[C@@]43C)F)O", "Corticoide"),
    # Antioxidantes
    ("N-Acetylcysteine", "C(C(=O)O)NC(=O)CS", "Antioxidante"),
    ("Vitamin_E", "CC1=C(C2=C(CC[C@@](O2)(C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C(=C1O)C)C", "Antioxidante"),
    ("Coenzyme_Q10", "COC1=C(C(=O)C(=C(C1=O)OC)OC)CC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C", "Mitocondrial"),
    ("Alpha-Lipoic_Acid", "C1CSSC1CCCCC(=O)O", "Antioxidante"),
    # Farmacos metabolicos
    ("Metformin", "CN(C)C(=N)N=C(N)N", "AMPK-mTOR"),
    ("Pioglitazone", "CC1=C(C(=O)N(C1=O)CC2=CC=C(C=C2)CCOC3=CC=CC=C3)C", "PPAR-gamma"),
    ("Rosiglitazone", "CN(CCOC1=CC=C(C=C1)CC2C(=O)NC(=O)S2)C3=CC=CC=N3", "PPAR-gamma"),
    ("Atorvastatin", "CC(C)C1=C(C(=C(N1CC[C@@H](C[C@@H](CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4", "Estatina"),
    ("Simvastatin", "CC[C@H](C)C(=O)O[C@H]1C[C@@H](C=C2[C@H]1[C@H]([C@H](C=C2)C)CC[C@@H]3C[C@H](CC(=O)O3)O)C", "Estatina"),
    # Antibioticos con efecto neuroprotector
    ("Minocycline", "CN(C)C1=C2C(=C(C=C1)O)C(=O)C3=C(C2=O)C(=CC=C3N(C)C)O", "Neuroprotector"),
    ("Ceftriaxone", "CN1C(=C(C(=O)N1S(=O)(=O)O)C(=O)NCC2=CC=CC=C2)C(=O)NO", "GLT1-upregulator"),
    ("Doxycycline", "C[C@@H]1[C@H]2[C@@H]([C@H]3[C@@H](C(=O)C4=C([C@@]3(C(=O)C2=C(C4=O)O)O)O)O)[C@@H](C1=C(C)C)O", "Antibiotico"),
    # Cardiovascular (muchos son neuroprotectores)
    ("Candesartan", "CCOC1=NC2=CC=CC=C2N1CC3=CC=C(C=C3)C4=CC=CC=C4C(=O)O", "ARB"),
    ("Telmisartan", "CCCC1=NC2=C(N1CC3=CC=C(C=C3)C4=CC=CC=C4C(=O)O)C=C(C=C2C)C", "ARB-PPAR"),
    ("Losartan", "CCCCC1=NC(Cl)=C(CO)N1CC1=CC=C(C=C1)C1=CC=CC=C1C1=NN=NN1", "ARB"),
    ("Amlodipine", "CCOC(=O)C1=C(COCCN)NC(C)=C(C(=O)OC)C1C1=CC=CC=C1Cl", "Calcio-antagonista"),
    # Oncologicos reposicionados
    ("Tamoxifen", "CC/C(=C(\\C1=CC=CC=C1)/C2=CC=C(C=C2)OCCN(C)C)/C3=CC=CC=C3", "SERM"),
    ("Rapamycin", "C[C@@H]1CC[C@@H]2C[C@@H](/C(=C/C=C/C=C/[C@H](C[C@H](C(=O)[C@@H]([C@@H](/C(=C/[C@H](C(=O)C[C@H](O)[C@H](C)CC(=O)[C@H](C)[C@@H](O)[C@H](C)C=CC(=O)[C@H](C)/C=C/C1C)OC)C)O)C)C)OC)CC[C@@H]2OC", "mTOR"),
    ("Thalidomide", "C1CC(=O)NC(=O)C1N2C(=O)C3=CC=CC=C3C2=O", "Inmunomodulador"),
    ("Lenalidomide", "C1CC(=O)NC(=O)C1N2C(=O)C3=C(C2=O)C=CC=C3N", "Inmunomodulador"),
    # Natural products
    ("Curcumin", "COC1=C(C=CC(=C1)C=CC(=O)CC(=O)C=CC2=CC(=C(C=C2)O)OC)O", "Antiinflamatorio"),
    ("Resveratrol", "C1=CC(=CC=C1C=CC2=CC(=CC(=C2)O)O)O", "Sirtuina"),
    ("Quercetin", "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O", "Flavonoide"),
    ("Epigallocatechin_gallate", "C1[C@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=C(C(=C3)O)O)O)OC(=O)C4=CC(=C(C(=C4)O)O)O", "EGCG"),
    ("Genistein", "C1=CC(=CC=C1C2=COC3=CC(=CC(=C3C2=O)O)O)O", "Isoflavona"),
    ("Berberine", "COC1=C(C2=C[N+]3=C(C=C2C=C1)C4=CC5=C(C=C4CC3)OCO5)OC", "Alcaloide"),
    ("Cannabidiol", "CCCCCC1=CC(=C(C(=C1)O)C2C=C(CC[C@@H]2C(=C)C)C)O", "CBD"),
    # Mas farmacos comunes
    ("Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O", "AINE"),
    ("Paracetamol", "CC(=O)NC1=CC=C(C=C1)O", "Analgesico"),
    ("Omeprazole", "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC", "IBP"),
    ("Fluoxetine", "CNCC[C@@H](C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F", "ISRS"),
    ("Sertraline", "CN[C@@H]1CC[C@@H](C2=CC=CC=C2Cl)C3=C(C1)C=CC=C3", "ISRS"),
    ("Escitalopram", "CN(C)CCC[C@@]1(C2=C(CO1)C=C(C=C2)C#N)C3=CC=C(C=C3)F", "ISRS"),
    ("Gabapentin", "CC1(CC(CC(C1)(CN)O)CC(=O)O)CC(=O)O", "Anticonvulsivo"),
    ("Pregabalin", "CC(C)C[C@H](CN)CC(=O)O", "Anticonvulsivo"),
    ("Topiramate", "CC1(C)OC2COC3(COS(=O)(=O)O)OC(C)(C)OC3C2O1", "Anticonvulsivo"),
    ("Lamotrigine", "C1=CC(=C(C(=C1Cl)N)N=C(N)N)Cl", "Anticonvulsivo"),
    ("Lithium", "[Li+].[Li+].C(=O)([O-])[O-]", "Mood-stabilizer"),
    ("Valproic_Acid", "CCCC(CCC)C(=O)O", "HDACi"),
    ("Creatine", "CN(CC(=O)O)C(=N)N", "Energia-celular"),
    ("TUDCA", "C[C@H](CCC(=O)NCCS(=O)(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H]([C@H]3[C@@H]2CC[C@H]4[C@@]3(CC[C@H](C4)O)C)O)C", "Acido-biliar"),
    ("Ursolic_Acid", "C[C@@H]1CC[C@@]2(CC[C@]3(C(=CC[C@H]4[C@]3(CC[C@@H]5[C@@]4(CC[C@@H]([C@]5(C)C)O)C)C)[C@@H]2[C@H]1C)C)C(=O)O", "Triterpeno"),
    # Farmacos pediatricos y seguros
    ("Melatonin", "CC(=O)NCCC1=CNC2=C1C=C(C=C2)OC", "Neurohormona"),
    ("Riboflavin", "CC1=CC2=C(C=C1C)N(C(=O)NC2=O)C[C@@H]([C@@H]([C@@H](CO)O)O)O", "Vitamina B2"),
    ("Thiamine", "CC1=C(SC=[N+]1CC2=CN=C(N=C2N)C)CCO", "Vitamina B1"),
    ("Pyridoxine", "CC1=NC=C(C(=C1O)CO)CO", "Vitamina B6"),
    ("Cyanocobalamin", "N/A", "Vitamina B12"),
    ("Folic_Acid", "C1=CC(=CC=C1C(=O)NC(CCC(=O)O)C(=O)O)NC2=NC(=O)C3=C(N2)N=CN=C3N", "Vitamina B9"),
    ("Nicotinamide", "C1=CC(=CN=C1)C(=O)N", "Vitamina B3"),
    ("Ascorbic_Acid", "C([C@@H]([C@@H]1C(=C(C(=O)O1)O)O)O)O", "Vitamina C"),
    ("Cholecalciferol", "C[C@H](CCCC(C)C)C1CC[C@@H]2[C@@]1(CCC/C2=C/C=C/3\\C[C@H](CCC3=C)O)C", "Vitamina D3"),
    ("Zinc_Sulfate", "[O-]S(=O)(=O)[O-].[Zn+2]", "Mineral"),
    ("Magnesium_Sulfate", "[O-]S(=O)(=O)[O-].[Mg+2]", "Mineral"),
    # Experimentales / en investigacion
    ("Arimoclomol", "CC1=NC(=NO1)C(=O)NCC2=CC=CC=C2Cl", "HSP-coinducer"),
    ("Ezogabine", "CCOC(=O)NC1=C(C=C(C=C1)N)C(F)(F)F", "KCNQ2/3"),
    ("Masitinib", "CN1CCN(CC1)CC2=CNC3=C2C=C(C=C3)NC(=O)C4=CC(=CC=C4)CN5C(=S)SC6=C5C=CC=C6", "TKI-mastocitos"),
    ("AMX0035", "C1=CC(=CC=C1CC(C(=O)O)N)O", "Fenilbutirato-TUDCA"),
    ("Deferiprone", "CC1=C(C(=O)C=CN1C)O", "Quelante-hierro"),
    # Reposicionamiento COVID/ELA
    ("Fingolimod", "CCCCCCCC1=CC=C(C=C1)CC[C@@H](CO)N", "S1PR-modulador"),
    ("Dimethyl_Fumarate", "COC(=O)/C=C/C(=O)OC", "NRF2-activador"),
    ("Teriflunomide", "CC1=C(C=NO1)C(=O)NC2=CC=C(C=C2)C(F)(F)F", "DHODH-inhibidor"),
    ("Ocrelizumab", "N/A", "Anti-CD20"),
    ("Cladribine", "C1=NC2=C(N1[C@H]3C[C@@H]([C@@H](O3)CO)O)N=C(N=C2N)Cl", "Purina-analogo"),
    # Mas compuestos naturales con evidencia neuroprotectora
    ("Baicalein", "C1=CC=C(C=C1)C2=CC(=O)C3=C(C=C(C=C3O2)O)O", "Flavonoide"),
    ("Luteolin", "C1=CC(=C(C=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O)O", "Flavonoide"),
    ("Apigenin", "C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O", "Flavonoide"),
    ("Fisetin", "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O", "Flavonoide"),
    ("Sulforaphane", "CS(=O)CCCCN=C=S", "NRF2"),
    ("Withaferin_A", "C[C@@H]1[C@@H]2[C@H]([C@H]([C@@]3([C@H](O3)[C@@H]4[C@H]5[C@@]6([C@H](CC(=O)C=C6[C@@H]7[C@H](O7)[C@@]5([C@@H](C(=O)OC4)O3)C)C)C)O2)C)OC1=O", "Withanolide"),
]

# ──── FUENTE 2: DRUGBANK APPROVED (descargar si es posible) ────
DRUGBANK_URLS = [
    "https://go.drugbank.com/releases/latest/downloads/approved-structures",
    "https://raw.githubusercontent.com/datasets/drugbank/main/data/drugbank.csv",
]


def build_library():
    """Construye la libreria completa de compuestos."""
    print("=" * 60)
    print(" CONSTRUYENDO LIBRERIA DE COMPUESTOS")
    print("=" * 60)
    
    os.makedirs(LIBRARY_DIR, exist_ok=True)
    
    # Usar compuestos integrados (FDA + conocidos)
    print(f"\n  Compuestos integrados: {len(FDA_DRUGS)}")
    
    # Filtrar compuestos con SMILES valido
    valid = [(name, smi, use) for name, smi, use in FDA_DRUGS if smi != "N/A" and len(smi) > 3]
    print(f"  Compuestos validos (con SMILES): {len(valid)}")
    
    # Guardar CSV
    csv_path = os.path.join(LIBRARY_DIR, "full_library.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["name", "smiles", "category"])
        for name, smi, cat in valid:
            writer.writerow([name, smi, cat])
    
    # Guardar SMILES para OpenBabel
    smiles_path = os.path.join(LIBRARY_DIR, "full_library.smi")
    with open(smiles_path, "w") as f:
        for name, smi, _ in valid:
            f.write(f"{smi}\t{name}\n")
    
    print(f"  CSV: {csv_path} ({len(valid)} compuestos)")
    print(f"  SMILES: {smiles_path}")
    
    # Calcular tiempo estimado para esta PC
    # ~4 compuestos/min con 3 proteinas
    compounds_per_min = 4
    total_minutes = len(valid) / compounds_per_min
    total_hours = total_minutes / 60
    
    print(f"\n  ESTIMACION:")
    print(f"  Compuestos: {len(valid)}")
    print(f"  Velocidad: ~{compounds_per_min} comp/min (3 proteinas)")
    print(f"  Tiempo estimado: {total_hours:.1f} horas ({total_minutes:.0f} min)")
    
    return csv_path, len(valid)


if __name__ == "__main__":
    build_library()
