import pandas as pd
import random
import math
random.seed(42)

# Read the file
biogrid = pd.read_csv("/Users/alexraggi/Desktop/BIOGRID-ALL-5.0.251.tab.txt", sep = "\t", low_memory = False, skiprows = 35)
# Remove non-human proteins

biogrid = biogrid[(biogrid["ORGANISM_A_ID"] == 9606) & (biogrid["ORGANISM_B_ID"] == 9606)]
# Remove self-interactions
biogrid = biogrid[biogrid["OFFICIAL_SYMBOL_A"] != biogrid["OFFICIAL_SYMBOL_B"]]

# Create a set of all the good proteins
with open("/Users/alexraggi/Desktop/ACADEMIC/TERMS/FALL 2025/BIEN 203/Assignments/Assignment 5/good-proteins.txt") as f:
    good_proteins = {line.split("\t")[0] for line in f}

# Only keep pairs of proteins where both proteins are considered good proteins
biogrid = biogrid[(biogrid["OFFICIAL_SYMBOL_A"].isin(good_proteins)) & (biogrid["OFFICIAL_SYMBOL_B"].isin(good_proteins))]

# Create the positive pairs
positive_pairs = {tuple(sorted((a, b))) for a, b in zip(biogrid["OFFICIAL_SYMBOL_A"], biogrid["OFFICIAL_SYMBOL_B"])}

# Create the negative pairs
list_positive_pairs = list(positive_pairs)
negative_pairs = set()

while len(negative_pairs) < len(positive_pairs):
    a, b = random.sample(list(good_proteins), 2)
    pair = tuple(sorted((a, b)))

    if pair not in positive_pairs:
        negative_pairs.add(pair)

list_negative_pairs = list(negative_pairs)

# Creating both positive and negative pairs DataFrame
positive_df = pd.DataFrame(list_positive_pairs, columns = ["Protein_A", "Protein_B"])
negative_df = pd.DataFrame(list_negative_pairs, columns = ["Protein_A", "Protein_B"])

# Adding the interaction Column
positive_df["Interaction"] = 1
negative_df["Interaction"] = 0

# Adding protocorr and depmapcorr values
protcorr = pd.read_csv("protcorr.tsv", sep = "\t", header = None)
protcorr.columns = ["Protein_A", "Protein_B", "Protcorr"]

depmapcorr = pd.read_csv("depmapcorr.tsv", sep = "\t", header = None)
depmapcorr.columns = ["Protein_A", "Protein_B", "Depmapcorr"]

# Dropping na
protcorr = protcorr.dropna(subset=["Protein_A", "Protein_B"])
depmapcorr = depmapcorr.dropna(subset=["Protein_A", "Protein_B"])

protcorr[["Protein_A", "Protein_B"]] = protcorr.apply(lambda row: sorted([row["Protein_A"], row["Protein_B"]]), axis=1, result_type="expand")
depmapcorr[["Protein_A", "Protein_B"]] = depmapcorr.apply(lambda row: sorted([row["Protein_A"], row["Protein_B"]]), axis=1, result_type="expand")

all_pairs = pd.concat([positive_df, negative_df], ignore_index = True)

# Merging both protocorr and depmapcorr into the initial dataframe
all_pairs = all_pairs.merge(protcorr, on = ["Protein_A", "Protein_B"], how = "left")
all_pairs = all_pairs.merge(depmapcorr, on = ["Protein_A", "Protein_B"], how = "left")

# Filling the missing values with 0s
all_pairs = all_pairs.fillna(0)

# Creating a dictionary where each key is a protein and its value is a set of all its interactions
neighbors = {}

for a, b in positive_pairs:
    neighbors.setdefault(a, set()).add(b)
    neighbors.setdefault(b, set()).add(a)

lower_degrees = []
higher_degrees = []

for a, b in zip(all_pairs["Protein_A"], all_pairs["Protein_B"]):
    degree_a = len(neighbors.get(a, set()) - {b})
    degree_b = len(neighbors.get(b, set()) - {a})

    lower_degrees.append(min(degree_a, degree_b))
    higher_degrees.append(max(degree_a, degree_b))

all_pairs["Lower_Degree"] = lower_degrees
all_pairs["Higher_Degree"] = higher_degrees

# Finding the Jaccard Index for every pair of protein

jaccard_values = []
for a, b in zip(all_pairs["Protein_A"], all_pairs["Protein_B"]):
    neighbors_a = neighbors.get(a, set()) - {b}
    neighbors_b = neighbors.get(b ,set()) - {a}
    intersection = neighbors_a & neighbors_b
    union = neighbors_a | neighbors_b
    if len(union) == 0:
        jaccard = 0
    else:
        jaccard = len(intersection) / len(union)
    jaccard_values.append(jaccard)
all_pairs["Jaccard"] = jaccard_values



gene2pubmed = pd.read_csv("gene2pubmed.tsv", sep = "\t", header = None)
gene2pubmed.columns = (["Gene_ID", "Article_ID"])

# Creating a Series showing the count of proteins per papaer
paper_counts = gene2pubmed["Article_ID"].value_counts()

large_papers = paper_counts[paper_counts >= 20].index
gene2pubmed = gene2pubmed[~gene2pubmed["Article_ID"].isin(large_papers)]

# gene -> set of PubMed article IDs
papers_by_gene = {}

for gene, article in zip(gene2pubmed["Gene_ID"], gene2pubmed["Article_ID"]):
    papers_by_gene.setdefault(gene, set()).add(article)

shared_papers = []

for a, b in zip(all_pairs["Protein_A"], all_pairs["Protein_B"]):
    papers_a = papers_by_gene.get(a, set())
    papers_b = papers_by_gene.get(b, set())
    common_papers = papers_a & papers_b
    shared_papers.append(len(common_papers))

all_pairs["Common_Papers"] = shared_papers

# Loadin the domain file
domains_by_protein = {}

with open("pfam-domains.tsv") as f:
    for line in f:
        parts = line.strip().split("\t")
        protein = parts[0]
        domains = parts[2:]
        
        domains_by_protein[protein] = set(domains)

positive_domain_counts = {}
negative_domain_counts = {}

# Count domain pairs in positive protein pairs
for a, b in positive_pairs:
    domains_a = domains_by_protein.get(a, set())
    domains_b = domains_by_protein.get(b, set())

    for d1 in domains_a:
        for d2 in domains_b:
            domain_pair = tuple(sorted((d1, d2)))
            positive_domain_counts[domain_pair] = positive_domain_counts.get(domain_pair, 0) + 1

# Count domain pairs in negative pairs
for a, b in negative_pairs:
    domains_a = domains_by_protein.get(a, set())
    domains_b = domains_by_protein.get(b, set())
    for d1 in domains_a:
        for d2 in domains_b:
            domain_pair = tuple(sorted((d1, d2)))
            negative_domain_counts[domain_pair] = negative_domain_counts.get(domain_pair, 0) + 1

all_domain_pairs = set(positive_domain_counts) | set(negative_domain_counts)


domain_llr = {}

for domain_pair in all_domain_pairs:
    positive_count = positive_domain_counts.get(domain_pair, 0)
    negative_count = negative_domain_counts.get(domain_pair, 0)
    llr = math.log((positive_count + 0.5) / (negative_count + 0.5))
    domain_llr[domain_pair] = llr


max_domain_llrs = []

for protein_a, protein_b in zip(all_pairs["Protein_A"], all_pairs["Protein_B"]):
    domains_a = domains_by_protein.get(protein_a, set())
    domains_b = domains_by_protein.get(protein_b, set())

    llrs = []

    for domain_a in domains_a:
        for domain_b in domains_b:
            domain_pair = tuple(sorted((domain_a, domain_b)))

            if domain_pair in domain_llr:
                llrs.append(domain_llr[domain_pair])

    if len(llrs) == 0:
        max_llr = 0
    else:
        max_llr = max(llrs)

    max_domain_llrs.append(max_llr)

all_pairs["Max_Domain_LLR"] = max_domain_llrs


from sklearn.model_selection import train_test_split
feature_columns = ["Protcorr", "Depmapcorr", "Lower_Degree", "Higher_Degree", "Jaccard", "Common_Papers", "Max_Domain_LLR"]
X = all_pairs[feature_columns]
y = all_pairs["Interaction"]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state = 42, stratify = y)


from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

logistic_model = make_pipeline(StandardScaler(), LogisticRegression())
logistic_model.fit(X_train, y_train)
logistic_prob = logistic_model.predict_proba(X_test)[:, 1]

from sklearn.metrics import roc_auc_score

logistic_auc = roc_auc_score(y_test, logistic_prob)

from sklearn.metrics import roc_curve

fpr_log, tpr_log, thresholds = roc_curve(y_test, logistic_prob)

from sklearn.ensemble import RandomForestClassifier

forest_model = RandomForestClassifier(n_estimators=100, random_state=42)
forest_model.fit(X_train, y_train)
forest_prob = forest_model.predict_proba(X_test)[:, 1]
forest_auc = roc_auc_score(y_test, forest_prob)

print("Logistic Regression AUC Score:", logistic_auc)
print("Forest AUC Score:", forest_auc)

fpr_forest, tpr_forest, thresholds_forest = roc_curve(y_test, forest_prob)

import matplotlib.pyplot as plt


plt.plot(fpr_log, tpr_log, label=f"Logistic Regression (AUC = {logistic_auc:.3f})")

plt.plot(fpr_forest, tpr_forest, label=f"Random Forest (AUC = {forest_auc:.3f})")

plt.plot([0, 1], [0, 1], "--", label="Random classifier")

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve")
plt.legend()

plt.show()

