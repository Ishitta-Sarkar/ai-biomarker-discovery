from src.biomarker_tools import BiomarkerDiscovery


discovery = BiomarkerDiscovery("data/gene_expression.csv")

print("AI Biomarker Discovery")
print("=" * 30)

print("\nGroup Mean Expression")
print("-" * 30)

group_means = discovery.calculate_group_means()

for gene, values in group_means.items():
    print(f"{gene}:")
    print("  Healthy Mean:", values["healthy_mean"])
    print("  Disease Mean:", values["disease_mean"])
    print("  Difference:", values["difference"])

print("\nRanked Candidate Biomarkers")
print("-" * 30)

ranked_genes = discovery.rank_candidate_biomarkers()

for rank, (gene, values) in enumerate(ranked_genes, start=1):
    print(f"{rank}. {gene} | Difference: {values['difference']}")
