import csv


class BiomarkerDiscovery:
    def __init__(self, file_path):
        self.file_path = file_path
        self.records = self.load_data()

    def load_data(self):
        records = []

        with open(self.file_path, "r", encoding="utf-8") as file:
            reader = csv.DictReader(file)
            for row in reader:
                records.append(row)

        return records

    def get_gene_names(self):
        columns = list(self.records[0].keys())
        return columns[2:]

    def calculate_group_means(self):
        genes = self.get_gene_names()

        healthy_samples = [
            record for record in self.records
            if record["group"].lower() == "healthy"
        ]

        disease_samples = [
            record for record in self.records
            if record["group"].lower() == "disease"
        ]

        results = {}

        for gene in genes:
            healthy_values = [float(sample[gene]) for sample in healthy_samples]
            disease_values = [float(sample[gene]) for sample in disease_samples]

            healthy_mean = sum(healthy_values) / len(healthy_values)
            disease_mean = sum(disease_values) / len(disease_values)

            results[gene] = {
                "healthy_mean": round(healthy_mean, 2),
                "disease_mean": round(disease_mean, 2),
                "difference": round(disease_mean - healthy_mean, 2)
            }

        return results

    def rank_candidate_biomarkers(self):
        group_means = self.calculate_group_means()

        ranked_genes = sorted(
            group_means.items(),
            key=lambda item: abs(item[1]["difference"]),
            reverse=True
        )

        return ranked_genes
