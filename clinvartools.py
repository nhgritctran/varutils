import requests
import pandas as pd
from tqdm import tqdm


class DataProcessing:

    def __init__(self):
        self.spdi_base_url = "https://api.ncbi.nlm.nih.gov/variation/v0/spdi/"

    def spdi2vcf(self, spdi_string, output_type):
        """
        take spdi format and convert to vcf format
        :param spdi_string: variant data in SPDI format
        :param output_type: accepts "locus", "allele_change", "var"
        :return: vcf format data
        """

        # make spdi api call
        base_url = self.spdi_base_url
        complete_url = base_url + spdi_string + "/vcf_fields"
        resp = requests.get(url=complete_url)
        data = resp.json()

        # process output
        chr_num = data["data"]["chrom"].split(".")[0].split("_")[1].lstrip("0")
        locus = "chr" + \
                chr_num + ":" + \
                str(data["data"]["pos"])
        allele_change = data["data"]["ref"] + ":" + \
                        data["data"]["alt"]
        var = "chr" + \
              chr_num + ":" + \
              str(data["data"]["pos"]) + ":" + \
              data["data"]["ref"] + ":" + \
              data["data"]["alt"]

        if output_type == "locus":
            return locus
        elif output_type == "allele_change":
            return allele_change
        elif output_type == "var":
            return var
        else:
            print("Invalid output type.")

    def process_clinvar_output(self, clinvar_file_path):
        """
        take ClinVar output .txt file and process data ClinVar output file must not be altered for this to work.
        :param clinvar_file_path: path to ClinVar .txt output
        :return: new dataframe with new columns;
                 column "hail_variant" contains variant data match that of Hail matrix table.
        """

        tqdm.pandas()

        # load data
        df = pd.read_csv(clinvar_file_path, sep="\t")

        # process data
        df["refseq"] = df["Name"].str.split("(").str[0]
        df["locus_change"] = df["Name"].str.split("(").str[1].str.split("c.").str[1]
        df["clinical_significance"] = df["Clinical significance (Last reviewed)"].str.split("(").str[0]
        df["ClinVar"] = df["clinical_significance"].copy()
        df["last_review_date"] = df["Clinical significance (Last reviewed)"] \
                                 .str.split(":").str[1] \
                                 .str.split(")").str[0]

        # group clinical significance phenotypes
        p_list = ["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"]
        b_list = ["Benign", "Likely benign", "Benign/Likely benign"]
        df["Phenotype"] = "Uncertain/Conflicting"
        df["Phenotype"].loc[df["clinical_significance"].isin(b_list)] = "Benign (B/LB)"
        df["Phenotype"].loc[df["clinical_significance"].isin(p_list)] = "Pathogenic (P/LP)"

        # remove NAs in Canonical SPDI column
        df = df.loc[~(df["Canonical SPDI"].isna())]

        # add column "snv" to specify if a mutation is snv or not
        df["snv"] = 0
        df["snv"].loc[df["locus_change"].str.contains(">")] = 1

        # create column hail_variant
        df["hail_variant"] = df["Canonical SPDI"].copy()
        df["hail_variant"] = df["hail_variant"].progress_apply(lambda x: self.spdi2vcf(x, "var"))

        df["grch38_locus"] = df["hail_variant"].str.split(":").str[0] + ":" + \
                             df["hail_variant"].str.split(":").str[1]

        df["alleles"] = df["hail_variant"].str.split(":").str[2] + ":" + \
                        df["hail_variant"].str.split(":").str[3]

        # select and change column names
        cols = ["Name",
                "Gene(s)",
                "dbSNP ID",
                "Condition(s)",
                "clinical_significance",
                "last_review_date",
                "Review status",
                "locus_change",
                "Protein change",
                "refseq",
                "Accession",
                "GRCh37Chromosome",
                "GRCh37Location",
                "GRCh38Chromosome",
                "GRCh38Location",
                "grch38_locus",
                "alleles",
                "snv",
                "hail_variant",
                "VariationID",
                "AlleleID(s)",
                "ClinVar",
                "Phenotype",
                "Canonical SPDI"]
        new_cols = ["name",
                    "gene",
                    "dbsnp_rsid",
                    "condition",
                    "clinical_significance",
                    "last_review_date",
                    "review_status",
                    "locus_change",
                    "protein_change",
                    "refseq",
                    "accession",
                    "grch37_chromosome",
                    "grch37_position",
                    "grch38_chromosome",
                    "grch38_position",
                    "grch38_locus",
                    "alleles",
                    "snv",
                    "hail_variant",
                    "variation_id",
                    "allele_id",
                    "ClinVar",
                    "Phenotype",
                    "spdi"]
        df = df[cols]
        df.columns = new_cols

        return df
