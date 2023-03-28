from concurrent.futures import as_completed
from concurrent.futures import ThreadPoolExecutor
from IPython.display import display, clear_output
from ratelimit import limits, sleep_and_retry
from tqdm import tqdm

import os
import requests


class VariationServices:
    def __init__(self):
        pass

    @staticmethod
    def spdi_to_vcf(spdi, output_type):
        """
        convert one spdi format to information from vcf format
        :param spdi: spdi format
        :param output_type: accepts "locus", "allele_change" and "var"
        :return: depending on output_type
        """
        base_url = "https://api.ncbi.nlm.nih.gov/variation/v0/spdi/"
        complete_url = base_url + spdi + "/vcf_fields"
        resp = requests.get(url=complete_url)
        data = resp.json()
        chr_num = data["data"]["chrom"].split(".")[0].split("_")[1].lstrip("0")

        locus = "chr" + \
                chr_num + ":" + \
                str(data["data"]["pos"])

        allele_change = data["data"]["ref"] + ":" + data["data"]["alt"]

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

    def rsid_to_vcf(self, rsid):
        """
        convert one rsid to vcf format (chr:pos:ref:alt)
        :param rsid: rsid of variant
        :return: all variants associated with the rsid, in vcf format
        """
        base_url = "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/"
        rsid_num = rsid.replace("rs", "")
        complete_url = base_url + str(rsid_num)
        resp = requests.get(url=complete_url)
        data = resp.json()

        assembly = data["primary_snapshot_data"]["placements_with_allele"][0] \
            ["placement_annot"]["seq_id_traits_by_assembly"][0]["assembly_name"]

        spdi = []
        if "GRCh38" in assembly:
            alleles = data["primary_snapshot_data"]["placements_with_allele"][0]["alleles"]
            for allele in alleles:
                if "=" not in allele["hgvs"]:
                    spdi.append(allele["allele"]["spdi"]["seq_id"] \
                                + ":" \
                                + str(allele["allele"]["spdi"]["position"]) \
                                + ":" \
                                + allele["allele"]["spdi"]["deleted_sequence"] \
                                + ":" \
                                + allele["allele"]["spdi"]["inserted_sequence"])
        else:
            spdi.append("")

        vcf_dict = {rsid: [self.spdi_to_vcf(s, "var") for s in spdi]}

        return vcf_dict

    def multi_rsid_to_vcf(self, rsid_list):
        """
        this function is multi-threaded version of rsid_to_vcf
        :param rsid_list: a list of rsid
        :return: list of dict {rsid: [variant(s) in vcf format]}
        """
        jobs = []
        with ThreadPoolExecutor(max_workers=os.cpu_count() - 1) as executor:
            print("Submitting jobs:")
            for rsid in tqdm(rsid_list):
                jobs.append(executor.submit(self.rsid_to_vcf,
                                            rsid))
            n = 1
            for job in as_completed(jobs):
                clear_output(wait=True)
                display(f"Completed {n} jobs!")
                n += 1

            result_dicts = [job.result() for job in jobs]

        return result_dicts
