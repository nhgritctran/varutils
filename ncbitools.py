from concurrent.futures import as_completed
from concurrent.futures import ThreadPoolExecutor
from IPython.display import display, clear_output
from ratelimit import limits, sleep_and_retry
from tqdm import tqdm

import bz2
import json
import os
import requests


class VariationServices:
    def __init__(self):
        pass

    @staticmethod
    @sleep_and_retry
    @limits(calls=5, period=1)
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

    @sleep_and_retry
    @limits(calls=5, period=1)
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

        vcf_dict = {"rsid": rsid,
                    "vcf": [self.spdi_to_vcf(s, "var") for s in spdi]}

        return vcf_dict

    def multi_rsid_to_vcf(self, rsid_list):
        """
        multi-threaded version of rsid_to_vcf
        :param rsid_list: a list of rsID
        :return: list of dict {rsid: [variant(s) in vcf format]}
        """
        jobs = []
        with ThreadPoolExecutor(max_workers=os.cpu_count()-1) as executor:
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

    @staticmethod
    def _parse_rs_obj(rs_obj):
        """
        parse data from rs_obj
        :param rs_obj: a line from dbSNP bz2 file
        :return: dictionary with variant information
        """
        # process rs_obj
        main_obj = rs_obj["primary_snapshot_data"]["placements_with_allele"][0]
        assembly_name = main_obj["placement_annot"]["seq_id_traits_by_assembly"][0]["assembly_name"]
        alleles = main_obj["primary_snapshot_data"]["placements_with_allele"][0]["alleles"]
        spdi = []
        for allele in alleles:  # get spdi of variants
            if "=" not in allele["hgvs"]:
                spdi.append(allele["allele"]["spdi"]["seq_id"] \
                            + ":" \
                            + str(allele["allele"]["spdi"]["position"]) \
                            + ":" \
                            + allele["allele"]["spdi"]["deleted_sequence"] \
                            + ":" \
                            + allele["allele"]["spdi"]["inserted_sequence"])

        # put information into dict
        rs_dict = {"rsid": rs_obj["refsnp_id"],
                   "assembly_name": assembly_name,
                   "spdi": spdi}

        return rs_dict

    def parse_dbsnp_bz2_file(self, bz2_file_path):
        jobs = []
        with ThreadPoolExecutor(max_workers=os.cpu_count()-1) as executor:
            print("Submitting jobs:")
            with bz2.BZ2File(bz2_file_path, "rb") as f:
                for line in tqdm(f):
                    rs_obj = json.loads(line.decode('utf-8'))
                    jobs.append(executor.submit(self._parse_rs_obj, rs_obj))

            n = 1
            for job in as_completed(jobs):
                clear_output(wait=True)
                display(f"Completed {n} jobs!")
                n += 1

            result_dicts = [job.result() for job in jobs]

        return result_dicts

    def bz2_rsid_to_spdi(self, rsid_file_path, bz2_file_path):
        """
        used to convert multiple rsID to spdi format using bz2 dbSNP chromosome data
        :param rsid_file_path: path to rsID file, should be in txt format, one rsid per line
        :param bz2_file_path: path to dbSNP bz2 file
        :return: spdi format of all variants of each rsid
        """
        # parse rsid file to list
        print("Parse input file.")
        with open(rsid_file_path, "r") as rsid_file:
            rsid_list = []
            for line in rsid_file:
                rsid_list.append(line.rstrip())
        refsnp_list = [int(rsid.replace("rs", "")) for rsid in rsid_list]
        refsnp_list.sort()
        print("Done!")
        print()

        # read file
        with bz2.BZ2File(bz2_file_path, "rb") as f:
            line_count = 0
            spdi_dict = {}

            print("Start mapping refSNP to SPDI")
            # loop refsnp list
            for refsnp in tqdm(refsnp_list):

                # loop dbSNP file and find match
                for i, line in tqdm(enumerate(f)):
                    if i >= line_count:  # only look at unchecked lines
                        rs_obj = json.loads(line.decode('utf-8'))

                        if refsnp == rs_obj["refsnp_id"]:  # if match found
                            rs_dict = self._parse_rs_obj(rs_obj)
                            line_count += 1
                            break

                        else:
                            line_count += 1

                spdi_dict[f"rs{refsnp}"] = rs_dict["spdi"]

        return spdi_dict
