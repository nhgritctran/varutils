# varutils - tools for working with genomic variants

# Description
This collection of tools were developed for All of Us jupyter notebook environment.

# Setup
Clone repository in Jupyter Notebook instance using:
```angular2html
!git clone https://github.com/nhgritctran/varutils
```
Simply import the desired tool afterwards, for example:

```angular2html
from varutils import clinvartools, hailtools
```

# Current tools

## clinvartools

### class DataProcessing
Process ClinVar txt output file downloaded from ClinVar and generate variant format that would work with Hail Matrix Table.

## hailtools

### class VariantMapping
Scan Hail matrix table for variants of interest.

## ncbitools

### class VariationServices
Utilize Variation Services API from NCBI. 

Current capabilities include mapping variant format spdi to VCF and rsid to VCF. VCF format is chr:position:ref:alt

Method multi_rsid_to_vcf can be used to map multiple rsid (recommended less than 10,000 rsid) to VCF.

# Dependencies

There is no extra dependencies needed besides those that already available on All of Us Researcher Workbench.