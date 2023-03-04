# varutils - tools for working with genomic variants.

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

# Dependencies

There is no extra dependencies needed besides those that already available on All of Us Researcher Workbench.