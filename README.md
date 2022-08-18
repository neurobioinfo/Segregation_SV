# Segregation_SV
A tool to flexibly merge Copy Number Variants and other Structural Variants from multiple callers in multiple samples, annotate with various external data sources, and perform family segregation

```
usage: segregation.py [-h] -i INPUTDATA -p PEDFILE [-r RUNNAME]
                      [-a CUSTOMANNODATA] [-o OUTFILE] [-cfg CFGFILE] [-eda]
                      [-min MINSIZEFILTER] [-max MAXSIZEFILTER]
                      [-e EMITSTATUS] [-v]

run segregation analysis on structural variants dataset

optional arguments:
  -h, --help            show this help message and exit
  -r RUNNAME, --run_name RUNNAME
                        run name prefix for all output files
                        ["Segregation_SV"]
  -a CUSTOMANNODATA, --custom_annotation_file CUSTOMANNODATA
                        /path/to/custom_annotation_file:anno_type:minimum_reci
                        procal_overlap_fraction ; (can be specified multiple
                        times; file must be .bed or .bedpe format)
  -o OUTFILE, --out_file OUTFILE
                        /path/to/output/file; where outputs will be written
                        [/dev/stdout]
  -cfg CFGFILE, --config_file CFGFILE
                        path/to/cfg_file (used to modify run parameters)
  -eda, --exclude_default_annotations
                        exclude default annotations
  -min MINSIZEFILTER, --min_size_filter MINSIZEFILTER
                        minimal size of raw variant to keep (smaller will be
                        excluded from any processing)
  -max MAXSIZEFILTER, --max_size_filter MAXSIZEFILTER
                        maximal size of raw variant to keep (larger will be
                        excluded from any processing
  -e EMITSTATUS, --emit_vars_from_samples_with_status EMITSTATUS
                        emit variants from samples with clinical status of
                        [affected,nonaffected,all]
  -v, --verbose         verbose output

mandatory arguments:
  -i INPUTDATA, --input_data INPUTDATA
                        /path/to/input_data_file:input_type:[batch_name] (can
                        be specified multiple times; accepted values are
                        "lumpy", "popsv", "conifer", "codex", "exomedepth",
                        "cnmops", "xhmm", "generic_bed" and "generic_bedpe")
  -p PEDFILE, --ped_file PEDFILE
                        /path/to/ped_file
```

# Principle
Segregation_SV works by first converting individual per-sample CNV calls from multiple callers into “consensus CNVs” based on the principle of reciprocal overlap (RO). This merging is done in 3 successive steps:

- Step 1: calls by the same caller in the same sample are merged if there is any overlap between them (>0% RO), to avoid call fragmentation. 
- Step 2: calls by different callers in the same sample are merged if there is >=50% RO. 
- Step 3: Once this has been done for all samples, calls are then merged between multiple samples if there is >80% RO between them. 

Finally for each consensus CNV, sample counts are broken down into affected vs nonaffected and family vs non-family, for each family as defined by an input pedigree file

# Notes:
- At each step, CNV merging proceeds iteratively from smallest to largest, to ensure that small CNVs are not accidentally excluded.
- Pedigree file format is described [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format)
- The config file (specified with command argument -cfg|--config_file) can be used to very precisely tailor the variant merging rules. Different RO values can be set for each variant merging step (steps 1-3 above), as well as for different variant types (deletion, duplication, inversion and translocation) within each merging step.
- Variant merging for translocations is not supported. Any such variants are simply passed along intact through all merging steps into the final outputs, though if >1 sample has a called translocation with _exactly_ the same boundaries, this will be reflected in the final sample counts for that variant. 
- Annotations included by default are: 
