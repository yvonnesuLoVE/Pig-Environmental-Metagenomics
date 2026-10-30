# Air and environmental sampling expand metagenomic viral surveillance in Cambodian pig slaughterhouses

## Overview
This repository contains the analysis code used for the study:

**"Air and environmental sampling expand metagenomic viral surveillance in Cambodian pig slaughterhouses"**

This study investigates how air and environmental sampling can expand viral surveillance beyond direct animal sampling in Cambodian pig slaughterhouses. Direct animal sampling describes selected hosts but can miss the wider environment through which viral material is shed, mixed and transported.

Across nine visits to two Cambodian pig slaughterhouses, we applied capture-enriched metagenomic sequencing to air, faeces, wastewater, carcass wash water, holding-pen swabs, slaughter-box swabs and ropes. Viral detections from these environmental samples were compared with those from concurrently collected pig nasal swabs to assess surveillance breadth and the contribution of complementary sample types.


## Key Findings

- Samples comprised air (n=26), faeces (n=146), wastewater (n=15), carcass wash water (n=18), holding-pen swabs (n=134), slaughter-box swabs (n=63), ropes (n=16) and pig nasal swabs (n=200).
- Combined environmental sampling recovered more than 96% of the viral taxa detected in pig nasal swabs at every visit, while also detecting additional swine-associated taxa and viral RNA associated with avian and human hosts.
- Air provided the most informative signal from a single environmental platform, recovering more than 77% of the taxa detected in pig nasal swabs, retaining many common swine-associated viruses and contributing the largest number of detections not recovered from any other environmental sample type.
- Other environmental samples captured different shedding routes and stages of slaughter, demonstrating the value of complementary sampling approaches.
- Slaughter-box swabs recovered the fewest taxa, showing that visible exposure to animal material did not guarantee surveillance breadth.

---

## Data Sources
- Air and environmental samples collected across nine visits to two Cambodian pig slaughterhouses.
- Faeces, wastewater, carcass wash water, holding-pen swabs, slaughter-box swabs and ropes representing different shedding routes and stages of slaughter.
- Concurrently collected pig nasal swabs used as the direct animal-sampling comparator.
- Capture-enriched metagenomic sequencing used to compare viral detections across sample types and visits.

---

## Repository Contents
- `Pig_Script.R`: analysis script covering CZ ID count processing, taxon filtering, reads-per-million calculations, presence/absence matrices, detection overlap and unique detections, alpha and beta diversity, paired pig–environment comparisons, slaughterhouse comparisons, DIAMOND annotation summaries and phylogenetic tree visualisation.
- `README.md`: study overview, key findings, data sources and code-use information.
- Required inputs are not included. The script uses count, detection, metadata and taxonomy workbooks (`Pig_*.xlsx`), raw CZ ID CSV exports in `Working_Data/Pig_Counts/` and `Working_Data/Pig_Contig_Taxonomy/`, annotated workbooks in `DIAMOND/`, and precomputed tree files in `Trees/`. The script header describes the expected layout.
- To use the code, start R in the project/data folder, install the packages listed in `required_packages`, run the setup section and then run the relevant numbered sections. `phyloseq`, `ComplexHeatmap` and `ggtree` are Bioconductor packages. File locations are configured once in the setup; active figure and table exports are written to `results/`.
- Preprocessing checkpoint exports remain commented as in the original analysis. To run these sections sequentially, enable the corresponding checkpoint exports before subsequent sections read their workbooks, or supply the existing intermediate workbooks. DIAMOND annotations and phylogenetic trees must be supplied separately.

---

## Implications
Air sampling expands the reach of non-contact viral surveillance in pig slaughterhouses. Interpreting several complementary environmental sample types together provides a broader view of viral detections than any single platform, supporting surveillance approaches that account for different shedding routes and stages of slaughter.

## Citation
If you use this repository in your research, please cite:

**[Your Full Citation Here]**
