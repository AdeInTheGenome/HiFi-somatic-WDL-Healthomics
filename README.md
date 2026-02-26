# HiFi Somatic WDL – SRS Adaptation

**SRS Research Adaptation — Not an Official PacBio Release**

This repository is a fork and research adaptation of the original  
**Pacific Biosciences HiFi Somatic WDL HealthOmics workflow**:

https://github.com/PacificBiosciences/HiFi-somatic-WDL

It was modified for the Somatic Reference Samples (SRS) Initiative to support long-read somatic variant analysis for engineered HG002 clones and validation datasets described in:

Daniels et al., *Somatic Reference Sample Development and Evaluation Using Unedited and CRISPR-Cas9 Edited Human Cell Lines* (in preparation)

---

## Purpose of This Fork

The original PacBio HiFi somatic workflow was adapted to meet the benchmarking and orthogonal validation needs of the SRS program.

Primary updates include:

- Integration of **Sentieon TNscope** for long-read somatic SNV/INDEL calling
- Modifications to support SRS tumor–normal experimental design
- Parameterization adjustments for AWS HealthOmics execution within the SRS cloud infrastructure
- Harmonization with SRS short-read pipelines for cross-technology comparison

This repository reflects the version used for SRS research analyses and manuscript generation. It is not intended as a production pipeline or clinical workflow.

---

## Acknowledgment

We thank the **HiFi Somatic workflow development team at Pacific Biosciences** for making the original workflow publicly available. Their work provided the foundation for long-read somatic benchmarking within the SRS initiative.

All credit for the original architecture and design belongs to the PacBio team. This repository contains only SRS-specific modifications layered on top of their original implementation.

For the official PacBio workflow, please refer to:
https://github.com/PacificBiosciences/HiFi-somatic-WDL
