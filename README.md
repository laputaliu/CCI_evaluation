# Evaluation of cell-cell interaction methods by integrating single-cell RNA sequencing data with spatial information

Recent advances in single-cell RNA-sequencing (scRNA-seq) enable the characterization of cell-cell interactions (CCIs). However, it’s hard to evaluate these methods since no ground truth is provided. Here, we designed a comprehensive workflow to evaluate the accuracy of inferred CCIs from scRNA-seq data based on the expected and observed spatial interaction tendencies. In our workflow, we calculated the metrics of distance enrichment score (DES) and commonly identified interactions to give each tool a comprehensive evaluation.

## CCI evaluation workflow

![image](https://github.com/laputaliu/CCI_evaluation/blob/main/fig/workflow.png)

First, generate known ligand-receptor pairs from CellChatDB, then select short-range and long-range interactions from known pairs for each dataset (top left). Next, perform spatial cell type annotation on ST data coupled with matched scRNA-seq data, and define near and far distributed cell type pairs based on the annotation (medium left). Then feed annotated scRNA-seq data to CCI tools and extract predicted results (bottom left). Finally, evaluate tools’ performances on both distance enrichment score and the metric of commonly identified interactions.

## CCI tools included in this study
- [CellChat](https://github.com/sqjin/CellChat) v1.0.0
- [CellPhoneDB](https://github.com/Teichlab/cellphonedb) v2.0
- [Connectome](https://github.com/msraredon/Connectome) v1.0.1
- [CytoTalk](https://github.com/huBioinfo/CytoTalk) v3.1.0
- [ICELLNET](https://github.com/soumelis-lab/ICELLNET) v0.99.3
- [iTALK](https://github.com/Coolgenome/iTALK) v0.1.0
- [NATMI](https://github.com/forrest-lab/NATMI/)
- [NicheNet](https://github.com/saeyslab/nichenetr) v1.0.0
- [scMLnet](https://github.com/SunXQlab/scMLnet) v0.1.0
- [SingleCellSignalR](https://github.com/SCA-IRCM/SingleCellSignalR_v1) v1.4.0

## Datasets included in this study
- <b>Human pancreatic ductal adenocarcinoma (PDAC) dataset</b>
  - The scRNA-seq and ST data of the PDAC dataset can be accessed through GEO under accession number [GSE111672](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111672).
- <b>Human squamous cell carcinoma (SCC) dataset</b>
  - The scRNA-seq and ST data of the SCC dataset can be accessed through GEO under accession number [GSE144240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144240).
- <b>Mouse cortex dataset</b>
  - The scRNA-seq data of the mouse cortex dataset can be downloaded from [Allen Brain Atlas](https://dx.doi.org/10.1038/nn.4216). The ST data of cortex dataset can be downloaded from [10x Genomics website](https://www.10xgenomics.com/resources/datasets/adult-mouse-brain-ffpe-1-standard-1-3-0).
- <b>Human heart dataset</b>
  - The scRNA-seq and ST data of the human heart dataset can be collected from [Spatial Research website](https://www.spatialresearch.org/resources-published-datasets/doi-10-1016-j-cell-2019-11-025/).
- <b>Human intestinal dataset</b>
  - The scRNA-seq and ST data of the human intestinal dataset can be accessed through GEO under accession numbers [GSE158328](https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/query/acc.cgi?acc=GSE158328) (ST) and [GSE158702](https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/query/acc.cgi?acc=GSE158702) (scRNA-seq)

## Repository

The running scripts of CCI tools and analysis codes used in this study are available in this GitHub repository.
- `prepare_dir_script.py` this Python script is used to prepare the working directories and generate script submitting commands for every CCI tool.
- `/scripts` folder contains all the scripts for running each CCI tool.
- `extract_tool_results.py` this Python script is used to extract all the predicted interactions from each tool's outputs.



