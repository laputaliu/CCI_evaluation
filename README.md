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

<br/>

## Evaluation pipeline
Users can evaluate CCI tools using their own data following the workflow below.

### Step 1: compute the d_rat and P-value for LR pairs using ST data
Users need to provie a count file for ST data, a "coordinates" file recording the coordinates of each spot, and a "LR pair" file recording the ligand-receptor interactions. 
- The input "coordinates" file should be a tab-delimited text file with spot barcodes in the first column, the x coordinates of spot in the second column (named as row), and the y coordinated of the spot in the third column (named as col). 
- The "LR pair" file is a LR database in the text format, it should be a tab-delimited text file recording the ligand, receptor and their interactions. You can simply used the "LR pair" file generated from the CellChatDB which is used in our study, you can find it at `./cci_database/cc_ip_multi_split.tsv`. Or you can also use your customized LR database. There are 4 columns needed in the "LR pair", the first column should be the interaction name, you can connect partipated L,R genes with the `_`; next a column named "interaction_name_2" for interaction names in another format (`ligand - receptor` for single subunit, `ligand - (receptor1+receptor2+...)` for multi-subunits); and another two columns recording ligand and receptor gene symbols, for multi-subunit receptors, the subunits need to be seprated by `,`.

<br/>

Example input "coordinates" file:
```
	row	col
CTATCGGGTCTCAACA-1	35	69
TCGAGACCAACACCGT-1	37	55
GATGCGTCCTGCATTC-1	37	57
GCATAGAGCACTCAGG-1	37	59
GCAGATTAGGGATATC-1	37	61
```

<br/>

Example input "LR pair" file (CellChatDB):
```
	pathway_name	annotation	interaction_name_2	ligand	receptor
TGFB1_TGFBR1_TGFBR2	TGFb	Secreted Signaling	TGFB1 - (TGFBR1+TGFBR2)	TGFB1	TGFBR1,TGFBR2
TGFB2_TGFBR1_TGFBR2	TGFb	Secreted Signaling	TGFB2 - (TGFBR1+TGFBR2)	TGFB2	TGFBR1,TGFBR2
TGFB3_TGFBR1_TGFBR2	TGFb	Secreted Signaling	TGFB3 - (TGFBR1+TGFBR2)	TGFB3	TGFBR1,TGFBR2
TGFB1_ACVR1B_TGFBR2	TGFb	Secreted Signaling	TGFB1 - (ACVR1B+TGFBR2)	TGFB1	ACVR1B,TGFBR2
GDF15_TGFBR2	GDF	Secreted Signaling	GDF15 - TGFBR2	GDF15	TGFBR2
GDNF_GFRA1	GDNF	Secreted Signaling	GDNF - GFRA1	GDNF	GFRA1
```

<br/>

when having all these input files prepared, you can submit computing d_rat and P-value of each interaction by using the script we provided, the example command below:
```
nohup python ./scripts/prepare_ip_dis_sinkhorn2.py -c st_count_file -p st_coord_file -o ./ip_dis_sinkhorn2 &
```
Since this step will cost a lot of time, for a quick start, you can use the ip_dis file that we generated before. you can find it in the `./example_data/ST_A3_GSM4797918/data/ip_dis_sinkhorn2/ip_distance_all.tsv`, which computed the d_rat and P-value using the example ST data from sample ST_A3_GSM4797918.







