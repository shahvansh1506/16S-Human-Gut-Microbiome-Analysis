# Archive and settings audit

The run inventory was retrieved from ENA and NCBI. The workflow uses the physical FASTQ files listed by ENA to choose a read layout. Some publication or SRA layout labels describe an earlier stage of sequencing. Counts below are archive runs and may differ from the number of participants in the paper.

| Owner | Project | What the archive shows | Decision and remaining check |
|---|---|---|---|
| Vansh | PRJNA687346 | 24 single-end MiSeq runs, about 0.14 GB. One BioSample ID is reused across specimen libraries. | Used unique library names as specimen keys. Pilot trimming and 220-base truncation reproduced the three saved samples. Checked the first 3,000 reads per file across the full set. |
| Vansh | PRJNA761547 | 322 paired runs. 187 names end in `_16S` and 135 end in `_ITS`. Selected FASTQs total about 1.05 GB. | Used all 187 bacterial 16S runs. Verified paired merging in subsets and used DADA2's `noqualErrfun` because archived quality scores had very little variation. The archive lists NovaSeq and the paper describes MiSeq PE300. |
| Bhavya | PRJNA604466 | 53 paired amplicon runs, about 0.94 GB. The project includes 14 shotgun WGS runs. | Used the 53 amplicon runs as already trimmed. Three common primer matches were found among 159,000 reads checked across the 106 files. |
| Bhavya | PRJNA610672 | 30 ENA FASTQs, one per run, about 0.55 GB. NCBI and the paper describe paired libraries. | Processed one file per run as supplied. The lengths and primer checks of two runs were consistent with premerged reads. The exact deposited preprocessing was unconfirmed. |
| Krupa | PRJEB20022 | 899 paired archive runs. NCBI SRA sample titles identify 859 subjects and 40 controls. | Used the forward reads of the 859 selected subject runs after the paired trial found no merges. The paper reports 897 healthy participants and describes a clinical panel. The reported output covers the selected archive runs. |
| Krupa | PRJNA668251 | 14 paired HiSeq runs, about 0.51 GB. | Used both mates and trimmed both expected primers. Both primer checks and a three-sample paired trial supported the configuration. |

The exploratory FASTQ checks used the first 1,500 reads per inspected file. They support starting settings and cannot establish homogeneity across every sample. The initial batch label `provisional_1` is a placeholder for each study. Use distinct batch labels if you confirm different sequencing runs or error profiles. The workflow learns a separate error model for each label. A `batch` is a physical sequencing batch and can contain many SRA run accessions.

The comparison reports within-sample bacterial phylum percentages using the same SILVA taxonomy reference. Study averages weight each retained sample equally. Study, primer region, processing history and cohort differ, so read-proportion differences are descriptive. Raw ASVs stay within the study for this comparison. Clinical groups require verified metadata and their own statistical design.

Sources

- PRJNA687346 [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA687346) and the paper *Association of Gut Microbiome and Vitamin D Deficiency in Knee Osteoarthritis Patients: A Pilot Study*.
- PRJNA761547 [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA761547) and the paper *A Sedentary Lifestyle Changes the Composition and Predicted Functions of the Gut Bacterial and Fungal Microbiota of Subjects from the Same Company*.
- PRJNA604466 [Belkova et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7479317/).
- PRJNA610672 [Thai cohort paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7409811/) and [ENA FASTQ report](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP251788&result=read_run).
- PRJEB20022 [Almonacid et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC5414997/), [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB20022) and [ENA study](https://www.ebi.ac.uk/ena/browser/view/PRJEB20022).
- PRJNA668251 [Vietnamese cohort paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7776958/).
- [DADA2 guidance on batch-specific error models and compatible amplicons](https://benjjneb.github.io/dada2/bigdata.html), and [paired-read workflow](https://benjjneb.github.io/dada2/bigdata_paired.html).
