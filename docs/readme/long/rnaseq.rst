
.. rnaseq_bulk:


RNA-Seq bulk
============

Summary
-------

This pipeline aims to perform classical RNASeq-bulk analyzes:


#. `TLDR <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#tldr-run-this-pipeline>`_
#. `design.tsv file <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#design-file>`_
#. `config.yaml file <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#configyaml>`_
#. `Classical use <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#classical-use>`_
#. `Quality controls <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#quality-controls>`_
#. `Quantification <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#quantification>`_
#. `Differential Gene Expression <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#differential-gene-expression>`_
#. `Aggregate factors <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#aggregate-factors>`_
#. `Ignore factors <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#ignore-factors>`_
#. `Perform only a subset of DGE <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#perform-only-a-subset-of-dge>`_
#. `Gene set enrichment analysis <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#gene-set-enrichment-analysis>`_
#. `Immune cell population deconvolution <>`_ Under construction
#. `Fusions <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#fusions>`_ Under construction
#. `Variant Calling <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#variant-calling>`_ Under construction

Each subsection can be selected one by one, or the whole pipeline may be executed.

TLDR: run this pipeline
-----------------------

.. code-block:: {sh}

   ## Go to your project directory
   cd /path/to/my/project/tmp

   ## Setup IO repositories
   ln -sfrv ../data_output data_output || mkdir -pv data_output
   ln -sfrv ../data_input data_input || mkdir -pv data_input

   ## Put design file in ${PWD}

   ## Run this pipeline
   bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh

Design file
-----------

The design file contains a description of your samples and experimental design. It should look like:

.. list-table::
   :header-rows: 1

   * - Sample_id
     - Upstream_file
     - Downstream_file
     - {ConditionA}
     - {ConditionB}
     - ...
   * - Sample1
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - BatchA
     - ...
   * - Sample2
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - BatchB
     - ...
   * - Sample3
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - BatchA
     - ...
   * - Sample4
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - BatchB
     - ...
   * - Sample5
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - BatchA
     - ...
   * - Sample6
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - BatchB
     - ...


``ConditionA`` and ``ConditionB`` are example names. Put the name of your own conditions instead of ConditionA and ConditionB. Add as many conditions as you with below. Keep in mind that *all* possible comparisons will be done (using ``~{condition_name}``\ , see below for accounting batch effects).

Header line **must** be the first line:


#. ``Sample_id`` **must** be in the header (first line). Column order does not matter.
#. ``Upstream_file`` **must** be in the header (first line). Column order does not matter.
#. ``Downstream_file`` **must** be in the header (first line). No single-ended analyzes here. Column order does not matter.
#. If ``BatchEffect`` is present in header (first line), and only if ``BatchEffect`` is present in header (first line), then the Differential Expression Analysis will take this column as batch effect (using ``~BatchEffect+{condition_name}``\ ). 

Only one batch effect can be considered. Do not call any of your conditions ``BatchEffect`` if you don't know what you are doing.

Each column will be considered as a factor. Each value in each cell will be considered as a level of its corresponding factor: 


#. A level with less than two replicates will be ignored. 
#. A factor with less than two levels will be ignored.

Config.yaml
-----------

This file describes all optional parameters set in the pipeline. 

Possible edition:


#. You may want to modify the list of genes of interest by adding a list after the provided example. 
#. You may change the path to your design file by changing the default ``design: design.tsv``.
#. You may want to analyze two factors as together (see: `aggregate factors <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#aggregate-factors>`_\ )
#. You may want to ignore factors present in your design file (see `ignore factors <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#ignore-factors>`_\ )
#. You may want to run a subset of possible DGE (see `subset DGE <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#perform-only-a-subset-of-dge>`_\ )

Do not change command line parameters if you do not know what you are doing. Do not set neither threading, temporary files, memory management, nor IO files in command line. This has already been done by the pipeline.

Changing anything else shall be done at your own risks.

Classical use
-------------

.. code-block:: {sh}

   ## Go to your project directory
   cd /path/to/my/project/tmp

   ## Setup IO repositories
   ln -sfrv ../data_output data_output || mkdir -pv data_output
   ln -sfrv ../data_input data_input || mkdir -pv data_input

   ## Put design file in ${PWD}

   ## Run basic quality controls, keep intermediary files
   bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh QC --nt

   ## Choose wether to continue the analysis or not.

   ## Edit config file
   ## 1. with your genes of interest
   ## 1. subset list of DGE
   ## 1. aggregate factors, etc.
   bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh quant --nt

   ## Have a look at mapping rates and basic quantification QC

   ## Run Differential Gene Expression
   bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh dge --nt

   ## Have a look at the DGE results

   ## Run fusions
   bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh fusions --nt

   ## Have a look at the results

   ## Clean TEMPORARY files and temporary files only
   bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh --delete-temp-output

Quality controls
----------------

Pipeline
^^^^^^^^


#. iRODS copy 
#. Fastp
#. FastqScreen
#. MultiQC

Use iRODS command to get your fastq files, then clean them with Fastp. Assess organism quality with FastqScreen, then aggregate quality reports with MultiQC.

Command line
^^^^^^^^^^^^

Command line argument order does not matter.

.. code-block:: {sh}

   ## Go to your project directory
   cd /path/to/my/project/tmp

   ## Setup IO repositories
   ln -sfrv ../data_output data_output || mkdir -pv data_output
   ln -sfrv ../data_input data_input || mkdir -pv data_input

   ## Put design file in ${PWD}

   ## GRCh38 / HG38
   bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh QC --nt
   ## GRCm38 / MM10
   bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh mm10 QC --nt

Results
^^^^^^^

The only output is a MultiQC report. By default, all other output are deleted, use ``--nt`` to keep cleaned fastq files.

.. code-block::

   data_output/
   ├── multiqc
       ├── MultiQC.QC_data
       ├── multiqc_data.json
       │   ├── multiqc_fastp.txt
       │   ├── multiqc_fastq_screen.txt
       │   ├── multiqc_general_stats.txt
       │   ├── multiqc.log
       │   └── multiqc_sources.txt
       └── MultiQC.QC.html

Quantification
--------------

Pipeline
^^^^^^^^


#. iRODS copy 
#. Fastp
#. FastqScreen
#. Salmon
#. MultiQC

Use iRODS command to get your fastq files, then clean them with Fastp. Assess organism quality with FastqScreen. Salmon estimates transcripts abundance, then aggregate quality reports with MultiQC. Salmon results are annotated with a basic GTF.

Command line
^^^^^^^^^^^^

Command line argument order does not matter.

.. code-block:: {sh}

   ## Go to your project directory
   cd /path/to/my/project/tmp

   ## Setup IO repositories
   ln -sfrv ../data_output data_output || mkdir -pv data_output
   ln -sfrv ../data_input data_input || mkdir -pv data_input

   ## Put design file in ${PWD}

   ## GRCh38 / HG38
   bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh quant --nt
   ## GRCm38 / MM10
   bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh mm10 quant --nt

Output
^^^^^^

The repository ``multiqc`` contains two multiqc reports: basic quality controls and salmon. MultiQC.Salmon.html contains all information present in MultiQC.QC.html, and adds the results of Salmon.

The repository ``Quantification`` contains three quantification files:


#. The file: ``Raw.genes.tsv`` contains the raw gene expression estimates, obtained from dummy summary of the expression of each transcripts counting for a same gene. This is **not** normalized. Useful for DESeq2, etc.
#. The file: ``TPM.genes.tsv`` contains TPM normalized gene expression estimates, obtained from dummy summary of the expression of each transcripts counting for a same gene. This is normalized by TPM. It does not take the factors and levels present in the ``design.tsv`` into account while normalizing.
#. The file: ``TPM.transcripts.tsv`` contains TPM normalized transcripts expression estimates, obtained Salmon. This is normalized by TPM. It does not take the factors and levels present in the ``design.tsv`` into account while normalizing.

.. code-block::

   data_output/
   ├── multiqc
   │   ├── MultiQC.QC_data
   │   │   ├── multiqc_data.json
   │   │   ├── multiqc_fastp.txt
   │   │   ├── multiqc_fastq_screen.txt
   │   │   ├── multiqc_general_stats.txt
   │   │   ├── multiqc.log
   │   │   └── multiqc_sources.txt
   │   ├── MultiQC.QC.html
   │   ├── MultiQC.Salmon_data
   │   │   ├── multiqc_data.json
   │   │   ├── multiqc_fastp.txt
   │   │   ├── multiqc_fastq_screen.txt
   │   │   ├── multiqc_general_stats.txt
   │   │   ├── multiqc.log
   │   │   └── multiqc_sources.txt
   │   └── MultiQC.Salmon.html
   └── Quantification
       ├── Raw.genes.tsv
       ├── TPM.genes.tsv
       └── TPM.transcripts.tsv

Differential Gene Expression
----------------------------

Pipeline
^^^^^^^^


#. iRODS copy 
#. Fastp
#. FastqScreen
#. Salmon
#. tximport
#. DESeq2
#. In-house scripts
#. MultiQC

Command line
^^^^^^^^^^^^

.. code-block:: {sh}

   ## Go to your project directory
   cd /path/to/my/project/tmp

   ## Setup IO repositories
   ln -sfrv ../data_output data_output || mkdir -pv data_output
   ln -sfrv ../data_input data_input || mkdir -pv data_input

   ## Put design file in ${PWD}

   ## GRCh38 / HG38
   bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh dge --nt
   ## GRCm38 / MM10
   bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh mm10 dge --nt

Results
^^^^^^^

The repository ``multiqc`` contains two multiqc reports: basic quality controls and salmon. MultiQC.Salmon.html contains all information present in MultiQC.QC.html, and adds the results of Salmon.

The repository ``Quantification`` contains three quantification files:


#. The file: ``Raw.genes.tsv`` contains the raw gene expression estimates, obtained from dummy summary of the expression of each transcripts counting for a same gene. This is **not** normalized. Useful for DESeq2, etc.
#. The file: ``TPM.genes.tsv`` contains TPM normalized gene expression estimates, obtained from dummy summary of the expression of each transcripts counting for a same gene. This is normalized by TPM. It does not take the factors and levels present in the ``design.tsv`` into account while normalizing.
#. The file: ``TPM.transcripts.tsv`` contains TPM normalized transcripts expression estimates, obtained Salmon. This is normalized by TPM. It does not take the factors and levels present in the ``design.tsv`` into account while normalizing.

The ``DEseq2`` repository contains one sub-folder for each comparison between two levels, for each factor. Aggregation, sub-setting and factor ignore may reduce of increase the number of sub-folders. See config file modification to know how to master the list of DESeq2 results. TODO

In each DESeq2 sub-folder:


#. The file ``Complete_XXX.tsv`` contains annotated results of DESeq2 with both differentially expressed genes under provided alpha threshold and the other ones. XXX being the name of the comparison.
#. The file ``SortedOnLogFC_XXX.tsv`` contains the differentially expressed genes only. They have been annotated and sorted on Log(FC). XXX being the name of the comparison.
#. The file ``SortedOnPadj_XXX.tsv`` contains the differentially expressed genes only. They have been annotated and sorted on adjusted P-Values. XXX being the name of the comparison.
#. The folder ``gene_plots`` contains a list of plots, one for each gene of interest. This is a highlight on its status, expression and per-condition weight.
#. The MultiQC report contains information about the samples belonging to the given comparison, and no other sample. It also contains Volcanoplot, PCAs, Expression weights, etc. Each of these plots is unique to this comparison, since expressions were normalized with DESeq2, which takes factors/levels into account.

.. code-block::

   data_output/
   ├── DEseq2
   │   ├── DGE_considering_factor_Batch_comparing_test_B1_vs_reference_B2
   │   │   ├── Complete_DGE_considering_factor_Batch_comparing_test_B1_vs_reference_B2.tsv
   │   │   ├── gene_plots
   │   │   │   └── ENSG00000141510.png
   │   │   ├── MultiQC.DEseq2_data
   │   │   │   ├── multiqc_data.json
   │   │   │   ├── multiqc_fastp.txt
   │   │   │   ├── multiqc_fastq_screen.txt
   │   │   │   ├── multiqc_general_stats.txt
   │   │   │   ├── multiqc.log
   │   │   │   └── multiqc_sources.txt
   │   │   ├── MultiQC.DEseq2.html
   │   │   ├── SortedOnLogFC_DGE_considering_factor_Batch_comparing_test_B1_vs_reference_B2.tsv
   │   │   └── SortedOnPadj_DGE_considering_factor_Batch_comparing_test_B1_vs_reference_B2.tsv
   ├── multiqc
   │   ├── MultiQC.QC_data
   │   │   ├── multiqc_data.json
   │   │   ├── multiqc_fastp.txt
   │   │   ├── multiqc_fastq_screen.txt
   │   │   ├── multiqc_general_stats.txt
   │   │   ├── multiqc.log
   │   │   └── multiqc_sources.txt
   │   ├── MultiQC.QC.html
   │   ├── MultiQC.Salmon_data
   │   │   ├── multiqc_data.json
   │   │   ├── multiqc_fastp.txt
   │   │   ├── multiqc_fastq_screen.txt
   │   │   ├── multiqc_general_stats.txt
   │   │   ├── multiqc.log
   │   │   └── multiqc_sources.txt
   │   └── MultiQC.Salmon.html
   └── Quantification
       ├── Raw.genes.tsv
       ├── TPM.genes.tsv
       └── TPM.transcripts.tsv

Aggregate factors
-----------------

Problem
^^^^^^^

Consider the following design:

.. list-table::
   :header-rows: 1

   * - Sample_id
     - Upstream_file
     - Downstream_file
     - Treatment
     - Status
   * - Sample1
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - Diseased
   * - Sample2
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - Relapse
   * - Sample3
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - Diseased
   * - Sample4
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - Relapse
   * - Sample5
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - Relapse
   * - Sample6
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - Diseased
   * - Sample7
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - Relapse
   * - Sample8
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - Diseased


It will produce the following comparisons:

.. code-block::

   DEseq2/
       ├── DGE_considering_factor_ConditionA_comparing_test_Untreated_vs_reference_Treated
       ├── DGE_considering_factor_ConditionB_comparing_test_Relapse_vs_reference_Diseased
       ├── DGE_considering_factor_ConditionA_comparing_test_Treated_vs_reference_Untreated
       └── DGE_considering_factor_ConditionB_comparing_test_Diseased_vs_reference_Relapse

Let us imagine we are interested in the effect of the treatment itself, AND the effect of the treatment on sample under "relapse" in relation to sample on relapse but without treatment.

Solution 1
^^^^^^^^^^

The first easy solution is to create an additional column, called as you wish (lets say Treatment_On_Relapse) and concatenate the two columns values. As a result, the design would be:

.. list-table::
   :header-rows: 1

   * - Sample_id
     - Upstream_file
     - Downstream_file
     - Treatment
     - Status
     - Treatment_On_Relapse
   * - Sample1
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - Diseased
     - Treated_Diseased
   * - Sample2
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - Relapse
     - Treated_Relapse
   * - Sample3
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - Diseased
     - Treated_Diseased
   * - Sample4
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - Relapse
     - Treated_Relapse
   * - Sample5
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - Relapse
     - Untreated_Relapse
   * - Sample6
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - Diseased
     - Untreated_Diseased
   * - Sample7
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - Relapse
     - Untreated_Relapse
   * - Sample8
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - Diseased
     - Untreated_Diseased


Now, the pipeline will produce the following results:

.. code-block::

   DEseq2/
       ├── DGE_considering_factor_ConditionA_comparing_test_Untreated_vs_reference_Treated
       ├── DGE_considering_factor_Treatment_On_Relapse_comparing_test_Treated_Relapse_vs_reference_Treated_Diseased
       ├── DGE_considering_factor_Treatment_On_Relapse_comparing_test_Untreated_Relapse_vs_reference_Treated_Diseased
       ├── DGE_considering_factor_Treatment_On_Relapse_comparing_test_Untreated_Diseased_vs_reference_Treated_Diseased
       ├── DGE_considering_factor_Treatment_On_Relapse_comparing_test_Untreated_Relapse_vs_reference_Untreated_Diseased
       ├── DGE_considering_factor_ConditionB_comparing_test_Relapse_vs_reference_Diseased
       ├── DGE_considering_factor_ConditionA_comparing_test_Treated_vs_reference_Untreated
       ├── DGE_considering_factor_Treatment_On_Relapse_comparing_test_Treated_Diseased_vs_reference_Treated_Relapse
       ├── DGE_considering_factor_Treatment_On_Relapse_comparing_test_Treated_Diseased_vs_reference_Untreated_Relapse
       ├── DGE_considering_factor_Treatment_On_Relapse_comparing_test_Treated_Diseased_vs_reference_Untreated_Diseased
       ├── DGE_considering_factor_Treatment_On_Relapse_comparing_test_Untreated_Diseased_vs_reference_Untreated_Relapse
       └── DGE_considering_factor_ConditionB_comparing_test_Diseased_vs_reference_Relapse

You do have the levels you are interested in.

Solution 2
^^^^^^^^^^

Within the ``config.yaml`` file, modify the value of 

.. code-block:: {yaml}

   deseq2: 
       design: 
           columns_to_aggregate: null

to:

.. code-block:: {yaml}

   deseq2: 
       design: 
           columns_to_aggregate:
               - - Treatment
                 - Status

It will have the very same consequences as described in `solution 1 <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#solution-1>`_ above. Please, note that the new factor will be named ``Treatment_Status``. You cannot rename factors with this solution.

Ignore factors
--------------

Problem
^^^^^^^

Sometimes, factors are not interesting anymore (not relevant in PCA, loss of interest by project research investigator, ...). We want this pipeline to ignore this factor. Consider the following ``design.tsv``\ :

.. list-table::
   :header-rows: 1

   * - Sample_id
     - Upstream_file
     - Downstream_file
     - Treatment
     - Status
   * - Sample1
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - Diseased
   * - Sample2
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - Relapse
   * - Sample3
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - Diseased
   * - Sample4
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - Relapse
   * - Sample5
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - Relapse
   * - Sample6
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - Diseased
   * - Sample7
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - Relapse
   * - Sample8
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - Diseased


Let's pretend we want to ignore the sample status Diseased/Relapse.

Solution 1
^^^^^^^^^^

Remove this column from your design. Then the pipeline won't be aware of it and will not consider this factor anymore.

Solution 2
^^^^^^^^^^

Within the ``config.yaml`` file, modify the value of 

.. code-block:: {yaml}

   deseq2: 
       design: 
           columns_to_ignore: null

to:

.. code-block:: {yaml}

   deseq2: 
       design: 
           columns_to_ignore:
               - Status

It will have the very same consequences as described in `solution 1 <https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#solution-1-1>`_ above.

Perform only a subset of DGE
----------------------------

Problem
^^^^^^^

Sometimes, we are not interested in all possible comparisons. For instance, the should be only one reference through all the comparisons. Let us consider the following ``design.tsv``\ :

.. list-table::
   :header-rows: 1

   * - Sample_id
     - Upstream_file
     - Downstream_file
     - Treatment
     - Status
   * - Sample1
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - Diseased
   * - Sample2
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - Relapse
   * - Sample3
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - Diseased
   * - Sample4
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Treated
     - Relapse
   * - Sample5
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - Relapse
   * - Sample6
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - Diseased
   * - Sample7
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - Relapse
   * - Sample8
     - /path/file.1.fq.gz
     - /path/file.2.fq.gz
     - Untreated
     - Diseased


We want "Untreated" level within "Treatment" factor to always be a reference level. We also want "Untreated_Relapse" to be our reference while comparing Treatment/Status interactions.

Solution
^^^^^^^^

Within the ``config.yaml`` file, modify the value of 

.. code-block:: {yaml}

   deseq2: 
       design: 
           columns_to_aggregate: null
           include_only: null

to:

.. code-block:: {yaml}

   deseq2: 
       design: 
           columns_to_aggregate:
               - - Treatment
                 - Status
           include_only:
               - _reference_Untreated_Relapse
               - _reference_Untreated

Now, the pipeline will produce the following results:

.. code-block::

   DEseq2/
       ├── DGE_considering_factor_ConditionA_comparing_test_Treated_vs_reference_Untreated
       └── DGE_considering_factor_Treatment_On_Relapse_comparing_test_Untreated_Diseased_vs_reference_Untreated_Relapse

Gene set enrichment analysis
----------------------------

Pipeline
^^^^^^^^

Pipeline
^^^^^^^^


#. iRODS copy 
#. Fastp
#. FastqScreen
#. Salmon
#. tximport
#. DESeq2
#. In-house scripts
#. ClusterProfiler
#. In-house scripts
#. MultiQC

Command line
^^^^^^^^^^^^

.. code-block:: {sh}

   ## Go to your project directory
   cd /path/to/my/project/tmp

   ## Setup IO repositories
   ln -sfrv ../data_output data_output || mkdir -pv data_output
   ln -sfrv ../data_input data_input || mkdir -pv data_input

   ## Put design file in ${PWD}

   ## GRCh38 / HG38
   bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh gsea --nt
   ## GRCm38 / MM10
   bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh mm10 gsea --nt

Results
^^^^^^^

The repository ``multiqc`` contains two multiqc reports: basic quality controls and salmon. MultiQC.Salmon.html contains all information present in MultiQC.QC.html, and adds the results of Salmon.

The repository ``Quantification`` contains three quantification files:


#. The file: ``Raw.genes.tsv`` contains the raw gene expression estimates, obtained from dummy summary of the expression of each transcripts counting for a same gene. This is **not** normalized. Useful for DESeq2, etc.
#. The file: ``TPM.genes.tsv`` contains TPM normalized gene expression estimates, obtained from dummy summary of the expression of each transcripts counting for a same gene. This is normalized by TPM. It does not take the factors and levels present in the ``design.tsv`` into account while normalizing.
#. The file: ``TPM.transcripts.tsv`` contains TPM normalized transcripts expression estimates, obtained Salmon. This is normalized by TPM. It does not take the factors and levels present in the ``design.tsv`` into account while normalizing.

The ``DEseq2`` repository contains one sub-folder for each comparison between two levels, for each factor. Aggregation, sub-setting and factor ignore may reduce of increase the number of sub-folders. See config file modification to know how to master the list of DESeq2 results. TODO

In each DESeq2 sub-folder:


#. The file ``Complete_XXX.tsv`` contains annotated results of DESeq2 with both differentially expressed genes under provided alpha threshold and the other ones. XXX being the name of the comparison.
#. The file ``SortedOnLogFC_XXX.tsv`` contains the differentially expressed genes only. They have been annotated and sorted on Log(FC). XXX being the name of the comparison.
#. The file ``SortedOnPadj_XXX.tsv`` contains the differentially expressed genes only. They have been annotated and sorted on adjusted P-Values. XXX being the name of the comparison.
#. The folder ``gene_plots`` contains a list of plots, one for each gene of interest. This is a highlight on its status, expression and per-condition weight.
#. The MultiQC report contains information about the samples belonging to the given comparison, and no other sample. It also contains Volcanoplot, PCAs, Expression weights, etc. Each of these plots is unique to this comparison, since expressions were normalized with DESeq2, which takes factors/levels into account.

.. code-block::

   data_output/
   ├── DEseq2
   │   ├── DGE_considering_factor_Batch_comparing_test_B1_vs_reference_B2
   │   │   ├── Complete_DGE_considering_factor_Batch_comparing_test_B1_vs_reference_B2.tsv
   │   │   ├── gene_plots
   │   │   │   └── ENSG00000141510.png
   │   │   ├── MultiQC.DEseq2_data
   │   │   │   ├── multiqc_data.json
   │   │   │   ├── multiqc_fastp.txt
   │   │   │   ├── multiqc_fastq_screen.txt
   │   │   │   ├── multiqc_general_stats.txt
   │   │   │   ├── multiqc.log
   │   │   │   └── multiqc_sources.txt
   │   │   ├── MultiQC.DEseq2.html
   │   │   ├── SortedOnLogFC_DGE_considering_factor_Batch_comparing_test_B1_vs_reference_B2.tsv
   │   │   └── SortedOnPadj_DGE_considering_factor_Batch_comparing_test_B1_vs_reference_B2.tsv
   ├── multiqc
   │   ├── MultiQC.QC_data
   │   │   ├── multiqc_data.json
   │   │   ├── multiqc_fastp.txt
   │   │   ├── multiqc_fastq_screen.txt
   │   │   ├── multiqc_general_stats.txt
   │   │   ├── multiqc.log
   │   │   └── multiqc_sources.txt
   │   ├── MultiQC.QC.html
   │   ├── MultiQC.Salmon_data
   │   │   ├── multiqc_data.json
   │   │   ├── multiqc_fastp.txt
   │   │   ├── multiqc_fastq_screen.txt
   │   │   ├── multiqc_general_stats.txt
   │   │   ├── multiqc.log
   │   │   └── multiqc_sources.txt
   │   └── MultiQC.Salmon.html
   └── Quantification
       ├── Raw.genes.tsv
       ├── TPM.genes.tsv
       └── TPM.transcripts.tsv

Fusions
-------

Pipeline
^^^^^^^^


#. iRODS copy (access iRODS collections and merge samples that were sequenced through multiple runs)
#. Fastp (trimm fastq reads)
#. STAR (mapping + fusions)
#. FusionAnnotator
#. FusionInspector
#. MultiQC

Command line
^^^^^^^^^^^^

.. code-block:: {sh}

   ## Go to your project directory
   cd /path/to/my/project/tmp

   ## Setup IO repositories
   ln -sfrv ../data_output data_output || mkdir -pv data_output
   ln -sfrv ../data_input data_input || mkdir -pv data_input
   ln -sfrv data_input/design.tsv design.tsv || echo "No design found. Create it, or let the pipeline guess sample pairs (risky)"

   ## Run basic quality controls, keep intermediary files, search, check and annotate fusions
   bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh fusions --nt

Variant Calling
---------------

Under construction
