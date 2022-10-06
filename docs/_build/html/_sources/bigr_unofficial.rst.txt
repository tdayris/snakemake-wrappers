.. _bigr_unofficial:

BiGR Unofficial
===============


Why BiGR Unofficial snakemake Wrappers ?
----------------------------------------

Well, we all use the same cluster at BiGR. We all follow the same project template, all use the same gitlab 
instance, all use the same pipelines. We all share knowledge and code snippets. Wouldn't it be nice if we also
has a single entry point for our developments projects ?

Sure.

Now let's create a new pipeline together. We need a simple fastp + fastq screen pipeline. Let us list the
requirements of such pipeline:

#. Single-command launch
#. This pipeline expects a R1 and a R2 file.
#. Default parameters must be fine with Human genome, and oncology studies.
#. If user provided fastq files in working directories, search them, use them.
#. If user provided a list of fastq files in a file called ``design.tsv``, use this file instead of searching for fastq.
#. If user provided optional parameters through a file called ``config.yaml``, use this file. Else, create one with default parameters.
#. Input datasets must be in a directory called "data_input" and should be deleted once they are not useful anymore.
#. Output reports and datasets must be in a directory called "data_output".
#. Temporary files must be in a directory called "tmp" or deleted once they are not useful anymore.
#. Pipeline must tell user wether it failed, it is still running, or it has finished.
#. Pipeline must install itself and must be re-runnable with same results.

All of these points are made easier on Flamingo with BiGR Unofficial snakemake wrappers, let's see how and why.

BiGR pipeline template
----------------------

We need a working environment. Let's create a working branch::
    
    cd /mnt/beegfs/userdata/t_dayris/
    git clone https://github.com/tdayris/snakemake-wrappers.git
    cd snakemake-wrappers
    git checkout Unofficial
    git checkout -B generic_fastp_fastq_screen_pipeline
    source bigr_pipelines/bash_aliases.sh


Let's now create our working environment::
    
    mkdir -pv bigr_pipelines/fastp_fastqscreen/{rules,script,config}
    touch bigr_pipelines/fastp_fastqscreen/{Snakefile,meta.yaml,readme.md,run.sh}
    tree bigr_pipelines/fastp_fastqscreen


This gives::

    bigr_pipelines/fastp_fastqscreen/
    ├── Snakefile      # <- This is the Snakemake entry point
    ├── meta.yaml      # <- Short codified help
    ├── readme.md      # <- The long help. Feel free to write anything you want
    ├── run.sh         # <- This is our bash entry point
    ├── rules          # <- There will be our snakemake rules
    ├── script         # <- Your in-house scripts
    └── config         # <- Default configuration for human/mouse/...

.. _snakemake_basics:


Snakemake basics with a simple rule
-----------------------------------

We need a fastp / fastq_screen pipeline. Let us start with fastp, the easier rule to write today. We need to
create a file in which to store the snakemake rule::

    touch bigr_pipelines/fastp_fastqscreen/rules/001.fastp.smk
    tree bigr_pipelines/fastp_fastqscreen/

    bigr_pipelines/fastp_fastqscreen/
    ├── Snakefile
    ├── meta.yaml
    ├── readme.md
    ├── run.sh
    ├── rules
    │   └── 001.fastp.smk
    ├── script
    └── config


I called it ``001.fastp.smk`` because it's the first rule in out pipeline, and the file will contain fastp only.
I feel like it is a good practice to name your snakefile based on their content, thus, for future developers,
it will help them to find the snakemake rules "in oder of execution". At least, coworkers asked for that.

Our snakefile `bigr_pipelines/fastp_fastqscreen/rules/001.fastp.smk` contains the following::

    rule fastp:          # The name of the rule, do not start with number. Keep it clear.
        input:           # The block that contains Input files
            sample=      # List of strings. Paths to fastq files
        output:          # The block that contains Output files
            trimmed=     # List of strings. Paths to trimmed fastq files
            failed=      # List of strings. Paths to discarded reads
            html=        # String. Path to html formatted trimming report.
            json=        # String. Path to json formatted trimming report.
        threads:         # Number of threads used
        resources:       # The block that contains cluster resources
            mem_mb=      # Integer or Function, the number of MB or ram needed.
            time_min=    # Integer or Function, the time (minutes) required to run fastp
            tmp=         # Path to overloaded temporary directory
        log:             # Path to logging file
        params:          # The block that contains non-file parameters
            extra=       # Optional parameters for fastp
            adapters=    # Optional adapters list for fastp
        wrapper:         # Path/url to the used wrapper

How in the hell did I come up with all of that ?

Well, I'm using `fastp <https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/fastp.html>`_ and there is already 
a documentation, installation, and running examples of Fastp with the Snakemake wrappers. No need for me to search
for installation paths, nor deployments scheme. Everything is shipped with the wrappers. I copied the example, and
now we will fill it together.

The input block
^^^^^^^^^^^^^^^

Let's fill the input block::

        sample=expand(
            "data_input/{sample}.{stream}.fq.gz",   # Input are in data_input
            sample="{sample}",                      # Regular expression for sample names
            stream=["1", "2"]                       # List of sample names
        )


Snakemake has a function called `expand <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#the-expand-function>`_. To make it simple, it returns an iterable (like a list) wich contains all
possible combinations of strings between "{}". If I write::

    expand("My example is {adjective} {status}", adjective=['very', 'truly'], status=['awesome', 'cool'])

Then it will result in::

    ["My example is very awesome", 
     "My example is truly awesome", 
     "My example is very cool", 
     "My example is truly cool"]

Got it ? Good.

In our case, the expand function contains a regular expression: ``{sample}``. In snakemake, this is called a 
`wildcard <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards>`_ . We would like this
regular expression to match sample names, and sample names only. We can constrain it. We will dot that later,
be confident, and let a TODO marker on your desktop.

For now, the expand function we wrote in input of our rule fastp will result in::

    [r"data_input/.+.1.fq.gz", r"data_input/.+.2.fq.gz"]

These are two regular expressions, not that, for now, "{sample}" is solved as ".+", but we shall precise it later.

The output block
^^^^^^^^^^^^^^^^

Let's fill the output block::

    output:
        trimmed=expand(  # We expect two trimmed fastq files, 1 for R1, 1 for R2.
            "001.fastp/trimmed/{sample}.{stream}.fastq",
            sample="{sample}",
            stream=["1", "2"]
        ),
        # failed= We are not interested in discarded reads, this parameter is flagged optional in wrapper doc.
        html="001.fastp/html/{sample}.html",  # We want only one report per sample
        json="001.fastp/json/{sample}.json",  # We want only one report per sample


It is very important that all output share the same wildcards as input. Otherwise, Snakemake won't be able to solve
the dependency graph and won't be able to know what to do. This is the same with GNUMake, CMake, BioMake, etc.

That one was easy, wasn't it ?

Threads and Resources
^^^^^^^^^^^^^^^^^^^^^

These two blocks are here for cluster execution. The name and values available in it are *fixed and required*. Please
always fill them. Otherwise default values will be provided, but my not match with your needs.

The wrapper documentation says Fastp allows multiple threads.

Usually, Fastp requires quite few RAM and does not take long to run.

Let's fill these blocks::

    threads: 4
    resources:
        mem_mb=1024 * 2, # Two GB should be enough
        time_min=20,     # 20 minutes should be enough
        tmp="tmp",       # "tmp" not "/tmp" ! This is very important.


If you ever analyse a very big pair of fastq files, with hundreads of milions of reads, then Fastp fill take longer
than 20 minutes to run and the cluster will delete your job. This is not cool. We would like to automatically
re-launch the rule with more time. Or even better, reserve time according to the fastq file size !

Both are possible, but we will come to that a bit later. This is still Snakemake basics ! Don't ask such interesting 
questions so soon !


Parameters and logging
^^^^^^^^^^^^^^^^^^^^^^

Let's fill these blocks::

    log:
        "logs/fastp/{sample}.log"
    params:
        extra=""      # Optional parameters, see fastp CLI documentation
        # adapters="" # Optional according to the fastp wrapper documentation.


It is very important that all output and lgos share the same wildcards as input. Otherwise, Snakemake won't be 
able to solve the dependency graph and won't be able to know what to do. This is the same with GNUMake, 
CMake, BioMake, etc.

Parameters are almost unique for each wrapper, please refer to the wrapper documentation.


Link to the wrapper
^^^^^^^^^^^^^^^^^^^

The easiest way to link to the wrapper is to put::

    wrapper:
        "bio/fastp"


How do I come up with that ? Well, in the repository `snakemake-wrappers`, fastp is under: ``bio/fastp``. Yes, it's
that simple.

This will tell Snakemake how to install and deploy Fastp, how to build the command line and how to launch this tool.

Thanks to the Unofficial environment, Snakemake won't search on the web but on Flamingo, and if the environment is
already available, then Snakemake won't re-install it. Environments are sha-signed, so you cannot destroy other
people's environments, you can only share them ! Enjoy !


Conclusion
^^^^^^^^^^

Our file ``rules/001.fastp.smk`` contains::
    

    rule fastp:
        input:
            sample=expand(
                "data_input/{sample}.{stream}.fq.gz",   # Input are in data_input
                sample="{sample}",                      # Regular expression for sample names
                stream=["1", "2"]                       # List of sample names
            )
        output:
            trimmed=expand(  # We expect two trimmed fastq files, 1 for R1, 1 for R2.
                "001.fastp/trimmed/{sample}.{stream}.fastq",
                sample="{sample}",
                stream=["1", "2"]
            ),
            # failed= We are not interested in discarded reads, this parameter is flagged optional in wrapper doc.
            html="001.fastp/html/{sample}.html",  # We want only one report per sample
            json="001.fastp/json/{sample}.json",  # We want only one report per sample
        threads: 4
        resources:
            mem_mb=1024 * 2, # Two GB should be enough
            time_min=20,     # 20 minutes should be enough
            tmp="tmp",       # "tmp" not "/tmp" ! This is very important.
        log:
            "logs/fastp/{sample}.log"
        params:
            extra=""      # Optional parameters, see fastp CLI documentation
            # adapters="" # Optional according to the fastp wrapper documentation.
        wrapper:
            "bio/fastp"



.. _snakemake_advanced:

Snakemake level up with Unofficial environment
----------------------------------------------

We saw earlier that wildcards were equal to ".+" by default and might be constrained. We saw earlier that resources
could be set upon input file size. We saw earlier that Snakemake can re-run rules that failed due to out-of-memory
error, or out-of-time error.

We said we would talk about it in the future.

The. Future. Is. Now.

Also, the future is written in: ``000.commons.smk``. It's the first file included in the main Snakefile,
and it contains only pure python functions and intructions.


Constrain Wildcards
^^^^^^^^^^^^^^^^^^^

Remember, ``"{samlple}"`` is solved as ``".+"`` (any character, at least once). We would like to have the list of samples
instead.

So, we need to acquire the list of samples, either from user input, of though disc search.

Unofficial provides the functions for that. Here, we search for pairs fastq files. The snakemake-unofficial API lists a 
function named ``search_fastq_pairs`` in the ``file_manager`` module. This function returns a dictionary formatted 
as follows::

    {Sample1: {Upstream_file: /path/to/R1.fq.gz, Downstream_file: /path/to/R2.fq.gz}, Sample2 ...}


This can be passed to the function ``get_design`` in the module ``file_manager`` to produce a design file.

As easy as is sound, we just need to import a single module, use two function and *BAM*! A design file and a
list of available pairs of fastq file per sample is available !

So, in `000.commons.smk`, let's write the following:



    # Manage paths easily in Python
    from pathlib import Path

    # Get workflow path (the one you're working on!)
    workflow_source_dir = Path(snakemake.workflow.srcdir(".."))

    # Get path to shared libraries in Unofficial Snakemake wrappers
    common = str(workflow_source_dir / ".." / "common" / "python")

    # Add this to the list of available libraries in Python
    import sys
    sys.path.append(common)

    # Finally ! Add the module file_manager
    from file_manager import *


The search of all fastq files and the build of the design file is made by the following::

    import os # To deal with operating system in Python
    design = get_design(
        dirpath=os.getcwd(),           # Where to search for fastq files
        search_func=search_fastq_pairs # What function to use in order to search for fastq files
    )


This saves on disk the result, if and only if another design *does not* exist. No file will be deleted.

The list of available sample files is available with::

    print(design.Sample_id)


We want to constrain our wildcards, this is easily done with the following command::

    wildcard_constraints:
        sample=r"|".join(design.Sample_id),
        stream=r"1|2",


Tadaaaam ! Our pipeline now finds fastq files by itself. By the way, if your sample name ends
with a number, Snakemake won't confuse it with sequencing strands. 

Let's imagine we have two samples called: `Sample1` and `Sample2`. The ``expand`` in `001.fastp.smk`
contains::

    sample=expand(
        "data_input/{sample}.{stream}.fq.gz",   # Input are in data_input
        sample="{sample}",                      # Regular expression for sample names
        stream=["1", "2"]                       # List of sample names
    )


and solves itself as::

    ["data_input/Sample1|Sample2.1.fq.gz", "data_input/Sample1|Sample2.2.fq.gz"]


No ambiguity is allowed in Snakemake.


Resources reservation made easier
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We said earlier that resources could be set upon input file size. These resources should be adjusted
in case we need more memory and/or time for our job on the cluster.

This can be done with functions respecting the `signature <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources>`_ 
described in official snakemake documentation. Unofficial Snakemake Wrapper has a set of function designed 
to make it easy and human readable. These functions are in ``reservation`` module. Just like before, they 
will be imported in your `000.commons.smk`::

    # Manage paths easily in Python
    from pathlib import Path

    # Get workflow path (the one you're working on!)
    workflow_source_dir = Path(snakemake.workflow.srcdir(".."))

    # Get path to shared libraries in Unofficial Snakemake wrappers
    common = str(workflow_source_dir / ".." / "common" / "python")

    # Add this to the list of available libraries in Python
    import sys
    sys.path.append(common)

    # Finally ! Add the module file_manager
    from file_manager import *

    # Import the functions for easy time/memory reservation
    from reservation import *

    import os # To deal with operating system in Python
    design = get_design(
        dirpath=os.getcwd(),           # Where to search for fastq files
        search_func=search_fastq_pairs # What function to use in order to search for fastq files
    )

And your ``fastp`` rule goes from::

    resources:
        mem_mb=1024 * 2, # Two GB should be enough
        time_min=20,     # 20 minutes should be enough
        tmp="tmp",       # "tmp" not "/tmp" ! This is very important.


to::

    resources:
        mem_mb=get_2gb_per_attempt,   # 2gb on first try, then 4gb, then 6gb, ...
        time_min=get_20m_per_attempt, # 20min on first try, then 40min, 1h, ...
        tmp="tmp",


Unofficial Snakemake Wrappers will handle the re-submission.


Configuration and user-defined parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Depending on a project content, on the biological question, or simply based on the organism, many
tools and command line may differ. For instance, read adapters may be provided depending on the
sequencing library that was designed. Trimming parameters differ from bulk-rnaseq and whole exome
sequencing.

These information will be stored in a configuration (yaml formatted) that your fellow bioinformatician
may change (or not). These parameters must be forewarded while executing the pipeline instructions.

Let us create the first configuration file::

    touch bigr_pipelines/fastp_fastqscreen/config/config.yaml
    tree bigr_pipelines/fastp_fastqscreen/

    bigr_pipelines/fastp_fastqscreen/
    ├── Snakefile
    ├── meta.yaml
    ├── readme.md
    ├── run.sh
    ├── rules
    │   ├── 000.common.smk
    │   └── 001.fastp.smk
    ├── script
    └── config
        └── config.yaml


Let us define some variables in this configuragion file::

    # Maximum number of threads
    max_threads: 20

    fastp:
        extra: "" # Optional parameters
        adapters: null # Optional adapters


This file is stored in the pipeline source. We cannot allow the user to modify it. We must provide
a copy of the default parameters for the user to change them.

Once again, Unofficial Snakemake Wrappers provide a function for that. In the file `000.common.smk`
let's write::

    # Manage paths easily in Python
    from pathlib import Path

    # Get workflow path (the one you're working on!)
    workflow_source_dir = Path(snakemake.workflow.srcdir(".."))

    # Get path to shared libraries in Unofficial Snakemake wrappers
    common = str(workflow_source_dir / ".." / "common" / "python")

    # Add this to the list of available libraries in Python
    import sys
    sys.path.append(common)

    # Finally ! Add the module file_manager
    from file_manager import *

    # Import the functions for easy time/memory reservation
    from reservation import *

    # Import Config file if needed, use user-defined one if present.
    default_config = read_yaml(workflow_source_dir / "config" / "config.yaml")
    configfile: get_config(default_config=default_config)

    import os # To deal with operating system in Python
    design = get_design(
        dirpath=os.getcwd(),           # Where to search for fastq files
        search_func=search_fastq_pairs # What function to use in order to search for fastq files
    )


We can now modify our file `001.fastp.smk` in order to include configurations::

    rule fastp:
        input:
            sample=expand(
                "data_input/{sample}.{stream}.fq.gz",   # Input are in data_input
                sample="{sample}",                      # Regular expression for sample names
                stream=["1", "2"]                       # List of sample names
            )
        output:
            trimmed=expand(  # We expect two trimmed fastq files, 1 for R1, 1 for R2.
                "001.fastp/trimmed/{sample}.{stream}.fastq",
                sample="{sample}",
                stream=["1", "2"]
            ),
            # failed= We are not interested in discarded reads, this parameter is flagged optional in wrapper doc.
            html="001.fastp/html/{sample}.html",  # We want only one report per sample
            json="001.fastp/json/{sample}.json",  # We want only one report per sample
        threads: 4
        resources:
            mem_mb=1024 * 2, # Two GB should be enough
            time_min=20,     # 20 minutes should be enough
            tmp="tmp",       # "tmp" not "/tmp" ! This is very important.
        log:
            "logs/fastp/{sample}.log"
        params:
            # If user deleted 'extra' in fastp, then leave it empty and do not raise error
            extra=config["fastp"].get('extra', '')
            # If user deleted 'adapters' in fastp, then leave it to `None` and do not raise error
            adapters=config["fastp"].get('extra', None)
        wrapper:
            "bio/fastp"


Fastq Screen the final example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fastq screen requires a new file in the repository ``rules``, we name it ``002.fastq_screen.smk``::

    touch bigr_pipelines/fastp_fastqscreen/rules/002.fastq_screen.smk
    tree bigr_pipelines/fastp_fastqscreen/

    bigr_pipelines/fastp_fastqscreen/
    ├── Snakefile
    ├── meta.yaml
    ├── readme.md
    ├── run.sh
    ├── rules
    │   ├── 000.common.smk
    │   ├── 001.fastp.smk
    │   └── 002.fastq_screen.smk
    ├── script
    └── config
        └── config.yaml


Thanks to the `wrapper documentation <https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/fastq_screen.html>`_, 
we know that the file ``rules/002.fastq_screen.smk`` contains::

    rule fastq_screen:
        input:
            "data_input/{sample}.{stream}.fq.gz"
        output:
            txt=temp("fastq_screen/{sample}.fastq_screen.txt"),
            png=temp("fastq_screen/{sample}.fastq_screen.png")
        threads: config.get("threads", 20)
        resources:
            mem_mb=get_15gb_per_attempt,
            time_min=get_1h_per_attempt,
            tmpdir="tmp"
        params:
            fastq_screen_config=config["fastq_screen"],
            subset=100000,     # We do not allow user to change the default number of mapped reads
            aligner='bowtie2'  # We do not allow user to change default mapper
        log:
            "logs/fastq_screen/{sample}.log"
        wrapper:
            "bio/fastq_screen"


In the file ``config/config.yaml`` we add the databases::

    # Maximum number of threads
    max_threads: 20

    # Fastp parameters
    fastp:
        extra: "" # Optional parameters
        adapters: null # Optional adapters
    
    # Fastq Screen info
    fastq_screen:
        database:
            Contaminants:
            bowtie2: /mnt/beegfs/database/bioinfo/Index_DB/Fastq_Screen/0.14.0/Adapters/Contaminants
            AThaliana:
            bowtie2: /mnt/beegfs/database/bioinfo/Index_DB/Fastq_Screen/0.14.0/Arabidopsis/Arabidopsis_thaliana.TAIR10
            Drosophila:
            bowtie2: /mnt/beegfs/database/bioinfo/Index_DB/Fastq_Screen/0.14.0/Drosophila/BDGP6
            Ecoli:
            bowtie2: /mnt/beegfs/database/bioinfo/Index_DB/Fastq_Screen/0.14.0/E_coli/Ecoli
            Human:
            bowtie2: /mnt/beegfs/database/bioinfo/Index_DB/Fastq_Screen/0.14.0/Human/Homo_sapiens.GRCh38
            Lambda:
            bowtie2: /mnt/beegfs/database/bioinfo/Index_DB/Fastq_Screen/0.14.0/Lambda/Lambda
            Mitochondria:
            bowtie2: /mnt/beegfs/database/bioinfo/Index_DB/Fastq_Screen/0.14.0/Mitochondria/mitochondria
            Mouse:
            bowtie2: /mnt/beegfs/database/bioinfo/Index_DB/Fastq_Screen/0.14.0/Mouse/Mus_musculus.GRCm38
            PhiX:
            bowtie2: /mnt/beegfs/database/bioinfo/Index_DB/Fastq_Screen/0.14.0/PhiX/phi_plus_SNPs
            Rat:
            bowtie2: /mnt/beegfs/database/bioinfo/Index_DB/Fastq_Screen/0.14.0/Rat/Rnor_6.0
            rRNA:
            bowtie2: /mnt/beegfs/database/bioinfo/Index_DB/Fastq_Screen/0.14.0/rRNA/GRCm38_rRNA
            Vectors:
            bowtie2: /mnt/beegfs/database/bioinfo/Index_DB/Fastq_Screen/0.14.0/Vectors/Vectors
            Worm:
            bowtie2: /mnt/beegfs/database/bioinfo/Index_DB/Fastq_Screen/0.14.0/Worm/Caenorhabditis_elegans.WBcel235
            Yeast:
            bowtie2: /mnt/beegfs/database/bioinfo/Index_DB/Fastq_Screen/0.14.0/Yeast/Saccharomyces_cerevisiae.R64-1-1
            SalmoSalar:
            bowtie2: /mnt/beegfs/database/bioinfo/Index_DB/Fastq_Screen/0.14.0/SalmoSalar/Salmo_salar.ICSASG_v2
        aligner_paths:
            bowtie2: bowtie2


The file `000.commons.smk` does not need to be changed.
