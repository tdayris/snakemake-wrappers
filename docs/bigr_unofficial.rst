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
    touch bigr_pipelines/fastp_fastqscreen/{Snakefile,meta.yam,readme.md,run.sh}
    tree bigr_pipelines/fastp_fastqscreen


.. _snakemake_basics:


Snakemake basics with a simple rule
-----------------------------------

We need a fastp / fastq_screen pipeline. Let us start with fastp, the easier rule to write today. We need to
create a file in which to store the snakemake rule::

    touch bigr_pipelines/fastp_fastqscreen/rules/001.fastp.smk


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
     "My example is very cool"]

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


Resources reservation made easier
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Profile: Slurm made easier
^^^^^^^^^^^^^^^^^^^^^^^^^^

