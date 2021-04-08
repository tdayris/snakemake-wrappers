.. _`EaCoN_Cytoscan`:

EACON_CYTOSCAN
==============

Analyse Cytoscans with EaCoN on Flamingo

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to your directory containing CEL files

  cd /path/to/CEL/dir

  # Copy/paster the following line

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/cytoscan_eacon/run.sh


Input/Output
------------


**Input:**

 
  
* CEL files
  
 


**Output:**

 
  
* Complete EaCoN anaysis directory
  
 
  
* HTML report
  
 





Used meta-wrappers
------------------

The following individual meta-wrappers are used in this pipeline:


* :ref:`meta/bio/eacon_cytoscan`

* :ref:`meta/bio/eacon_post_process`


Please refer to each meta-wrapper in above list for additional configuration parameters and information about the executed code.






Authors
-------


* Thibault Dayris

* Bastien Job

* Gérôme Jules-Clément
