.. _`design_builder`:

DESIGN_BUILDER
==============

Build design file from sample list

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 


Input/Output
------------


**Input:**

 
  
* A list of samples, one per line, stored in a file called "samples.txt".
  
 


**Output:**

 
  
* A TSV-formatted design file with Sample_id, Upstream_file, Downstream_file
  
 






Used wrappers
-------------

The following individual wrappers are used in this pipeline:


* :ref:`iRODS/search_samples`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.




Notes
-----

The samples must be in iRODS.

Under development. Do no use.




Authors
-------


* Thibault Dayris
