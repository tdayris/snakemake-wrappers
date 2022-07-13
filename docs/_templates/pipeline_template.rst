.. _`{{name}}`:

{{ name|upper }}
{{ name | length * '=' }}

{{ description }}

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash {% for cmd in usage %}

  {{ cmd | replace("\#", "#") }}{% endfor %}

{% if input and output %}
Input/Output
------------
{# Parse the input and output section of .yaml #}
{% for iotitle, io in ({"Input":input,"Output": output}).items() %}
**{{ iotitle }}:**

 {% for foo in io %}
  {% if foo is mapping %}
   {% for key, value in foo.items() %}
* ``{{ key }}``: {{ value }}
   {% endfor %}
  {% else %}
* {{ foo }}
  {% endif %}
 {% endfor %}

{% endfor %}
{% endif %}
{% if usedmetawrappers|length %}

Used meta-wrappers
------------------

The following individual meta-wrappers are used in this pipeline:

{% for umw in usedmetawrappers %}
* :ref:`{{ umw }}`
{% endfor %}

Please refer to each meta-wrapper in above list for additional configuration parameters and information about the executed code.

{% endif %}
{% if usedwrappers|length %}

Used wrappers
-------------

The following individual wrappers are used in this pipeline:

{% for uw in usedwrappers %}
* :ref:`{{ uw }}`
{% endfor %}

Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.

{% endif %}
{% if notes %}

Notes
-----

{{ notes }}
{% endif %}

{% if pipeline_snakefile_path %}

Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

{{ pipeline_snakefile_path }}

{% endif %}


Authors
-------

{% for author in authors %}
* {{ author }}
{% endfor %}
