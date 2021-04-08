.. _`{{name}}`:

{{ name|upper }}
{{ name | length * '=' }}

{{ description }}

Usage
-----

In order to run the pipeline, use the following commands

.. list-table::
  :widths: 10 80
  :header-rows: 1
  :align: left

  .. code-block:: bash {% for cmd in usage %}

      {{ cmd | replace("\#", "#")}}{% endfor %}

{% if notes %}

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

Notes
-----

{{ notes }}
{% endif %}


Authors
-------

{% for author in authors %}
* {{ author }}
{% endfor %}
