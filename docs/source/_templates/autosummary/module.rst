{{ fullname }}
{{ "=" * fullname|length }}

.. automodule:: {{ fullname }}
   :members:
   :undoc-members:
   :show-inheritance:

{% if modules %}
Submodules
----------

.. autosummary::
   :toctree: .
   :recursive:

{% for item in modules %}
   {{ item }}
{% endfor %}
{% endif %}

{% if classes %}
Classes
-------

.. autosummary::
   :toctree: .

{% for item in classes %}
   {{ item }}
{% endfor %}
{% endif %}

{% if functions %}
Functions
---------

.. autosummary::
   :toctree: .

{% for item in functions %}
   {{ item }}
{% endfor %}
{% endif %}

{% if exceptions %}
Exceptions
----------

.. autosummary::
   :toctree: .

{% for item in exceptions %}
   {{ item }}
{% endfor %}
{% endif %}
