.. _documentation-label:

.. currentmodule:: jacketse.jacket

.. _interfaces-label:

Documentation
-------------

If you have a normal Python code there are better ways to do this with autodoc.  For that case, see an example at `<https://raw.githubusercontent.com/WISDEM/CCBlade/master/docs/documentation.rst>`_.  For OpenMDAO classes there aren't any autodoc plugins (yet), so there isn't much you can do.  Fortunately, the definition of the assemblies can in many cases serve as useful documentaiton (assuming you have done a good job documenting through desc tags, units, etc.)  You can dump this out with a literalinclude

.. literalinclude:: ../src/jacketse/jacket.py
    :language: python
    :start-after: JacketSE(Assembly)
    :end-before: def configure(self)
    :prepend: class JacketSE(Assembly):

Beyond that you should provide links to important modules in your code.  This will allow you to refer to them in other modules and the source will be linked.


Referenced Utilization Modules
============================
.. module:: jacketse.Utilization
.. class:: TwrUtilization
.. class:: JcktUtilization
.. class:: UtilAssembly

Referenced loads Modules
============================
.. module:: jacketse.loads
.. class:: WindInputs
.. class:: WaterInputs
.. class:: LoadOutputs
.. class:: JcktLoadPre
.. class:: JcktLoadPost
.. class:: JcktLoad

Referenced Jacket Modules
============================

.. module:: jacketse.jacket
.. class:: Legs
.. class:: Xbraces
.. class:: HBraces
.. class:: MudBraces
.. class:: TP
.. class:: Tower
.. class:: Piles
.. class:: Soil
.. class:: SPIstiffness
.. class:: PileEmbdL

Referenced TowerSE Modules
====================================

.. module:: towerse.tower_supplement

Referenced CommonSE Modules
====================================

.. module:: commonse.Tube
.. module:: commonse.SoilC
.. module:: commonse.EmbedLength
.. module:: commonse.environment
.. module:: commonse.Material



