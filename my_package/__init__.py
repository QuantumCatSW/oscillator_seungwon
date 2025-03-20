# my_package/__init__.py
# In legacy version of python, any directory without an __init__.py file is not recognized as a package.
# In modern Python, a folder can be recognized as a namespace package without an __init__.py, 
# but it can still cause complications when mixing older tools or certain packaging scenarios.

from .module import *
from .pulse_manager import *
