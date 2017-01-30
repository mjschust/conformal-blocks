from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "Conformal Blocks",
    ext_modules = cythonize('fusion_prod/*.pyx', compiler_directives={'embedsignature': True}),  # accepts a glob pattern
)
