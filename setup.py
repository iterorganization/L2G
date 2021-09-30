from setuptools import setup, find_packages
from setuptools.extension import Extension
import sys
import os

try:
    from Cython.Build import cythonize
except ImportError:
    sys.exit("Cython not found. Cython is needed to build the extension "
             "modules.")

try:
    import numpy as np
except ImportError:
    sys.exit("Numpy not found. Numpy is needed to build the extension "
             "modules.")

if "L2G_CPP_ROOT_DIR" not in os.environ:
    sys.exit("L2G_CPP_ROOT_DIR is not set! L2G_cpp is required!")

if "EMBREE_ROOT_DIR" not in os.environ:
    sys.exit("EMBREE_ROOT_DIR is not set! Embree is required!")

useOpenMP = True

def getIncludeDirectories():
    """Required include directories:
        FLT
        embree3
    """
    out = []
    error = False
    error_msg = []
    # Get

    out.append(os.path.join(os.environ['L2G_CPP_ROOT_DIR'], 'include', 'flt'))
    out.append(os.path.join(os.environ['EMBREE_ROOT_DIR'], 'include'))

    out.append(np.get_include())

    return out

def getLibraries():
    """Required libraries to link:
        flt
    """
    out = []
    if sys.platform == "win32":
        out.append("flt")
        out.append("embree3")
    else:
        out.append("flt")
    return out

def getLibraryDirs():
    out = []

    if sys.platform == "win32":
        out.append(os.path.join(os.environ['L2G_CPP_ROOT_DIR'], 'lib'))
        out.append(os.path.join(os.environ['EMBREE_ROOT_DIR'], 'lib'))
    else:
        out.append(os.path.join(os.environ['L2G_CPP_ROOT_DIR'], 'lib'))

    return out

def getExtraCompileArgs():
    out = []

    if sys.platform == "win32":
        out.append("/O2")
        if useOpenMP:
            out.append("/openmp")

    else:
        out.append("-O2")
        out.append("-march=native")
        out.append("-Wall")
        if useOpenMP:
            out.append("-fopenmp")
    return out

def getExtraLinkArgs():
    out = []
    if sys.platform == "win32":
        pass
        #if useOpenMP:
        #    out.append("/openmp")
    else:
        if useOpenMP:
            out.append("-fopenmp")
    return out

def getDataFiles():
    out = []

    if sys.platform == "win32":
        # Gather the required DLLs
        out.append(os.path.join(os.environ['L2G_CPP_ROOT_DIR'], 'bin', 'flt.dll'))
        out.append(os.path.join(os.environ['EMBREE_ROOT_DIR'], 'bin', 'embree3.dll'))
        out.append(os.path.join(os.environ['EMBREE_ROOT_DIR'], 'bin', 'tbb12.dll'))

    return out

def getLongDescription():
    if os.path.exists("README"):
        return open("README", "r").read()
    return ""


cython_ext_module = Extension("L2G.core",
    sources=["L2G/core/flt.pyx"],
    include_dirs=getIncludeDirectories(),
    libraries=getLibraries(),
    library_dirs=getLibraryDirs(),
    extra_compile_args=getExtraCompileArgs(),
    extra_link_args=getExtraLinkArgs(),
    language="c++",
    # define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")], # Causes error with np.ndarray[].shape
)

cython_module = cythonize(cython_ext_module, annotate=True)

setup(
    name="L2G",
    version="1.0.0",
    description="Python module for running FLT",
    long_description=getLongDescription(),
    long_description_content_type="text/markdown",
    ext_modules=cython_module,
    author="Gregor Simic",
    author_email="simic.gregor@gmail.com",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Programming Language :: Cython",
        "Programming Language :: C++",
    ],
    packages = ["L2G"],
    data_files = [('', getDataFiles())],
    python_requires=">=3.6",
    scripts = ['bin/runL2G', 'bin/submitL2G', 'bin/L2G.sbatch',
               'bin/get_elm_data', 'bin/get_owl_conlen_graph']
)