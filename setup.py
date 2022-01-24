from setuptools import setup, find_packages
from setuptools.extension import Extension
import sys
import os

from typing import List
from pathlib import Path

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

def get_include_directories():
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

def get_libraries():
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

def get_library_dirs():
    out = []

    if sys.platform == "win32":
        out.append(os.path.join(os.environ['L2G_CPP_ROOT_DIR'], 'lib'))
        out.append(os.path.join(os.environ['EMBREE_ROOT_DIR'], 'lib'))
    else:
        out.append(os.path.join(os.environ['L2G_CPP_ROOT_DIR'], 'lib'))

    return out

def get_extra_compile_args():
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

def get_extra_link_args():
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

def get_pyx_files(path: Path) -> list:
    dir_exclusion = ["__pycache__"]
    pyx_files = []

    for p in path.glob("*"):
        if p.is_dir():
            if p.stem in dir_exclusion:
                continue
            pyx_files += get_pyx_files(p)
        if p.suffix == ".pyx":
            pyx_files.append(p)
    return pyx_files

def prepare_cython_extensions(pyx_files: List[Path], root_path):
    extensions = []
    for p in pyx_files:
        rel_path = p.relative_to(root_path)
        source_file_path = str(rel_path)
        # Correct module name so that it will not have src at the beginning.
        # Collect the parts of the relative path
        parts = rel_path.with_suffix("").parts
        # Remove the "src" and replace it with ModuleName
        module_name = ".".join(list(parts))
        print(f"Cythonizing {source_file_path} as module {module_name}...")
        # Create extensions
        ext = Extension(
            module_name,
            sources=[source_file_path],
            include_dirs=get_include_directories(),
            libraries=get_libraries(),
            library_dirs=get_library_dirs(),
            extra_compile_args=get_extra_compile_args(),
            extra_link_args=get_extra_link_args(),
            language="c++",
        )

        # Cythonize returns a list
        extensions += cythonize(ext, annotate=True)

    return extensions

# Collect the pyx files and create extensions from them
setup_path = Path(__file__).resolve().parent
package_path = setup_path / "l2g"
pyx_files = get_pyx_files(package_path)
extensions = prepare_cython_extensions(pyx_files, setup_path)

setup(
    name="l2g",
    version="1.0.0",
    description="Python module for running FLT",
    long_description=getLongDescription(),
    long_description_content_type="text/markdown",
    ext_modules=extensions,
    author="Gregor Simic",
    author_email="simic.gregor@gmail.com",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Programming Language :: Cython",
        "Programming Language :: C++",
    ],
    packages = ["l2g", "l2g.comp", "l2g.comp.core", "l2g.equil", "l2g.hlm",
                "l2g.hlm", "l2g.plot",  "l2g.utils"],
    data_files = [('', getDataFiles())],
    python_requires=">=3.6",
    scripts = ['bin/runL2G', 'bin/submitL2G', 'bin/L2G.sbatch',
               'bin/get_elm_data', 'bin/get_owl_conlen_graph',
               'bin/apply_elm_hlm', 'bin/apply_exp_hlm',
               'bin/imas2eqdsk', 'bin/mkEqdskMovie', 'bin/mkImasMovie']
)