from setuptools import setup
from setuptools.extension import Extension
import sys
import os

from pathlib import Path

try:
    from Cython.Build import cythonize
except ImportError:
    sys.exit("Cython not found. Cython is needed to build the extension " +
             "modules.")

try:
    import numpy as np
except ImportError:
    sys.exit("Numpy not found. Numpy is needed to build the extension " +
             "modules.")

# Find L2G_CPP_ROOTDIR
L2G_CPP_ROOTDIR: str
EMBREE_VERSION: str
EMBREE_ROOTDIR: str

if (L2G_CPP_ROOTDIR := os.environ.get("L2G_CPP_ROOT_DIR") or os.environ.get("EBROOTL2G_CPP")) is None:
    sys.exit("L2G_CPP_ROOT_DIR or EBROOTL2G_CPP are not set! L2G_cpp is required!")

if (EMBREE_ROOTDIR := os.environ.get("EMBREE_ROOT_DIR") or os.environ.get("EBROOTEMBREE")) is None:
    sys.exit("EMBREE_ROOT_DIR or EBROOTEMBREE are not set! Embree is required!")

if (EMBREE_VERSION:=os.environ.get("EMBREE_VERSION") or os.environ.get("EBVERSIONEMBREE")) is None:
    # Try to get it from EBVERSIONEMBREE
    sys.exit("EMBREE_VERSION is not set in the envrionment!")

if len(EMBREE_VERSION) > 1:
    EMBREE_VERSION = EMBREE_VERSION[0]

useOpenMP = True

def get_include_directories() -> list[str]:
    """Get the include directories of the C++ code and Embree.
    """
    out: list[str] = []

    out.append(os.path.join(L2G_CPP_ROOTDIR, 'include'))
    out.append(os.path.join(EMBREE_ROOTDIR, 'include'))
    out.append(np.get_include())

    return out

def get_libraries() -> list[str]:
    """Get libraries required for linking against flt.so
    """
    out: list[str] = []
    out.append("flt")
    if sys.platform == "win32":
        out.append("embree3")
    return out

def get_library_dirs() -> list[str]:
    """Get paths to library dirs.
    """
    out: list[str] = []

    out.append(os.path.join(L2G_CPP_ROOTDIR, 'lib'))
    if sys.platform == "win32":
        out.append(os.path.join(EMBREE_ROOTDIR, 'lib'))

    return out

def get_extra_compile_args() -> list[str]:
    """Specify extra compilation arguments. Use maximum optimization but avoid
    using march=native as it shown that on some CPUs it can generate segfaulty
    code.
    """
    out: list[str] = []

    if sys.platform == "win32":
        out.append("/O2")
        if useOpenMP:
            out.append("/openmp")
        # Conversion warning
        out.append("/wd4244")
        # dll-interface warning
        out.append("/wd4251")
        # The C++ code is compatible with Embree 3/4
        out.append(f"/DEMBREE_VERSION={EMBREE_VERSION}")
    else:
        out.append("-O3")
        out.append("-march=native")
        out.append("-Wall")
        if useOpenMP:
            out.append("-fopenmp")
        # The C++ code is compatible with Embree 3/4
        out.append(f"-DEMBREE_VERSION={EMBREE_VERSION}")
    return out

def get_extra_link_args() -> list[str]:
    """Extra link flags. For OpenMP or parts of the code in cython that might
    use openmp.
    """
    out: list[str] = []

    if useOpenMP:
        if sys.platform == "win32":
           # out.append("/openmp")
           pass
        else:
            out.append("-fopenmp")
    return out

def getDataFiles() -> list[str]:
    """Get a list of data files to package inside the wheel. This is important
    on Microsoft sstems where the dll files should be bundled next to the 
    cythonized dll files.

    """
    out: list[str] = []

    # if sys.platform == "win32":
    #     # Gather the required DLLs
    #     out.append(os.path.join(L2G_CPP_ROOTDIR, 'bin', 'flt.dll'))
    #     out.append(os.path.join(EMBREE_ROOTDIR, 'bin', 'embree3.dll'))
    #     out.append(os.path.join(EMBREE_ROOTDIR, 'bin', 'tbb12.dll'))

    return out

def getLongDescription() -> str:
    if os.path.exists("README"):
        return open("README", "r").read()
    return ""

def get_pyx_files(path: Path) -> list[str]:
    """Collect all the pyx files
    """
    dir_exclusion = ["__pycache__"]
    pyx_files: list[str] = []

    for p in path.glob("*"):
        if p.is_dir():
            if p.stem in dir_exclusion:
                continue
            pyx_files += get_pyx_files(p)
        if p.suffix == ".pyx":
            pyx_files.append(p)
    return pyx_files

def prepare_cython_extensions(pyx_files: list[Path], root_path: Path) -> list[Extension]:
    """Automatically setup the list of extensions. Basically every pyx file in
    the package path directory get's properly defined as an Extension and then
    cythonized.
    """
    extensions: list[Extension] = []
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

        # Cythonize returns a list.
        # Add embedsignature so that the functions are visible by Sphinx.
        extensions += cythonize(ext, annotate=True,
                                compiler_directives={'embedsignature': True})

    return extensions

# Collect the pyx files and create extensions from them
setup_path = Path(__file__).resolve().parent
package_path = setup_path / "l2g"
pyx_files = get_pyx_files(package_path)
extensions = prepare_cython_extensions(pyx_files, setup_path)

setup(
    ext_modules=extensions,
    packages = ["l2g", "l2g.comp", "l2g.equil", "l2g.hlm", "l2g.plot",
                "l2g.settings", "l2g.mesh",
                "l2g.external"],
    data_files = [('', getDataFiles())],
    scripts = ['bin/flat', 'bin/submitFLAT', 'bin/FLAT.sbatch',
               'bin/imas2eqdsk', 'bin/mkEqdskMovie', 'bin/mkImasMovie',
               'bin/mkImasMovieFromPsi',
               'bin/med2mat', 'bin/mat2med', 'bin/plotMFlux', "bin/plotIP",
               "bin/get_disruption_profile_from_imas",
               "bin/torbeam_mapper", "bin/torbeam_plotter",
               "bin/get_disruption_profile_from_imas",
               "bin/mkDisruptionMovie"]
)
