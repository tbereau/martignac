from pathlib import Path

from setuptools import setup


def _get_version() -> str:
    """Read martignac/VERSION.txt and return its contents."""
    path = Path("martignac").resolve()
    version_file = path / "VERSION.txt"
    return version_file.read_text().strip()


version = _get_version()


with open("README.md", encoding="utf-8") as readme_file:
    readme = readme_file.read()


requirements = [
    "signac>=2.0.0",
    "signac-flow>=0.25.1",
    "numpy>=1.25.0,<2.0.0",
    "regex",
    "pymbar>=4.0.1",
    "jupyter>=1.0.0",
    "jupyterlab",
    "pandas>=2.0.2",
    "h5py",
    "signac-dashboard",
    "alchemlyb",
    "networkx",
    "confuse",
    "mdanalysis>=2.6.1",
    "python-decouple",
    "marshmallow-dataclass>=8.6.0",
    "cachetools>=5.3.2",
    "isort>=5.8.0",
    "insane @ git+https://github.com/Tsjerk/Insane.git@f382b447d6dd6c3a89f884347cceaede02bf0166",
    "anyio>=3.0.0,<4.0",
    "GromacsWrapper",
    "mkdocs",
    "mkdocs-material",
    "mkdocstrings[python]",
]

setup(
    name="martignac",
    version=version,
    description="Computational workflow for Martini coarse-grained simulations",
    long_description=readme,
    long_description_content_type="text/markdown",
    author="Tristan Bereau",
    author_email="bereau@uni-heidelberg.de",
    url="https://lin0.thphys.uni-heidelberg.de:4443/bereau/martignac",
    project_urls={
        "Documentation": "https://lin0.thphys.uni-heidelberg.de:4443/bereau/martignac",
        "Issues": "https://lin0.thphys.uni-heidelberg.de:4443/bereau/martignac/-/issues",
    },
    packages=["martignac"],
    package_dir={"martignac": "martignac"},
    entry_points={"console_scripts": ["martignac = martignac.__main__:main"]},
    include_package_data=True,
    python_requires=">=3.9",
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Programming Language :: Python",
        "Topic :: Software Development",
    ],
    keywords=[
        "martignac",
        "Python",
        "Martini",
        "Coarse-grained simulations",
        "CG simulations",
    ],
)
