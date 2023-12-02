from setuptools import setup, find_packages
from io import open # Should be Py2.7 compatible

classifiers = [
    "Development Status :: Development",
    "License :: MIT License",
    "Natural Language :: English",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
]


with open("README.md", encoding="utf-8") as f:
    long_description = f.read()

with open("requirements.txt", encoding="utf-8") as f:
    requirements = f.readlines()

install_requires = [r.strip() for r in requirements]

with open("verypy/__init__.py") as f:
    for line in f:
        if line.startswith("__version__"):
            module_version = line.split('"')[1]

setup(
    name="verypy",
    version=module_version,
    description="A python library with implementations of classical heuristics for the capacitated vehicle routing problem",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yorak/VeRyPy",
    license="MIT License",
    platforms="any",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'VeRyPy = verypy.VeRyPy:main',
        ],
    },
    include_package_data=True,
    install_requires=install_requires,
)
