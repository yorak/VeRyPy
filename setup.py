from setuptools import setup

classifiers = [
    "Development Status :: Development",
    "License :: MIT License",
    "Natural Language :: English",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.8",
]


with open("README.md") as f:
    long_description = f.read()

with open("requirements.txt") as f:
    requirements = f.readlines()

install_requires = [r.strip() for r in requirements]

setup(
    name="verypy",
    version="0.0.1",
    description="A python library with implementations of classical heuristics for the capacitated vehicle routing problem",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yorak/VeRyPy",
    license="MIT License",
    platforms="any",
    packages=["verypy"],
    include_package_data=True,
    install_requires=install_requires,
)