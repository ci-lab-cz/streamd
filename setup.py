import setuptools
from os import path
import streamd

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name="streamd",
    version=streamd.__version__,
    author="Aleksandra Ivanova, Olena Mokshyna, Pavel Polishchuk",
    description="Streamd Python module to facilitate molecular dynamics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ci-lab-cz/streamd",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: POSIX :: Linux ",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    python_requires='>=3.6',
    extras_require={
        'rdkit': ['rdkit>=2017.09'],
    },
    entry_points={'console_scripts':
                      ['run_md = streamd.run_md:main',
                       'run_gbsa = streamd.run_gbsa:main',
                       'run_prolif = streamd.prolif.run_prolif:main',
                       'prolif_drawmap = streamd.prolif.prolif2png:main',
                       'prolif_draw_by_frame = streamd.prolif.prolif_frame_map:main']},
    include_package_data=True
)
