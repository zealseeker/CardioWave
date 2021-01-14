import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="CardioWave",
    version=__import__("cdwave").__version__,
    author="Hongbin Yang",
    author_email="yanyanghong@163.com",
    description="A tool to derive parameters from waveform data for cardiotoxicity research",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/zealseeker/cardiowave",
    packages=setuptools.find_packages(),
    entry_points={'console_scripts': ['cardiowave=cdwave.viewer.gui:main']},
    include_package_data=True,
    package_data = {
        'cdwave': ['viewer/*', 'param_annot.json']
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    zip_safe=False,
)
