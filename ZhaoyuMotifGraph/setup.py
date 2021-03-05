import setuptools

with open("README.md", "r") as fh:
  long_description = fh.read()

setuptools.setup(
  name="ZhaoyuMotifGraph",
  version="0.0.7",
  author="Matt Bright,Zhaoyu Han",
  author_email="skyfallhalo@hotmail.com",
  description="package for motif graph",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/MattB-242/Motif_Complex",
  packages=setuptools.find_packages(),
  classifiers=[
  "Programming Language :: Python :: 3",
  "License :: OSI Approved :: MIT License",
  "Operating System :: OS Independent",
  ],
)
