import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="stdpairs", # Replace with your own username
    version="1.0.6",
    author="Byeongsu Yu",
    author_email="byeongsu.yu@gmail.com",
    description="A library of SageMath doing symbolic computation over a monomial ideal of an affine (non-normal) semigroup ring",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/byeongsuyu/StdPairs",
    project_urls={
        "Bug Tracker": "https://github.com/byeongsuyu/StdPairs/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    python_requires=">=3.7",
)