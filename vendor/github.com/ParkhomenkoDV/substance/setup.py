import time

from setuptools import find_packages, setup

with open("README.md", "rt", encoding="utf-8") as file:
    long_description = file.read()

with open("requirements.txt", "rt") as file:
    install_requires = file.readlines()

setup(
    name="substance",
    version=time.strftime("%Y.%m.%d.%H", time.localtime()),
    description="lib",
    long_description=long_description,
    long_description_content_type="text/markdown",  # если long_description = .md
    author="Daniil Andryushin",
    author_email="",
    url="https://github.com/ParkhomenkoDV/substance.git",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=install_requires,
    package_data={
        "substance": [],  # доп. файлы библиотеки
    },
)
