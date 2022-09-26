from setuptools import setup, find_packages



with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
setup(
    name='auto-abdab',
    version='0.0.1',
    description='auto-abdab: Pipeline to automatically generate disease-specific antibody databases',
    license='BSD 3-clause license',
    maintainer='Emily Jin, Fabian Spoendlin, Gemma Gordon, Hongyu Qian, Jesse Murray',
    long_description=long_description,
    long_description_content_type='text/markdown',
    include_package_data=True,
    packages=find_packages(include=('auto_abdab', 'auto_abdab.*')),
    package_data={
        '': ['*.txt']
    },
    install_requires=[
        'pandas'
    ],
)