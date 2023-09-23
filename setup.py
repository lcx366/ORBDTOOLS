from setuptools import setup,find_packages 

setup(
    name='orbdtools',
    version='0.1.1',
    description='A set of routines for data processing related to ORBit Determination(ORBD) of space objects',
    author='Chunxiao Li',
    author_email='lcx366@126.com',
    url='https://github.com/lcx366/ORBDTOOLS',
    license='MIT',
    long_description_content_type='text/markdown',
    long_description=open('README.md', 'rb').read().decode('utf-8'),
    keywords = ['Arc Matching','Arc Association','Initial OD','Cataloging OD','Precise OD'],
    python_requires = '>=3.8',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'License :: OSI Approved :: MIT License',
        ],
    packages = find_packages(),
    include_package_data=True,
    install_requires=[
        'astropy>=4.3.1',
        'skyfield',
        'sgp4',
        'spacetrack',
        'numpy',
        'scipy',
        'loess',
        'wget',
        'pandas',
        'colorama'
        ],
)
