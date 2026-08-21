from setuptools import find_packages, setup

setup(
    name='mmcFoamReader',
    packages=find_packages(include=['mmcFoamReader']),
    version='1.0.0',
    description='Library to read mmcFoam data into python for post processing',
    author='Jan Wilhelm Gaertner',
    license='GNUv3',
	install_requires=['numpy==1.25.0','scipy'],
	tests_require=['pytest-runner','pytest==6.2.5'],
	test_suite='tests',
)
