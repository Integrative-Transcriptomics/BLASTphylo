from setuptools import setup

setup(
    name='server',
    packages=['server'],
    include_package_data=True,
    install_requires=[
        'flask',
    ],
    package_data={'server': ['data/ncbi_normalisation.json', 'data/searchbar_entries.json',
                                 'test_example/comparison_example/*', 'test_example/mapping_example/*']},
)