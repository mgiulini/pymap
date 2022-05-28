"""Setup.py file."""
from setuptools import setup

setup(
    name='pymap',
    version='1.0.0',
    description='pymap',
    long_description="",
    long_description_content_type='text/markdown',
    license='Apache License 2.0',
    author='Marco Giulini',
    author_email='mrcgiulini@gmail.com',
    url='https://github.com/mgiulini/pymap',
    include_package_data=True,
    zip_safe=False,
    # Some weblinks
    project_urls={
        'webpage': 'https://github.com/mgiulini/pymap',
        'Documentation': 'https://github.com/mgiulini/pymap/README',
        'Changelog': '',
        'Issue Tracker': 'https://github.com/mgiulini/pymap/issues',
        'Discussion Forum': 'https://github.com/mgiulini/pymap/issues',
        },
    keywords=[
        'Information Theory',
        'Entropy',
        ],
    python_requires='>=3.9, <3.11',
    )
