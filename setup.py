try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'name': 'bci',
    'version': '0.1',
    'packages': ['bci'],
    'install_requires': ['nose', 'numpy', 'MDSplus'],
        # 'filters'],  # <- Not on PyPI
    'author': 'Evan M. Davis',
    'author_email': 'emd@mit.edu',
    'url': '',
    'description': 'Python tools for DIII-D bi-color interferometer signals.'
}

setup(**config)
