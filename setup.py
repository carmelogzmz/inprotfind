from setuptools import setup, find_packages

setup(
    name='inprotfind',  # Cambia esto al nombre de tu librería
    version='0.1',
    description='Librería para análisis de secuencias proteicas utilizando MMseqs2',
    author='Carmelo Gómez-Martínez',
    author_email='carmelogzmz@gmail.com',
    url='https://github.com/carmelogzmz/inprotfind',  # Cambia esto por la URL de tu repositorio
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'inprotfind': [
            'query_examples/*',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',  # Ajusta según la versión mínima de Python requerida
    install_requires=[
        'pandas',
        'matplotlib',
        'colorama',
        'biopython',
        'pyarrow',
        'fastparquet',
        'requests',
        'tqdm',
        # Añade aquí cualquier otra dependencia externa que tu código requiera
    ],
    entry_points={
    'console_scripts': [
        'inprotfind=inprotfind.inprotfind:main_function',  # Cambia esto según tu configuración
    ],
},
)