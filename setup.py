import setuptools

with open("readme.md", "r", encoding="utf-8") as ff:
    long_description = ff.read()

setuptools.setup(
    name='qlcp',
    version='0.24.913',
    author='Dr Jie Zheng & Dr Lin-Qiao Jiang',
    author_email='jiezheng@nao.cas.cn',
    description='Quick Light Curve Pipeline v2024', # short description
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/drjiezheng/qlcp',
    # packages=setuptools.find_packages(),
    packages=['qlcp', ],
    include_package_data=True,
    package_data={"qlcp": ["default.*", "bright.*",]},
    python_requires='>=3.10',
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        'Operating System :: POSIX',
        'Operating System :: MacOS',
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    install_requires=['numpy', 'scipy', 'matplotlib', 'astropy', 
        'tqdm', 'qastutil', 'qmatch']
)
