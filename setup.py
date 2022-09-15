#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from setuptools import setup, find_packages
setup(
  author="Charles Turner",
  author_email='charlesturner0987@gmail.com',
  classifiers=[
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3.9',
  ],
  install_requires=["numpy","pandas","shutil","os","warmings","scipy","lmfit"
  ,"matplotlib","openpyxl"],
  description="This is the cleaned up organic-alkalinity-sausage-machine, now known more sensibly as OrgAlkCalc",
  py_modules=["SausageMachine"],
  package_dir={'': 'src'},
  license="MIT license",
  include_package_data=True,
  name='OrgAlkCalc',
  version='0.1.0',
  zip_safe=False,
)