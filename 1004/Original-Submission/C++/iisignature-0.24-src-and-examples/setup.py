from setuptools import setup, Extension
import numpy, os, sys

version = "0.24"
headers = ["bch","calcSignature","logsig","logSigLength","makeCompiledFunction",
            "rotationalInvariants","readBCHCoeffs"]

args = []
link_args=[]
if os.name == 'posix':
    args=['-std=c++11']

if sys.platform == 'darwin' and "MACOSX_DEPLOYMENT_TARGET" not in os.environ:
    if "IISIGNATURE_MACOSX_DONOTBECLEVER" not in os.environ:
        args.append('-mmacosx-version-min=10.9')
        link_args.append('-mmacosx-version-min=10.9')
        #os.environ["MACOSX_DEPLOYMENT_TARGET"]="10.9"

xtn = Extension('iisignature', ['src/pythonsigs.cpp'], 
                extra_compile_args=args,
                extra_link_args=link_args,
                define_macros=[("VERSION",version)],
                include_dirs=[numpy.get_include()],
                depends=["src/"+i+".hpp" for i in headers])

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='iisignature',
      version=version,
      ext_modules=[xtn],
      description='Iterated integral signature calculations',
      long_description=readme(),
      classifiers=[
          # How mature is this project? Common values are
          #   3 - Alpha
          #   4 - Beta
          #   5 - Production/Stable
          'Development Status :: 4 - Beta',

          # Indicate who your project is intended for
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Intended Audience :: Financial and Insurance Industry',
          'Topic :: Scientific/Engineering :: Artificial Intelligence',
          'Topic :: Scientific/Engineering :: Mathematics',
          
          # Pick your license as you wish (should match "license" above)
          'License :: OSI Approved :: MIT License',
          
          # Specify the Python versions you support here. In particular, ensure
          # that you indicate whether you support Python 2, Python 3 or both.
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 3',
          ],
      url='https://github.com/bottler/iisignature',
      author='Jeremy Reizenstein',
      author_email='j.f.reizenstein@warwick.ac.uk',
      keywords = ["signature", "logsignature"],
      license='MIT',
      packages=["iisignature_data"],
      test_suite="tests.my_module_suite",
      #data_files=[("iisignature_data",["src/bchLyndon20.dat"])],
      package_data={"iisignature_data": ["bchLyndon20.dat"]},
      zip_safe=False,           
      install_requires=['numpy>1.7'], #For now, this is redundant as we've died above here on an import statement
      setup_requires=['numpy>1.7'],
)
#should add source url from github like at https://www.python.org/dev/peps/pep-0426/
