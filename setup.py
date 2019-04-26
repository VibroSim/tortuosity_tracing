import shutil
import os.path
from setuptools import setup
from setuptools.command.install_lib import install_lib
from setuptools.command.install import install
import setuptools.command.bdist_egg
import sys
#import glob

tortuosity_tracing_package_files=[ "pt_steps/*" ]


console_scripts=["tiffs_to_xlg","svg_histograms"]
#gui_scripts = []  # Could move graphical scrips into here to eliminate stdio window on Windows (where would error messages go?)

console_scripts_entrypoints = [ "%s = tortuosity_tracing.bin.%s:main" % (script,script.replace("-","_")) for script in console_scripts ]

#gui_scripts_entrypoints = [ "%s = limatix.bin.%s:main" % (script,script.replace("-","_")) for script in gui_scripts ]


setup(name="tortuosity_tracing",
      description="Tracing of crack tortuosity",
      author="Chevonne McInnis and Stephen D. Holland",
      # url="http://",
      zip_safe=False,
      packages=["tortuosity_tracing",
                "tortuosity_tracing.bin"],
      #data_files=[ ("share/tortuosity_tracing/pt_steps",pt_steps_files),]
      package_data={"tortuosity_tracing": tortuosity_tracing_package_files},
      entry_points={"limatix.processtrak.step_url_search_path": [ "limatix.share.pt_steps = tortuosity_tracing:getstepurlpath" ],
                    "console_scripts": console_scripts_entrypoints,
                    #"gui_scripts": gui_scripts_entrypoints 
                })


