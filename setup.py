import shutil
import os
import os.path
import subprocess
import re
from setuptools import setup
from setuptools.command.install_lib import install_lib
from setuptools.command.install import install
import setuptools.command.bdist_egg
import sys
#import glob




class install_lib_save_version(install_lib):
    """Save version information"""
    def run(self):
        install_lib.run(self)
        
        for package in self.distribution.command_obj["build_py"].packages:
            install_dir=os.path.join(*([self.install_dir] + package.split('.')))
            fh=open(os.path.join(install_dir,"version.txt"),"w")
            fh.write("%s\n" % (version))  # version global, as created below
            fh.close()
            pass
        pass
    pass



# Extract GIT version
if os.path.exists(".git"):
    # Check if tree has been modified
    modified = subprocess.call(["git","diff-index","--quiet","HEAD","--"]) != 0
    
    gitrev = subprocess.check_output(["git","rev-parse","HEAD"]).strip()

    version = "git-%s" % (gitrev)

    # See if we can get a more meaningful description from "git describe"
    try:
        versionraw=subprocess.check_output(["git","describe","--tags","--match=v*"],stderr=subprocess.STDOUT).decode('utf-8').strip()
        # versionraw is like v0.1.0-50-g434343
        # for compatibility with PEP 440, change it to
        # something like 0.1.0+50.g434343
        matchobj=re.match(r"""v([^.]+[.][^.]+[.][^-.]+)(-.*)?""",versionraw)
        version=matchobj.group(1)
        if matchobj.group(2) is not None:
            version += '+'+matchobj.group(2)[1:].replace("-",".")
            pass
        pass
    except subprocess.CalledProcessError:
        # Ignore error, falling back to above version string
        pass

    if modified and version.find('+') >= 0:
        version += ".modified"
        pass
    elif modified:
        version += "+modified"
        pass
    pass
else:
    version = "UNKNOWN"
    pass

print("version = %s" % (version))

tortuosity_tracing_package_files=[ "pt_steps/*" ]


console_scripts=["tiffs_to_xlg","svg_histograms",'point_spacing_calc']
#gui_scripts = []  # Could move graphical scrips into here to eliminate stdio window on Windows (where would error messages go?)

console_scripts_entrypoints = [ "%s = tortuosity_tracing.bin.%s:main" % (script,script.replace("-","_")) for script in console_scripts ]

#gui_scripts_entrypoints = [ "%s = limatix.bin.%s:main" % (script,script.replace("-","_")) for script in gui_scripts ]


setup(name="tortuosity_tracing",
      description="Tracing of crack tortuosity",
      author="Chevonne McInnis and Stephen D. Holland",
      version=version,
      # url="http://",
      zip_safe=False,
      packages=["tortuosity_tracing",
                "tortuosity_tracing.bin"],
      #data_files=[ ("share/tortuosity_tracing/pt_steps",pt_steps_files),]
      cmdclass={"install_lib": install_lib_save_version },
      package_data={"tortuosity_tracing": tortuosity_tracing_package_files},
      entry_points={"limatix.processtrak.step_url_search_path": [ "limatix.share.pt_steps = tortuosity_tracing:getstepurlpath" ],
                    "console_scripts": console_scripts_entrypoints,
                    #"gui_scripts": gui_scripts_entrypoints 
                })


