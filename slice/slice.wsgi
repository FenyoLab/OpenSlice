import sys
sys.path.append(".")

import os
# Change working directory so relative paths (and template lookup) work again
os.chdir(os.path.dirname(__file__))

import bottle
import slice 
#bottle.debug(True)
application = bottle.default_app()
