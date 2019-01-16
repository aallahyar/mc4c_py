import getpass
from datetime import datetime
import socket


def getRuntime():
    runtime = dict()
    runtime['version']=getVersion()
    runtime['datetime']=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    runtime['hostname']=socket.gethostname()
    runtime['username']=getpass.getuser()
    return runtime


def getVersion():
    version = 'unknown'
    try:
        import subprocess
        version = subprocess.check_output(["git", "describe", "--always"]).split()[0]
    except:
        pass
    return version


def printArgs(args):
    argdict = vars(args)
    print 'tool =', str(argdict['func']).split()[1]
    for arg in sorted(argdict.keys()):
        if arg != 'func':
            print arg, '=', argdict[arg]
