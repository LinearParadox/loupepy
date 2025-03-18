import platform
from urllib.request import urlretrieve
import hashlib
import os
import typing
OPERATING_SYSTEM = platform.system()


def _md5_checksum(url: str, path: str) -> str:
    '''
    Downloads a file from a URL and returns the MD5 checksum of the file
    '''
    urlretrieve(url, path)
    with open(path, 'rb') as f:
        hash = hashlib.md5()
        for chunk in iter(lambda: f.read(4096), b""):
            hash.update(chunk)
    return hash.hexdigest()


def download_loupeR(path: str) -> None: 
    '''
    Downloads the loupeR binary to the specified path
    '''
    if platform.architecture()[0] != '64bit':
        raise OSError('Only 64-bit operating systems are supported')
    if OPERATING_SYSTEM.startswith('linux'):
        tmp_path = path + '.tmp'
        md5=_md5_checksum("https://github.com/10XGenomics/loupeR/releases/download/v1.1.4/louper-linux-x64", tmp_path)
        if md5 != "b3fd93fd88a43fbcf3f6e40af3186eaa":
            raise OSError('MD5 checksum failed, please try again. If this problem persists, please raise an issue on the github.')
        os.rename(tmp_path, "")
        os.chmod(path, 0o755)
    elif OPERATING_SYSTEM.startswith('win'):
        tmp_path = os.environ['TEMP'] + '/louper.exe'
        md5 = _md5_checksum("https://github.com/10XGenomics/loupeR/releases/download/v1.1.4/louper-windows-x64.exe", tmp_path)
        if md5 != "f5d1e99138e840169a19191d10bb25ab":
            raise OSError('MD5 checksum failed, please try again. If this problem persists, please raise an issue on the github.')
        os.rename(tmp_path, path)
    elif OPERATING_SYSTEM.startswith('darwin'):
        tmp_path = os.environ['TMPDIR'] + '/louper.app'
        md5 = _md5_checksum("https://github.com/10XGenomics/loupeR/releases/download/v1.1.4/louper-macos-x64", tmp_path)
        if md5 != "f5d1e99138e840169a19191d10bb25ab":
            raise OSError('MD5 checksum failed, please try again. If this problem persists, please raise an issue on the github.')
        os.rename(tmp_path, path)
    else:
        raise OSError('Operating system not supported')
    print('Successfully downloaded loupeR to ' + path)
    
def setup(path=None):
    if path is None:
        if OPERATING_SYSTEM.startswith('linux'):
            path = os.environ['HOME'] + '/.local/bin/loupe/louper'
        elif OPERATING_SYSTEM.startswith('win'):
            path = os.environ['LOCALAPPDATA'] + '/Programs/louper/louper.exe'
        elif OPERATING_SYSTEM.startswith('darwin'):
            path = os.environ['HOME'] + '/Applications/loupe/louper.app'
    download_loupeR(path)
    eula(path)

def eula(path: str) -> None:
    '''
    Prompts the user to agree to the EULA, as it is in the R version of the tool
    '''
    if not os.path.exists(f'{path}/eula/eula_agreement'):
        print("This tool is independently maintained, but utilizes 10x genomics tools to perform conversions")
        print("By using this tool, you agree to the 10X Genomics End User License Agreement (https://www.10xgenomics.com/legal/terms-of-use/).")
        print("Do you agree to the terms of the EULA? (yes/no)")
        response = input()
        if response.lower() != 'yes':
            raise OSError('You must agree to the EULA to use this tool')
        os.makedirs(f'{path}/eula')
        with open(f'{path}/eula/eula_agreement', 'w') as f:
            f.write('agreed')

def eula_reset(path: str) -> None:
    '''
    Resets the EULA agreement
    '''
    if os.path.exists(f'{path}/eula/eula_agreement'):
        os.remove(f'{path}/euale/eula_agreement')
        print('EULA agreement reset')                        