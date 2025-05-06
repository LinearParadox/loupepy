import platform
from os import PathLike
from urllib.request import urlretrieve
import hashlib
from pathlib import Path
import os
OPERATING_SYSTEM = platform.system()

def _get_checksum() -> tuple[str, str]:
    '''
    Returns the checksum of the loupe converter binary
    '''
    if OPERATING_SYSTEM.startswith('win'):
        return ("https://github.com/10XGenomics/loupeR/releases/"
        "download/v1.1.4/louper-windows-x64.exe","f5d1e99138e840169a19191d10bb25ab")
    elif OPERATING_SYSTEM.startswith('darwin') or OPERATING_SYSTEM.startswith('apple'):
        return ("https://github.com/10XGenomics/loupeR/"
        "releases/download/v1.1.4/louper-macos-x64", "ea65a2ec372d623c54d45c51793014e2")
    else:
        return ("https://github.com/10XGenomics/loupeR/"
                       "releases/download/v1.1.4/louper-linux-x64","b3fd93fd88a43fbcf3f6e40af3186eaa")


def _md5_checksum() -> str:
    '''
    Downloads a file from a URL and returns the MD5 checksum of the file. Returns the link to download
    '''
    link, checksum = _get_checksum()
    with open(urlretrieve(link)[0], 'rb') as f:
        md5_hash = hashlib.md5()
        for chunk in iter(lambda: f.read(4096), b""):
            md5_hash.update(chunk)
    computed_md5 = md5_hash.hexdigest()
    if computed_md5 != checksum:
        raise OSError(f"Checksum mismatch: {computed_md5} != {checksum}. If this continues please raise an issue on github")
    return link

def _get_install_path() -> Path:
    '''
    A function to return the path to the users config directory.
    '''
    home = Path(os.path.expanduser("~"))
    if OPERATING_SYSTEM.startswith("win"):
        base = Path(os.environ.get("APPDATA", home / "AppData" / "Roaming" / "loupepy"))
    elif OPERATING_SYSTEM == "darwin":
        base = home / "Library" / "Application Support" / "loupepy"
    else:
        base = Path(os.environ.get("XDG_DATA_HOME", home / ".local" / "share" / "loupepy"))
    return base

def _download_loupe_converter(path: Path) -> None:
    '''
    Downloads the loupeR binary to the specified path
    '''
    if platform.architecture()[0] != '64bit':
        raise OSError('Only 64-bit operating systems are supported by loupe converter')
    path = path / 'loupe_converter'
    link = _md5_checksum()
    urlretrieve(link, path)
    path.chmod(0o755)
    print(path)

def setup(path: Path | None | str = None) -> None:
    '''
    Downloads the loupe converter binary to the specified path.
    If no path is specified, it will be downloaded to the default location for the operating system.
    '''
    if path is None:
        path = _get_install_path()
    if path is str:
        path = Path(path)
    eula()
    _download_loupe_converter(path) # type: ignore





def eula() ->  None:
    '''
    Prompts the user to agree to the EULA, as it is in the R version of the tool
    '''
    print("This tool is independently maintained,"
        " but utilizes 10x genomics tools to perform conversions")
    print("By using this tool, you agree to the 10X Genomics End User License Agreement "
        "(https://www.10xgenomics.com/legal/end-user-software-license-agreement).")
    print("Do you agree to the terms of the EULA? (yes/no)")
    response = input()
    if response.lower() != 'yes':
        raise OSError('You must agree to the EULA to use this tool')
    path = _get_install_path()
    path.mkdir(parents=False, exist_ok=True)
    path = path / 'eula'
    path.touch(exist_ok=True)






def eula_reset() -> None:
    '''
    Resets the EULA agreement
    '''
    path = _get_install_path()
    path = path / 'loupepy' / 'eula'
    path.unlink()
    path=path.parent / 'loupe_converter'
    path.unlink(missing_ok=True)
    path.parent.rmdir()

