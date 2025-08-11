from loupepy.setup import eula_reset, setup, _get_checksum, _md5_checksum  # type: ignore
import os
import platform
import pytest

def test_eula_and_reset(monkeypatch, tmp_path_factory):
    output_dir = tmp_path_factory.mktemp("fake_directory")
    monkeypatch.setattr('builtins.input', lambda _: "y")
    setup(output_dir)
    assert os.path.exists(output_dir / "eula")
    assert os.path.exists(output_dir / "loupe_converter")
    eula_reset(output_dir)
    assert not os.path.exists(output_dir / "eula")
    assert not os.path.exists(output_dir / "loupe_converter")

@pytest.mark.parametrize("fake_system", ["Windows", "Linux", "Darwin"])
def test_checksum(monkeypatch, fake_system):
    """
    Test the _get_checksum function for different operating systems.
    """
    monkeypatch.setattr(platform, "system", lambda: fake_system)
    for n in range(0,3):
        try:
            link = _md5_checksum()
            break
        except OSError:
            continue
    assert link == _get_checksum()[0]