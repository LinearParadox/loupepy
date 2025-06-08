from loupepy.setup import eula_reset, setup  # type: ignore
import os

def test_eula_and_reset(monkeypatch, tmp_path_factory):
    output_dir = tmp_path_factory.mktemp("fake_directory")
    monkeypatch.setattr('builtins.input', lambda _: "y")
    setup(output_dir)
    assert os.path.exists(output_dir / "eula")
    assert os.path.exists(output_dir / "loupe_converter")
    eula_reset(output_dir)
    assert not os.path.exists(output_dir / "eula")
    assert not os.path.exists(output_dir / "loupe_converter")