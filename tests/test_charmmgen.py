"""Tests for locating and fetching charmmgen, which CHARMM output depends on.

antechamber writes CHARMM files by invoking ``$AMBERHOME/bin/charmmgen`` **by absolute
path**; a copy merely on PATH is never used. Modern AmberTools dropped charmmgen, so an
installation without ACPYPE's bundled AmberTools has none, and the run used to announce
"Writing CHARMM files" and then write nothing at all, silently.
"""

import hashlib
import io
import os
import tarfile
from pathlib import Path

import pytest

from acpype import utils
from acpype.errors import AcpypeError
from acpype.utils import CHARMMGEN_SHA256, amber_bin_dir, charmmgen_path, charmmgen_platform, download_charmmgen

REPO_TARBALLS = {key: Path("charmmgen_%s.tgz" % key) for key in ("linux", "macos")}


def make_tarball(payload: bytes = b"#!/bin/sh\necho charmmgen\n") -> bytes:
    """Build a gzipped tar holding one charmmgen at the path the real ones use."""
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        info = tarfile.TarInfo("src/acpype/amber_linux/bin/charmmgen")
        info.size = len(payload)
        info.mode = 0o755
        archive.addfile(info, io.BytesIO(payload))
    return buffer.getvalue()


@pytest.fixture
def served(monkeypatch: pytest.MonkeyPatch):
    """Serve a crafted archive from a fake urlopen, and record the URL requested."""

    def serve(blob: bytes) -> list[str]:
        requested: list[str] = []

        def fake_urlopen(url, timeout=None):
            # BytesIO is already a context manager, which is all urlopen is used as.
            requested.append(url)
            return io.BytesIO(blob)

        monkeypatch.setattr(utils, "urlopen", fake_urlopen)
        monkeypatch.setitem(CHARMMGEN_SHA256, "linux", hashlib.sha256(blob).hexdigest())
        monkeypatch.setattr(utils, "charmmgen_platform", lambda: "linux")
        return requested

    return serve


@pytest.mark.parametrize(
    ("system", "machine", "expected"),
    [
        ("darwin", "arm64", "macos"),
        ("darwin", "x86_64", "macos"),
        ("linux", "x86_64", "linux"),
        ("linux", "AMD64", "linux"),
        ("linux", "aarch64", None),
        ("win32", "AMD64", None),
    ],
)
def test_charmmgen_platform(monkeypatch: pytest.MonkeyPatch, system: str, machine: str, expected: str | None) -> None:
    """Only macOS, whose build is universal, and Linux x86_64 have a build."""
    monkeypatch.setattr(utils.sys, "platform", system)
    monkeypatch.setattr(utils.platform, "machine", lambda: machine)

    assert charmmgen_platform() == expected


def test_amber_bin_dir_follows_the_antechamber_on_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """An antechamber on PATH wins over the bundle, because that is the one that runs.

    set_for_pip wires the bundle in only when nothing is on PATH, so preferring the
    bundle here would check one AmberTools while antechamber resolved charmmgen
    against another.
    """
    monkeypatch.setattr(utils, "bundled_amber_dir", lambda: str(tmp_path / "bundle"))
    monkeypatch.setattr(utils, "which", lambda name: str(tmp_path / "conda" / "bin" / name))

    assert amber_bin_dir() == str(tmp_path / "conda" / "bin")


def test_amber_bin_dir_falls_back_to_the_bundle(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """With nothing on PATH it is the bundled AmberTools, as a wheel install has."""
    monkeypatch.setattr(utils, "which", lambda name: None)
    monkeypatch.setattr(utils, "bundled_amber_dir", lambda: str(tmp_path / "bundle"))

    assert amber_bin_dir() == str(tmp_path / "bundle" / "bin")


def test_amber_bin_dir_without_any_ambertools(monkeypatch: pytest.MonkeyPatch) -> None:
    """Neither on PATH nor bundled is reported as nothing found."""
    monkeypatch.setattr(utils, "which", lambda name: None)
    monkeypatch.setattr(utils, "bundled_amber_dir", lambda: None)

    assert amber_bin_dir() is None


def test_charmmgen_path_looks_where_antechamber_looks(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """It reports the binary beside antechamber, and ignores one merely on PATH."""
    binDir = tmp_path / "bin"
    binDir.mkdir()
    monkeypatch.setattr(utils, "amber_bin_dir", lambda: str(binDir))
    assert charmmgen_path() is None

    binary = binDir / "charmmgen"
    binary.write_text("#!/bin/sh\n")
    binary.chmod(0o755)
    assert charmmgen_path() == str(binary)


def test_charmmgen_path_without_any_ambertools(monkeypatch: pytest.MonkeyPatch) -> None:
    """No AmberTools at all is reported as no charmmgen, not an error."""
    monkeypatch.setattr(utils, "amber_bin_dir", lambda: None)

    assert charmmgen_path() is None


def test_download_installs_beside_antechamber(served, tmp_path: Path) -> None:
    """The binary lands in the AmberTools bin, executable, with the pinned URL fetched."""
    blob = make_tarball()
    requested = served(blob)

    path = download_charmmgen(destDir=str(tmp_path))

    assert path == str(tmp_path / "charmmgen")
    assert Path(path).read_bytes() == b"#!/bin/sh\necho charmmgen\n"
    assert os.access(path, os.X_OK)
    assert requested == [utils.CHARMMGEN_URL.format(ref=utils.CHARMMGEN_REF, key="linux")]


def test_download_refuses_a_tampered_archive(served, tmp_path: Path) -> None:
    """A digest that does not match the recorded one is refused, and nothing is written."""
    served(make_tarball())
    CHARMMGEN_SHA256["linux"] = "0" * 64

    with pytest.raises(AcpypeError, match="does not match its recorded checksum"):
        download_charmmgen(destDir=str(tmp_path))

    assert not (tmp_path / "charmmgen").exists()


def test_download_refuses_an_archive_without_the_binary(served, tmp_path: Path) -> None:
    """An archive holding no charmmgen is an error rather than an empty install."""
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        info = tarfile.TarInfo("README")
        info.size = 2
        archive.addfile(info, io.BytesIO(b"hi"))
    served(buffer.getvalue())

    with pytest.raises(AcpypeError, match="holds no charmmgen"):
        download_charmmgen(destDir=str(tmp_path))


def test_download_on_an_unsupported_platform(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """A platform with no build says so, and points at the Docker image."""
    monkeypatch.setattr(utils, "charmmgen_platform", lambda: None)

    with pytest.raises(AcpypeError, match="no charmmgen build for"):
        download_charmmgen(destDir=str(tmp_path))


def test_download_when_the_destination_is_read_only(served, tmp_path: Path) -> None:
    """An unwritable AmberTools yields a copy-and-paste command instead of a failure."""
    served(make_tarball())
    readOnly = tmp_path / "bin"
    readOnly.mkdir()
    readOnly.chmod(0o500)
    try:
        with pytest.raises(AcpypeError, match=r"(?s)is not writable.*curl -sL"):
            download_charmmgen(destDir=str(readOnly))
    finally:
        readOnly.chmod(0o700)


def test_download_needs_an_ambertools(monkeypatch: pytest.MonkeyPatch, served) -> None:
    """With no AmberTools anywhere there is nowhere to install to, and it says so."""
    served(make_tarball())
    monkeypatch.setattr(utils, "amber_bin_dir", lambda: None)

    with pytest.raises(AcpypeError, match="no AmberTools found"):
        download_charmmgen()


@pytest.mark.parametrize("key", ["linux", "macos"])
def test_recorded_digests_match_the_repository_tarballs(key: str) -> None:
    """The digests baked into the source are the ones of the tarballs served from the tag."""
    tarball = REPO_TARBALLS[key]
    if not tarball.is_file():
        pytest.skip(f"{tarball} is not in this installation")

    assert hashlib.sha256(tarball.read_bytes()).hexdigest() == CHARMMGEN_SHA256[key]


@pytest.mark.parametrize("key", ["linux", "macos"])
def test_repository_tarballs_hold_exactly_one_charmmgen(key: str) -> None:
    """Each tarball offers the downloader one unambiguous binary to pick."""
    tarball = REPO_TARBALLS[key]
    if not tarball.is_file():
        pytest.skip(f"{tarball} is not in this installation")

    with tarfile.open(tarball) as archive:
        names = archive.getnames()

    assert [n for n in names if n.endswith("bin/charmmgen")] == [f"src/acpype/amber_{key}/bin/charmmgen"]


def test_download_ignores_the_appledouble_sibling(served, tmp_path: Path) -> None:
    """Both tarballs were made on macOS and carry a ``._charmmgen``; it must not be installed."""
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        for name, payload in (("._charmmgen", b"apple double junk"), ("charmmgen", b"#!/bin/sh\nreal\n")):
            info = tarfile.TarInfo(f"src/acpype/amber_linux/bin/{name}")
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))
    served(buffer.getvalue())

    path = download_charmmgen(destDir=str(tmp_path))

    assert Path(path).read_bytes() == b"#!/bin/sh\nreal\n"
