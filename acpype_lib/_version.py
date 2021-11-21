from subprocess import run, STDOUT, PIPE

out = run("git describe --tags --always", shell=True, stderr=STDOUT, stdout=PIPE)

if out.returncode == 0:

    version = out.stdout.decode()[0:10]
else:
    try:
        from importlib.metadata import version as ver

        version = ver("acpype")
    except Exception:
        import pkg_resources

        version = str(pkg_resources.get_distribution("acpype").version)
