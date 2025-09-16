import nox

nox.options.sessions = ["lint", "typecheck"]


def editable_install(session: nox.Session):
    session.install("-r", "cmake/sweepgen-requirements.txt")
    session.install("-e", ".")


@nox.session(python="3.10", tags=["qa", "code-quality"])
def lint(session: nox.Session):
    """Lint code using flake8"""

    session.install("flake8")
    session.run("flake8", "src/sweepgen")


@nox.session(python="3.10", tags=["qa", "code-quality"])
def typecheck(session: nox.Session):
    """Run MyPy for static type checking"""
    editable_install(session)
    session.install("mypy")
    session.run("mypy", "src/sweepgen")


@nox.session
def user_manual(session: nox.Session):
    editable_install(session)

    session.chdir("user_manual")
    session.install("-r", "requirements.txt")

    if "--clean" in session.posargs:
        session.run("make", "clean", external=True)

    env = {}

    session_args = session.posargs
    if "--fail-on-warnings" in session_args:
        env["SPHINXOPTS"] = "-W --keep-going"

    session.run("make", "html", external=True, env=env)
