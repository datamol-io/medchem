import typer

app = typer.Typer(help="The Medchem CLI", add_completion=False)


@app.command(help="A dummy CLI #1")
def dummy1():
    print("dummy1")


@app.command(help="A dummy CLI #2")
def dummy2():
    print("dummy2")


if __name__ == "__main__":
    app()  # pragma: no cover
