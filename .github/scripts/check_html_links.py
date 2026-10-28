
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


ROOT = Path.cwd().resolve()


class LinkParser(HTMLParser):
    """Collect local file references from HTML markup."""

    def __init__(self) -> None:
        super().__init__()
        self.references: list[str] = []

    def handle_starttag(
        self, tag: str, attrs: list[tuple[str, str | None]]
    ) -> None:
        for attribute, value in attrs:
            if attribute in {"href", "src"} and value:
                self.references.append(value)


def check_html_file(html_file: Path) -> list[str]:
    parser = LinkParser()
    parser.feed(html_file.read_text(encoding="utf-8"))

    errors = []

    for reference in parser.references:
        url = urlsplit(reference)

        # Ignore external URLs, mailto:, data:, fragments, etc.
        if url.scheme or url.netloc or not url.path:
            continue

        path = Path(unquote(url.path))

        # Absolute URLs are relative to the published site root.
        target = ROOT / path.as_posix().lstrip("/") if path.is_absolute() else (
            html_file.parent / path
        )

        target = target.resolve()

        if not target.is_relative_to(ROOT):
            errors.append(f"{html_file}: path escapes repository: {reference}")
            continue

        if target.is_dir():
            target = target / "index.html"

        if not target.is_file():
            errors.append(f"{html_file}: missing target: {reference}")

    return errors


def main() -> int:
    html_files = sorted(
        path
        for path in ROOT.rglob("*.html")
        if ".git" not in path.parts
        and ".github" not in path.parts
        and "_site" not in path.parts
    )

    if not html_files:
        print("ERROR: no HTML files found")
        return 1

    errors = [
        error
        for html_file in html_files
        for error in check_html_file(html_file)
    ]

    if errors:
        print("\n".join(errors))
        return 1

    print(f"Validated local links in {len(html_files)} HTML file(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
