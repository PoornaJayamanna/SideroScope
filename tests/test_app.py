"""Smoke tests: run the Streamlit app headlessly and check it does not raise."""

from pathlib import Path

from streamlit.testing.v1 import AppTest

APP = str(Path(__file__).resolve().parent.parent / "app.py")


def new_app(mode=None) -> AppTest:
    at = AppTest.from_file(APP, default_timeout=120)
    if mode:
        at.session_state["mode"] = mode
    return at.run()


def test_landing_page_shows_three_ways_to_start():
    at = new_app()
    assert not at.exception
    labels = [b.label for b in at.button]
    assert {"Open the example", "Choose a file", "Paste text"} <= set(labels)


def test_clicking_the_example_button_opens_the_example():
    at = new_app()
    next(b for b in at.button if b.label == "Open the example").click().run()
    assert not at.exception
    assert at.session_state["mode"] == "example"
    assert len(at.tabs) == 5


def test_example_runs_with_the_default_fasta():
    at = new_app("example")
    assert not at.exception
    assert not at.error


def test_example_runs_with_the_genbank_file():
    at = new_app("example")
    at.sidebar.selectbox[0].set_value("Synthetic demo with annotations (GenBank)").run()
    assert not at.exception


def test_back_button_returns_to_the_landing_page():
    at = new_app("example")
    next(b for b in at.sidebar.button if b.label == "Back to start").click().run()
    assert at.session_state["mode"] is None
    assert "Open the example" in [b.label for b in at.button]


def test_invalid_pasted_input_shows_an_error():
    at = new_app("paste")
    at.text_area[0].set_value("MKVLEEPQFHIL").run()
    assert not at.exception
    assert len(at.error) == 1


def test_a_pasted_sequence_is_analysed():
    at = new_app("paste")
    at.text_area[0].set_value("ATG" + "GCT" * 150 + "TAA").run()
    assert not at.exception
    assert not at.error
    assert len(at.tabs) == 5


def test_an_empty_upload_page_waits_for_a_file():
    at = new_app("upload")
    assert not at.exception
    assert not at.error
