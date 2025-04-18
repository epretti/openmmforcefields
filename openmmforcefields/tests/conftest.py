"""Default configuration and objects for tests"""

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--runespaloma",
        action="store_true",
        default=False,
        help="run espaloma tests",
    )

    parser.addoption(
        "--rungaff",
        action="store_true",
        default=False,
        help="run gaff tests",
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "espaloma: mark test as slow to run")
    config.addinivalue_line("markers", "gaff: mark test as using gaff")


def pytest_collection_modifyitems(config, items):
    skip_espaloma = pytest.mark.skip(reason="need --runespaloma option to run")

    if not config.getoption("--runespaloma"):
        for item in items:
            if "espaloma" in item.keywords:
                item.add_marker(skip_espaloma)

    skip_gaff = pytest.mark.skip(reason="need --rungaff option to run")

    if not config.getoption("--rungaff"):
        for item in items:
            if "gaff" in item.keywords:
                item.add_marker(skip_gaff)
