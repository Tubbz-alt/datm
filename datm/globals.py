"""
Script to hold variables meant to be shared throughout the repo.
"""
import logging
from pathlib import Path

from .utils import absolute_submodule_path

logger = logging.getLogger(__name__)

# Define some Path objects to folders within the repo
DIR_REPO = Path(absolute_submodule_path("datm/"))
DIR_DATA = DIR_REPO / "data"
DIR_DATA_RAW = DIR_DATA / "raw"
DIR_DATA_INT = DIR_DATA / "interim"
DIR_DATA_PROC = DIR_DATA / "processed"
DIR_LOGS = DIR_REPO / "logs"
DIR_NOTEBOOKS = DIR_REPO / "notebooks"
